#include <cstdint>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <h5.hpp>
#include <md.hpp>



struct simulation_config {
    md::index     chain_length;
    md::index     particle_count;
    md::scalar    core_diameter;
    md::scalar    core_repulsion;
    md::scalar    core_mobility;
    md::scalar    bond_length;
    md::scalar    bond_spring;
    md::scalar    flow_strength;
    md::scalar    box_size;
    md::scalar    temperature;
    md::scalar    timestep;
    md::step      simulation_length;
    md::step      sampling_interval;
    md::step      logging_interval;
    std::string   output_filename;
    std::uint64_t seed;
};

void run_simulation(simulation_config const& config);


int main()
{
    try {
        simulation_config const config = {
            .chain_length      = 100,
            .particle_count    = 900,
            .core_diameter     = 0.024,
            .core_repulsion    = 100,
            .core_mobility     = 1,
            .bond_length       = 0.02,
            .bond_spring       = 2.0e4,
            .flow_strength     = 1000,
            .box_size          = 0.46,
            .temperature       = 1,
            .timestep          = 1e-6,
            .simulation_length = 10000000,
            .sampling_interval = 10000,
            .logging_interval  = 10000,
            .output_filename   = "output.h5",
            .seed              = 0,
        };
        run_simulation(config);
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}


struct index_range {
    md::index start;
    md::index end;
};

template<>
struct h5::buffer_traits<std::vector<index_range>> {
    using buffer_type = std::vector<index_range>;
    using value_type = md::index;
    static constexpr int rank = 2;

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size(), 2};
    }

    static value_type* data(buffer_type& buffer)
    {
        return &buffer.data()->start;
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return &buffer.data()->start;
    }
};

template<>
struct h5::buffer_traits<md::array_view<md::point>> {
    using buffer_type = md::array_view<md::point>;
    using value_type = md::scalar;
    static constexpr int rank = 2;

    static h5::shape<rank> shape(buffer_type const& buffer)
    {
        return {buffer.size(), 3};
    }

    static value_type* data(buffer_type& buffer)
    {
        return &buffer.data()->x;
    }

    static value_type const* data(buffer_type const& buffer)
    {
        return &buffer.data()->x;
    }
};


class flow_forcefield : public md::forcefield {
public:
    flow_forcefield(std::vector<md::index> const& targets, md::scalar flow)
        : _targets{targets}, _flow{flow}
    {
    }

    md::scalar compute_energy(md::system const& system) override
    {
        (void) system;
        return 0;
    }

    void compute_force(md::system const& system, md::array_view<md::vector> forces) override
    {
        auto const positions = system.view_positions();

        for (auto const i : _targets) {
            auto dz = positions[i].z - 1;
            dz -= std::round(dz);

            if (std::fabs(dz) < 0.02) {
                forces[i].z += _flow;
            }
        }
    }

private:
    std::vector<md::index> _targets;
    md::scalar _flow;
};


void run_simulation(simulation_config const& config)
{
    h5::file store(config.output_filename, "w");

    std::mt19937_64 random(config.seed);
    md::system system;

    //
    // Particles
    //

    std::vector<index_range> chains;
    std::vector<md::index> particles;

    for (md::index i = 0; i < config.chain_length; i++) {
        system.add_particle({
            .mobility = config.core_mobility,
        });
    }
    chains.push_back(index_range{0, config.chain_length});

    for (md::index i = 0; i < config.particle_count; i++) {
        system.add_particle({
            .mobility = config.core_mobility,
        });
        particles.push_back(system.particle_count() - 1);
    }

    store.dataset<h5::i32, 2>("metadata/chain_ranges").write(chains);
    store.dataset<h5::i32, 1>("metadata/particles").write(particles);

    //
    // Forcefield
    //

    system.add_forcefield(
        md::make_neighbor_pairwise_forcefield<md::periodic_box>(
            md::softcore_potential<2, 3>{
                .energy   = config.core_repulsion,
                .diameter = config.core_diameter,
            }
        )
        .set_unit_cell({
            .x_period = config.box_size,
            .y_period = config.box_size,
            .z_period = config.box_size,
        })
        .set_neighbor_distance(config.core_diameter)
    );

    auto bonds = system.add_forcefield(
        md::make_bonded_pairwise_forcefield(
            md::spring_potential{
                .spring_constant      = config.bond_spring,
                .equilibrium_distance = config.bond_length,
            }
        )
    );

    for (auto const [start, end] : chains) {
        bonds->add_bonded_range(start, end);
    }

    system.add_forcefield(
        flow_forcefield(particles, config.flow_strength)
    );

    //
    // Initialization
    //

    for (auto const [start, end] : chains) {
        auto positions = system.view_positions();

        std::uniform_real_distribution<md::scalar> coord{0, config.box_size};
        std::normal_distribution<md::scalar> normal;

        md::point const center = {
            coord(random), coord(random), coord(random)
        };

        md::point pos;
        md::vector delta;
        for (md::index i = start; i < end; i++) {
            positions[i] = pos;
            delta += pos - center;

            md::vector step = {
                normal(random), normal(random), normal(random)
            };
            step *= config.bond_length / step.norm();

            pos += step;
        }
        delta /= md::scalar(end - start);

        for (md::index i = start; i < end; i++) {
            positions[i] -= delta;
        }
    }

    for (auto const i : particles) {
        auto positions = system.view_positions();

        std::uniform_real_distribution<md::scalar> coord{0, config.box_size};
        positions[i] = {
            coord(random), coord(random), coord(random)
        };
    }

    //
    // Simulation
    //

    auto callback = [&](md::step step) {
        if (step % config.logging_interval == 0) {
            auto const mean_energy =
                system.compute_energy() / md::scalar(system.particle_count());
            std::clog
                << step
                << '\t'
                << mean_energy
                << '\n';
        }

        if (step % config.sampling_interval == 0) {
            auto const step_key = std::to_string(step);
            std::vector<std::string> steps;
            if (auto steps_data = store.dataset<h5::str, 1>("snapshots/steps")) {
                steps_data.read_fit(steps);
            }
            steps.push_back(step_key);

            auto const positions = system.view_positions();
            store.dataset<h5::f32, 2>("snapshots/" + step_key + "/positions").write(
                positions,
                {.compression = 1, .scaleoffset = 4}
            );

            store.dataset<h5::str, 1>("snapshots/steps").write(steps);
        }
    };

    callback(0);

    md::simulate_brownian_dynamics(system, {
        .temperature = config.temperature,
        .timestep    = config.timestep,
        .steps       = config.simulation_length,
        .callback    = callback,
    });
}
