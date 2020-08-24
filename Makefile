CXXFLAGS = \
  -std=c++2a \
  -Wpedantic \
  -Wall \
  -Wextra \
  -Wconversion \
  -Wsign-conversion \
  $(OPTFLAGS) \
  $(INCLUDES) \

OPTFLAGS = \
  -O2 \
  -march=native \
  -mtune=native

INCLUDES = \
  -isystem./submodules/github.com/snsinfu/h5/include \
  -isystem./submodules/github.com/snsinfu/micromd/include

LIBS = \
  -lhdf5

PRODUCT = \
  main

OBJECTS = \
  src/main.o

ARTIFACTS = \
  $(OBJECTS) \
  $(PRODUCT)


.PHONY: all clean

all: $(PRODUCT)
	@:

clean:
	rm -f $(ARTIFACTS)

$(PRODUCT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)
