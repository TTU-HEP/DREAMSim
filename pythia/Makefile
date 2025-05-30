CXX = c++
CXXFLAGS = -O2 -std=c++17
ROOT_FLAGS := $(shell root-config --cflags --libs)

# PYTHIA8 paths
PYTHIA8_INCLUDE = $(PYTHIA8)/include
PYTHIA8_LIB = $(PYTHIA8)/lib
PYTHIA8_FLAGS = -I$(PYTHIA8_INCLUDE) -L$(PYTHIA8_LIB) -lpythia8

# Target
TARGET = test
SRC = partonShower.cc

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOT_FLAGS) $(PYTHIA8_FLAGS)

clean:
	rm -f $(TARGET)
