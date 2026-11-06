CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -pedantic -I.

.PHONY: all clean

all: build/cmvf_barcodes

build/cmvf_barcodes: cmvf_barcodes.cpp $(wildcard mvf_module/*.hpp)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) cmvf_barcodes.cpp -o $@

clean:
	rm -rf build
