CXX = g++
CXXFLAGS = -O3 -std=c++11


all: bin bin/construct_index bin/query_distance

bin:
	mkdir -p bin

bin/construct_index: construct_labelling_main.cpp highway_cover_labelling.h
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/query_distance: query_distance_main.cpp highway_cover_labelling.h
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^


clean:
	rm -rf bin
