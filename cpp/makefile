# CellSpace v0.0.5

CXX = g++
CXXFLAGS = -pthread -std=gnu++11

BOOST_DIR = /usr/local/bin/boost_1_63_0/

OBJS = normalize.o dict.o args.o proj.o parser.o data.o model.o starspace.o utils.o
INCLUDES = -I$(BOOST_DIR)

opt: CXXFLAGS += -O3 -funroll-loops
opt: starspace

normalize.o: src/utils/normalize.cpp src/utils/normalize.h
	$(CXX) $(CXXFLAGS) -g -c src/utils/normalize.cpp

dict.o: src/dict.cpp src/dict.h src/utils/args.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/dict.cpp

args.o: src/utils/args.cpp src/utils/args.h
	$(CXX) $(CXXFLAGS) -g -c src/utils/args.cpp

model.o: data.o src/model.cpp src/model.h src/utils/args.h src/proj.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/model.cpp

proj.o: src/proj.cpp src/proj.h src/matrix.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/proj.cpp

data.o: parser.o src/data.cpp src/data.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/data.cpp -o data.o

utils.o: src/utils/utils.cpp src/utils/utils.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/utils/utils.cpp -o utils.o

parser.o: dict.o src/parser.cpp src/parser.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/parser.cpp -o parser.o

starspace.o: src/starspace.cpp src/starspace.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -g -c src/starspace.cpp

starspace: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(INCLUDES) -g src/main.cpp -o CellSpace

clean:
	rm -rf *.o CellSpace CellSpace.dSYM
