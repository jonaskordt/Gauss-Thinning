CC = g++
CCFLAGS = -std=c++11 -O3 -Wno-format -shared -fPIC

INC_DIRS = /usr/include /usr/include/eigen3 ../libigl/include
INC_FLAGS = -L/usr/lib $(addprefix -I,$(INC_DIRS)) $(shell python3 -m pybind11 --includes) 
EXAMPLES = architecture boat coffee cone cylinder  dome drill face fandisk  guitar lilium mask nut swing dog fertility bumpy bunny cone_high

mainParallel: mainParallel.cpp
	$(CC) $(CCFLAGS) -fopenmp $(INC_FLAGS) mainParallel.cpp -o mainParallel$(shell python3-config --extension-suffix) 

example: example.cpp
	$(CC) $(CCFLAGS) -fopenmp $(INC_FLAGS) example.cpp -o example$(shell python3-config --extension-suffix) 