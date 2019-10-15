# Specify compiler commands
CXX = mpic++ -O3 -Wall -fopenmp -std=c++11
# Specify name of file containing main function, without ext.
BINARY = finDrVS
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))


$(info    Source files: $(SRCFILES))
$(info    Object files: $(OBJFILES))
all: folders compile

# Create required folders
folders:
	mkdir -p obj
	mkdir -p obj/Serialization
	mkdir -p obj/VinaInstance

# Link everything together 
compile: objs
	$(CXX) $(OBJFILES) -o $(BINARY) -lm

# Generate all object files
objs: $(OBJFILES)

testobjs: $(OBJTEST)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm -rf $(BINARY)
	rm -rf obj
