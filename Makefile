# Specify compiler commands
CXX = mpic++ -O3 -Wall -fopenmp -std=c++11
# Specify name of file containing main function, without ext.
BINARY = Dockr
# Source files; ** wildcard does not work on my Make so just one level depth
SRCFILES = $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
# Filter out Test files
SRCFILES := $(filter-out $(wildcard src/*/*Test.cpp), $(filter-out $(wildcard src/*Test.cpp),$(SRCFILES)))
# Object files
# Change ext. and path
OBJFILES = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCFILES)))

# Testing
TEST = DockrTest
SRCTEST = $(wildcard src/*Test.cpp) $(wildcard src/*/*Test.cpp)
OBJTEST = $(patsubst src/%, obj/%, $(patsubst %.cpp,%.o,$(SRCTEST)))

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

# Test
test: testobjs
	$(CXX) $(OBJTEST) $(filter-out $(wildcard obj/Dockr.o),$(OBJFILES)) -o $(TEST) -lgtest -lgtest_main -lpthread
	./$(TEST) # --gtest_filter=Serialization.Pairs

testobjs: $(OBJTEST)

# Generate object file
obj/%.o: src/%.cpp
	$(CXX) -c -o $@ $<

clean:
	rm -rf $(BINARY)
	rm -rf obj
