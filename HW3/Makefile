CXX = g++
NVCC = nvcc

CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -march=native -g
NVCCFLAGS = -std=c++20 -O3 -arch=sm_60 -I./inc -diag-suppress=2189

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES_CPP = $(wildcard $(SRCDIR)/*.cpp)
SOURCES_CU  = $(wildcard $(SRCDIR)/*.cu)
OBJECTS_CPP = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES_CPP))
OBJECTS_CU  = $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(SOURCES_CU))

EXEC = main_test
EXECUTABLES = $(EXEC)

all: $(EXECUTABLES)

$(EXEC): $(OBJECTS_CPP) $(OBJECTS_CU)
	@echo "Linking $@..."
	$(NVCC) -o $@ $^ -lcudart

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling C++ $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cu | $(OBJDIR)
	@echo "Compiling CUDA $<..."
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXEC)

test: $(EXEC)
	@echo "Running test case:"
	./$(EXEC) seq1.fa seq2.fa

debug: $(EXEC)
	cgdb ./$(EXEC)
