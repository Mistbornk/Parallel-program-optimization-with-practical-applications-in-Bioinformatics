CXX = g++
CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -march=native -g

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

EXEC = main_test
EXECUTABLES = $(EXEC)

all: $(EXECUTABLES)

$(EXEC): $(OBJECTS)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXEC)

test: $(EXEC)
	@echo "Running test case:"
	./$(EXEC) seq1.fa seq2.fa

debug: $(BINDIR)/$(EXEC)
	cgdb ./$(BINDIR)/$(EXEC)
