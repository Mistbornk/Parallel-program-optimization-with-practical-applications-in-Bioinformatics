CXX = g++
CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -g

SRCDIR = src
OBJDIR = obj

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

EXEC = FM_Index
EXECUTABLES = $(EXEC)

all: $(EXECUTABLES)

$(EXEC): $(OBJDIR)/main.o $(OBJDIR)/FM_Index.o
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXECUTABLES)

debug: $(EXEC)
	cgdb ./$(EXEC_MATRIX)