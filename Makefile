CXX = g++
CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -g

SRCDIR = src
OBJDIR = obj

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXEC = matrix

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@echo "Linking $(EXEC)..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXEC)

debug: $(EXEC)
	cgdb ./$(EXEC)
