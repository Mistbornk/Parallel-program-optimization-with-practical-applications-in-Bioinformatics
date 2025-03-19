CXX = g++
CXXFLAGS = -std=c++20 -Wall -O3 -I./inc -g

SRCDIR = src
OBJDIR = obj

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

EXEC_MATRIX = matrix_test
EXEC_THREADPOOL = threadpool_test
EXECUTABLES = $(EXEC_MATRIX) $(EXEC_THREADPOOL)

all: $(EXECUTABLES)

$(EXEC_MATRIX): $(OBJDIR)/matrix_test.o
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXEC_THREADPOOL): $(OBJDIR)/threadpool_test.o $(OBJDIR)/threadPool.o
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) $(EXECUTABLES)

debug: $(EXEC_MATRIX)
	cgdb ./$(EXEC_MATRIX)
