# Compiler
CXX = g++

# Compiler Flags
CXXFLAGS = -m64 -std=c++17 -Ofast -march=native -mtune=native \
           -funroll-loops -ftree-vectorize -fstrict-aliasing -fno-semantic-interposition \
           -fvect-cost-model=unlimited -fno-trapping-math -fipa-ra \
           -fassociative-math -ffast-math \
           -lgmp -lgmpxx -fopenmp -flto 

# Linker Flags
LDFLAGS = -lgmp -lgmpxx -fopenmp -flto

TARGET = kangaroo_solver
CXXFLAGS += -fPIC

# Source files
SRCS =  kangaroo.cpp 
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

# Targets
all: $(TARGET)
	@$(MAKE) clean_intermediates

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependencies
-include $(DEPS)

# Clean intermediate files only
clean_intermediates:
	rm -f $(OBJS) $(DEPS) && chmod +x $(TARGET)

# Clean everything )
clean: clean_intermediates
	rm -f $(TARGET)

.PHONY: all clean clean_intermediates
