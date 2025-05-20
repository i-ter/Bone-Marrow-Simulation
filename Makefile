# Compiler
CXX = clang++

# Compiler Flags
CXXFLAGS = -std=c++17 -I /opt/homebrew/include -Wall -Wextra -g

# Linker Flags
LDFLAGS = -L /opt/homebrew/lib -lsfml-graphics -lsfml-window -lsfml-system

# OpenMP Flags
OMP_INC := $(shell brew --prefix libomp)/include
OMP_LIB := $(shell brew --prefix libomp)/lib
OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(OMP_INC)
OMP_LDFLAGS = -L$(OMP_LIB) -lomp

SIM = main
MOVIE_GEN = gen_movie_main

# Source Files
SRCS = main.cpp
MOVIE_SRCS = gen_movie_main.cpp movie_generator.cpp

# Build all targets
all: $(SIM) $(MOVIE_GEN)

# Main simulation target
$(SIM): $(SRCS) cell_config.h
	$(CXX) $(CXXFLAGS) $(OMP_CFLAGS) $(SRCS) -o $(SIM) $(LDFLAGS) $(OMP_LDFLAGS)

# Movie generator target
$(MOVIE_GEN): $(MOVIE_SRCS)
	$(CXX) $(CXXFLAGS) $(OMP_CFLAGS) $(MOVIE_SRCS) -o $(MOVIE_GEN) $(LDFLAGS) $(OMP_LDFLAGS)

# Clean Rule (Removes the compiled files)
clean:
	rm -f $(SIM) $(MOVIE_GEN)

