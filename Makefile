# Compiler
CXX = g++

# Compiler Flags
CXXFLAGS = -std=c++17 -I /opt/homebrew/include -Wall -Wextra -g

# Linker Flags
LDFLAGS = -L /opt/homebrew/lib -lsfml-graphics -lsfml-window -lsfml-system

SIM = main
MOVIE_GEN = gen_movie

# Source Files
SRCS = main.cpp
MOVIE_SRCS = gen_movie_main.cpp MovieGenerator.cpp

# Build all targets
all: $(SIM) $(MOVIE_GEN)

# Main simulation target
$(SIM): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(SIM) $(LDFLAGS)

# Movie generator target
$(MOVIE_GEN): $(MOVIE_SRCS)
	$(CXX) $(CXXFLAGS) $(MOVIE_SRCS) -o $(MOVIE_GEN) $(LDFLAGS)

# Clean Rule (Removes the compiled files)
clean:
	rm -f $(SIM) $(MOVIE_GEN)

