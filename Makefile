CC=g++
CFLAGS=-std=c++11 -Wall -Wextra -pedantic -lm -g

all: simulation

simulation: simulation.cpp
	$(CC) $(CFLAGS) simulation.cpp -o simulation

clean:
	rm simulation
