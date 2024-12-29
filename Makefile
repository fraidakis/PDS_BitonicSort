# Makefile for MPI Bitonic Sort

# Compiler and flags
CC = mpicc
CFLAGS = -Wall -O3

# Executable name
EXEC = bitonic_sort

# Source and include directories
SRC_DIR = src
INC_DIR = inc

# Source files
SRCS = $(SRC_DIR)/main.c $(SRC_DIR)/bitonic_sort.c $(SRC_DIR)/utilities.c

# Header files (dependencies)
HEADERS = $(INC_DIR)/config.h

# Default target
all: $(EXEC)

# Compile and link in a single step
$(EXEC): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $^ -o $@

# Run the program with arguments q and p
run: $(EXEC)
	@if [ "$(q)" = "" ] || [ "$(p)" = "" ]; then \
		echo "Usage: make run q=<value> p=<value>"; \
		exit 1; \
	fi; \
	procs=$$(echo "2^$(p)" | bc); \
	mpiexec --oversubscribe -n $${procs} ./$(EXEC) $(q) $(p)

# Clean up build files
clean:
	rm -f $(EXEC)

# Phony targets
.PHONY: all clean run
