# Parallel and Distributed Systems - Exercise 2: MPI Bitonic Sort

This repository contains an implementation of the Bitonic Sort algorithm designed for distributed execution using MPI (Message Passing Interface). The project leverages parallel processing to sort large datasets across multiple processes efficiently.

## Overview

The **Bitonic Sort** algorithm is a parallel sorting method that sorts data through a sequence of compare-and-swap operations. This implementation extends the algorithm to run on a distributed system using the **MPI** framework, enabling scalability across multiple processes. Key features include:

- Local sorting using **qsort**.
- Hypercube communication for inter-process data exchange.
- Performance benchmarking with built-in timers.

## Requirements

Before running the program, ensure you have the following installed:

- **MPI Implementation**: Compatible with OpenMPI or MPICH.
- **Linux/Unix Environment**: Required for running the shell commands.
- **Make**: For building the project.

## How to Compile and Run

### Compilation

To compile the program, use the provided `Makefile`. Ensure that the required source files (`main.c`, `bitonic_sort.c`, and `utilities.c`) and header files are in their respective directories.

```bash
make
```

### Running the Program

The program uses two parameters, `q` and `p`, to control the size of the dataset and the number of processes:

- `q`: Log2 of the number of elements to sort.
- `p`: Log2 of the number of processes.

#### Run Command

```bash
make run q=<value> p=<value>
```

- Replace `<value>` with the desired values for `q` and `p`.
- Example:

```bash
make run q=16 p=4
```

This command will sort \(2^{16}\) elements using \(2^4 = 16\) processes.

## Directory Structure

- **`src`**: Contains the source code (`main.c`, `bitonic_sort.c`, `utilities.c`).
- **`inc`**: Contains header files (`config.h`, `utilities.h`, `bitonic_sort.h`).
- **`Makefile`**: Used for building and cleaning the project.

## Key Functions

### Utilities

- **Argument Parsing**: Ensures valid `q` and `p` values.
- **Local Initialization**: Fills the array with random integers.
- **Verification**: Compares distributed sorting with a sequential reference.
- **Timing and Debugging**: Measures performance and prints debugging info.

### Bitonic Sort

- **Local Sorting**: Uses `qsort` for initial sorting.
- **Hypercube Exchange**: Performs compare-and-swap operations across processes.
- **Elbow Sort**: Merges cyclic bitonic sequences.

## Debugging and Verification

To enable debugging or sequential verification, uncomment the relevant lines in `config.h`:

```c
// #define DEBUG // Enables detailed output for debugging
// #define VERIFY_SEQUENTIAL // Enables comparison with sequential sorting
```

## Cleanup

To remove the compiled executable, use:

```bash
make clean
```

## Output

The program outputs the sorted data and timing statistics for each process. Use the debug mode to inspect intermediate states.

---
