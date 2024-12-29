// Header File: utilities.h
#ifndef UTILITIES_H
#define UTILITIES_H

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

extern struct timer communication_timer, elbow_sort_timer, local_sort_timer, total_timer;

/**
 * Parses command-line arguments.
 * @param argc - Number of arguments.
 * @param argv - Argument array.
 * @param q - Pointer to store the value of q.
 * @param p - Pointer to store the value of p.
 * @return True if arguments are valid, false otherwise.
 */
bool parse_arguments(int argc, char **argv, int *q, int *p);

/**
 * Validates MPI parameters and ensures correct process count.
 * @param numProcesses - Total number of MPI processes running.
 * @param p - Log2 of the number of processes user wants to run.
 * @param rank - Rank of the current process.
 * @return True if parameters are valid, false otherwise.
 */
bool validate_parameters(int numProcesses, int p, int rank);

/**
 * Initializes local array with random integers.
 * @param array - Pointer to the array.
 * @param size - Size of the array.
 * @param rank - Rank of the current process (used for seeding).
 */
void initialize_data(int *localArray, int size, int rank);

/**
 * Initializes a sequential array and gathers all data
 * @param rank - Rank of the current process.
 * @param numProcesses - Total number of MPI processes.
 * @param local_size - Size of the local array.
 * @param send_buffer - The local array.
 * @return Pointer to the sequential array.
 */
int* initialize_sequential_array(int rank, int numProcesses, int local_size, int *send_buffer);

/**
 * Performs local sorting of an array using qsort.
 * @param local_array - The array to sort.
 * @param local_size - Number of elements in the array.
 * @param rank - Rank of the current MPI process (used to determine sorting order).
 */
void localSort(int *local_array, int local_size, int rank);

/**
 * Calculates number of chunks for communication based on local size.
 * @param local_size - Number of elements in the local array.
 * @return Number of chunks.
 */
int getNumChunks(int local_size);

/**
 * Determines the sorting direction of process based on rank and dimension.
 * @param rank - Rank of the current MPI process.
 * @param dimension - Current dimension of the hypercube network.
 * @return True for descending, false for ascending.
 */
bool getSortingDirection(int rank, int dimension);

/**
 * Compares and swaps two arrays based on the sorting order.
 * @param a - First array.
 * @param b - Second array.
 * @param keep_min - If true, keeps the minimum values in 'a'; otherwise, keeps the maximum.
 * @param size - Number of elements in each array.
 */
void compare_and_swap(int *a, int *b, bool keep_min, int size);

/**
 * Gathers and prints timing statistics for all processes.
 * @param rank - Rank of the current process.
 * @param numProcesses - Total number of MPI processes.
 * @param total_size - Total number of elements across all processes.
 */
void gather_and_print_statistics(int rank, int numProcesses, long total_size);

/**
 * Verifies the sorting against a sequential reference.
 * @param send_buffer - Local array after sorting.
 * @param sequential_array - Sequentially sorted reference array.
 * @param total_size - Total number of elements across all processes.
 * @param local_size - Number of elements in the local array.
 * @param rank - Rank of the current process.
 * @param numProcesses - Total number of MPI processes.
 */
void verify_sequential_sort(int *send_buffer, int *sequential_array, long total_size, int local_size, int rank, int numProcesses);

/**
 * Prints the local array of each process for debugging purposes.
 * @param local_array - Local array 
 * @param local_size - Number of elements in the local array.
 * @param rank - Rank of the current MPI process.
 * @param numProcesses - Total number of MPI processes.
 * @param dimension - Current dimension of the hypercube network.
 * @param distance - Distance to the partner process in the current step.
 */
void debug_print(int *local_array, int local_size, int rank, int numProcesses, int dimension, int distance);

/**
 * Compares integers in ascending order (used for qsort).
 * @param a - First integer pointer.
 * @param b - Second integer pointer.
 * @return Negative if *a < *b, zero if *a == *b, positive if *a > *b.
 */
int compareAscending(const void *a, const void *b);

/**
 * Compares integers in descending order (used for qsort).
 * @param a - First integer pointer.
 * @param b - Second integer pointer.
 * @return Negative if *a > *b, zero if *a == *b, positive if *a < *b.
 */
int compareDescending(const void *a, const void *b);

/**
 * Structure for timing operations.
 */
struct timer
{
    double start;
    double end;
    double elapsed;
};

/**
 * Initializes a timer by setting its start time.
 * @param t - Pointer to the timer structure.
 */
void start_timer(struct timer *t);

/**
 * Stops the timer and calculates elapsed time.
 * @param t - Pointer to the timer structure.
 */
void stop_timer(struct timer *t);

#endif // UTILITIES_H
