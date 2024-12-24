// ! Hints on : 5th videolecture , minute 1:25:00 - 1:33:00  && 1:47:00 - 1:54:00

/**
 * Distributed Bitonic Sort Implementation
 *
 * This program demonstrates a parallel implementation of Bitonic sort
 * using MPI (Message Passing Interface) for distributed computing.
 *
 * Key Concepts:
 * - Uses 2^p processes
 * - Each process handles 2^q random integers
 * - Total number of elements: N = 2^(p+q)
 * - Implements the Bitonic sorting network algorithm
 */

//! bash run.sh 25 2
// Runs for 2^25 elements and 2^2 processes

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

double get_time_in_seconds(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

// Comparison function for qsort (stdlib sorting)
int compare_asc(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int compare_des(const void *a, const void *b)
{
    return (*(int *)b - *(int *)a);
}

/**
 * Swap elements if they are in the wrong order based on the sorting direction
 * @param a Pointer to first element
 * @param b Pointer to second element
 * @param dir Sorting direction (TRUE for ascending, FALSE for descending)
 */
void compare_and_swap(int *a, int *b, bool dir)
{
    if ((dir && *a > *b) || (!dir && *a < *b))
    {
        int temp = *a;
        *a = *b;
        *b = temp;
    }
}

/**
 * Sort of a bitonic sequence
 * @param arr Array to merge
 * @param start Starting index
 * @param size Number of elements to merge
 * @param dir Sorting direction
 */
void bitonic_merge(int *arr, int start, int size, bool dir)
{
    if (size > 1)
    {
        int mid = size / 2;

        // Perform bitonic split
        for (int i = start; i < start + mid; i++)
        {
            compare_and_swap(&arr[i], &arr[i + mid], dir);
        }

        // Recursively merge the two halves
        bitonic_merge(arr, start, mid, dir);
        bitonic_merge(arr, start + mid, mid, dir);
    }
}

void bitonic_merge_partial(int *arr, int start, int size, bool dir)
{
    int mid = size / 2;

    // Perform bitonic split
    for (int i = start; i < start + mid; i++)
    {
        compare_and_swap(&arr[i], &arr[i + mid], dir);
    }
}

void debug_print(int *local_array, int local_size, int rank, int numProcesses, int size, int step)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Gather and print the global array for debugging (rank 0 only)
    int *global_array = NULL;

    if (rank == 0)
    {
        global_array = malloc(numProcesses * local_size * sizeof(int));
    }

    MPI_Gather(local_array, local_size, MPI_INT,
               global_array, local_size, MPI_INT,
               0, MPI_COMM_WORLD);

    if (rank == 0)
    {

        int processor = 0;
        for (int i = 0; i < numProcesses * local_size; i++)
        {
            if (i % local_size == 0)
            {
                printf("\nProcess %d: ", processor);
                processor++;
            }
            printf("%d ", global_array[i]);
        }
        printf("\n\n");
        free(global_array);

        for (int i = 0; i < numProcesses; i++)
        {
            printf("\nDir of proccess %d: %s", i, (((i >> (int)(log2(size) + 1)) & 1) == ((i >> (int)log2(size)) & 1)) ? "ASC" : "DES");
        }
    }

    if (rank == 0)
    {
        printf("\n\nAfter: Size %d, Step %d\n", size, step);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * Get the direction of the bitonic sequence based on the rank and dimension
 * !Check if the bit at the dimension position is the same as the bit at the dimension + 1 position
 * @param rank Process rank
 * @param dimension Hypercube dimension
 * @return TRUE for ascending, FALSE for descending
 */
bool getSortingDirection(int rank, int dimension)
{
    return (((rank >> (dimension + 1)) & 1) == ((rank >> dimension) & 1));
}


int main(int argc, char **argv)
{

    // Check for correct number of arguments
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <q> <p>\n", argv[0]);
        return 1;
    }

    // Parse command-line arguments
    int q = atoi(argv[1]); // log2 of local array size
    int p = atoi(argv[2]); // log2 of number of processes

    double start, end; // Variables to store execution time

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank and total number of processes
    int rank, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    // Validate input parameters
    if (numProcesses != (1 << p))
    {
        if (rank == 0)
        {
            fprintf(stderr, "Error: Number of processes must be 2^p. Size = %d", numProcesses);
        }
        MPI_Finalize();
        return 1;
    }

    // Calculate local array size
    int local_size = 1 << q; // 2^q elements per process
    int total_size = local_size * numProcesses;

    // Allocate memory for local array
    int *local_array = malloc(local_size * sizeof(int));
    int *recv_buffer = malloc(local_size * sizeof(int));

    // Seed random number generator differently for each process
    srand(time(NULL) + rank); // Seed based on current time and rank

    // Generate random integers for local array
    for (int i = 0; i < local_size; i++)
    {
        local_array[i] = rand() % 10000; // Random integers between 0 and 9999
    }

    int *sequenteal_array = NULL;

    if (rank == 0)
    {
        sequenteal_array = malloc(numProcesses * local_size * sizeof(int));
    }

    MPI_Gather(local_array, local_size, MPI_INT,
               sequenteal_array, local_size, MPI_INT,
               0, MPI_COMM_WORLD);

    start = get_time_in_seconds(); // Start timer

    // Sort local array using qsort as initial local sorting
    if (rank % 2 != 0)
        qsort(local_array, local_size, sizeof(int), compare_des);
    else
        qsort(local_array, local_size, sizeof(int), compare_asc);

    MPI_Barrier(MPI_COMM_WORLD);



    int max_hypercube_dimension = log2(numProcesses);

    for (int dimension = 1; dimension <= max_hypercube_dimension; dimension++)
    { 
        for (int dist = 1 << (dimension - 1); dist > 0; dist >>= 1)
        { 

            // debug_print(local_array, local_size, rank, numProcesses, 1 << dimension, dist);

            // Determine partner process
            int partner = rank ^ dist; // XOR operation to find partner process (bitwise complement)

            // Determine sorting direction (ascending or descending)
            bool dir = getSortingDirection(rank, dimension);
            bool partenr_dir = getSortingDirection(partner, dimension);

            // Exchange data with partner process
            MPI_Sendrecv(local_array, local_size, MPI_INT, partner, 0,
                         recv_buffer, local_size, MPI_INT, partner, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Perform bitonic comparison and merge
            int *merged_array = malloc(2 * local_size * sizeof(int));

            // Copy local and received data to merged array
            memcpy(merged_array, local_array, local_size * sizeof(int));
            memcpy(merged_array + local_size, recv_buffer, local_size * sizeof(int));

            if (dist == 1)
                bitonic_merge(merged_array, 0, 2 * local_size, dir);
            else
                bitonic_merge_partial(merged_array, 0, 2 * local_size, dir);

            // Update local array with merged results

            if (dir || dir == partenr_dir)
            {
                if (rank < partner)
                    memcpy(local_array, merged_array, local_size * sizeof(int));
                else
                    memcpy(local_array, merged_array + local_size, local_size * sizeof(int));
            }
            else
            {
                if (rank > partner)
                    memcpy(local_array, merged_array, local_size * sizeof(int));
                else
                    memcpy(local_array, merged_array + local_size, local_size * sizeof(int));
            }

            free(merged_array);

            // Debugging barrier to synchronize processes
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if (rank == numProcesses - 1)
    {
        end = get_time_in_seconds();       // Stop timer
        double elapsed_time = end - start; // Calculate elapsed time
        printf("MPI execution time : %f second\n", elapsed_time);
    }

    // debug_print(local_array, local_size, rank, numProcesses, 0, 0);

    if (rank == 0)
    {
        start = get_time_in_seconds(); // Start timer

        qsort(sequenteal_array, numProcesses * local_size, sizeof(int), compare_asc);

        end = get_time_in_seconds();       // Stop timer
        double elapsed_time = end - start; // Calculate elapsed time
        printf("Sequential execution time : %f second\n", elapsed_time);
    }

    // Gather results to process 0 for verification
    int *global_array = NULL;
    if (rank == 0)
    {
        global_array = malloc(total_size * sizeof(int));
    }

    MPI_Gather(local_array, local_size, MPI_INT,
               global_array, local_size, MPI_INT,
               0, MPI_COMM_WORLD);

    // Verify sorting on process 0
    if (rank == 0)
    {

        // Create a copy for verification
        int *verify_array = malloc(total_size * sizeof(int));
        memcpy(verify_array, global_array, total_size * sizeof(int));

        // Sort the verification copy using standard qsort
        qsort(verify_array, total_size, sizeof(int), compare_asc);

        // Compare sorted results
        int is_sorted = 1;
        for (int i = 0; i < total_size; i++)
        {
            if (global_array[i] != verify_array[i])
            {
                is_sorted = 0;
                break;
            }
        }

        printf("Sorting %d elements across %d processes %s\n",
               total_size, numProcesses, is_sorted ? "SUCCESSFUL" : "FAILED");

        free(verify_array);
        free(global_array);
    }

    // Clean up
    free(local_array);
    free(recv_buffer);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}