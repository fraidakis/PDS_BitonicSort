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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "mpi.h"

// #define DEBUG
// #define COMPARE_SECUENTIAL
#define VERIFY

struct timer
{
    double start;
    double end;
    double elapsed;
} local_sort_timer, elbow_sort_timer, communication_timer, total_timer, sequential_timer;

void debug_print(int *local_array, int local_size, int rank, int numProcesses, int dimension, int distance);

// Comparison function for qsort (stdlib sorting)
int compare_asc(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int compare_des(const void *a, const void *b)
{
    return (*(int *)b - *(int *)a);
}

void compare_and_swap(int *a, int *b, bool keep_min, int size)
{
    for (int i = 0; i < size; i++)
    {
        if ((a[i] > b[i]) == keep_min)
        {
            int temp = a[i];
            a[i] = b[i];
            b[i] = temp;
        }
    }
}

// Function to find the "elbow" of the bitonic sequence
// Returns the index of the maximum or minimum element
int find_elbow(int *arr, int size, bool *isAscendingFirst)
{
    int i = 0;
    while (i < size - 1 && arr[i] == arr[i + 1])
        i++;

    if (i != size - 1)
    {

        *isAscendingFirst = arr[i] < arr[i + 1]; // Determine the initial trend (ascending or descending)

        for (int i = 0; i < size - 1; i++)
            if ((*isAscendingFirst && arr[i] > arr[i + 1]) || (!*isAscendingFirst && arr[i] < arr[i + 1]))
                return i; // Return the index of the elbow
    }

    return -1; // If no elbow is found, sequence is entirely increasing or decreasing
}

// Function to merge a cyclic bitonic sequence using one elbow
void elbowSort(int *arr, int size, bool isAscending)
{
    elbow_sort_timer.start = MPI_Wtime();

    bool isAscendingFirst;
    int elbow = find_elbow(arr, size, &isAscendingFirst);

    if (elbow == -1)
    {
        // If no elbow is found, the array is already sorted
        if (isAscending != isAscendingFirst)
        {
            // If the array is not in the correct order, reverse it
            for (int i = 0; i < size / 2; i++)
            {
                int temp = arr[i];
                arr[i] = arr[size - i - 1];
                arr[size - i - 1] = temp;
            }
        }

        return;
    }

    int left = elbow;               // Start at the elbow (left of the transition)
    int right = (elbow + 1) % size; // Start at the first element after the elbow
    int sorted_index = 0;
    int *sorted = malloc(size * sizeof(int));

    // Merge the two parts
    while (sorted_index < size)
    {
        if (isAscendingFirst)
        {
            if (arr[left] >= arr[right])
            {
                sorted[sorted_index++] = arr[left];
                left = (left - 1 + size) % size; // Wrap around to the end if needed
            }
            else
            {
                sorted[sorted_index++] = arr[right];
                right = (right + 1) % size; // Wrap around to the start if needed
            }
        }

        else
        {
            if (arr[left] <= arr[right])
            {
                sorted[sorted_index++] = arr[left];
                left = (left - 1 + size) % size; // Wrap around to the end if needed
            }
            else
            {
                sorted[sorted_index++] = arr[right];
                right = (right + 1) % size; // Wrap around to the start if needed
            }
        }

        // Stop if left and right have fully wrapped around
        if ((left + 1) % size == right)
        {
            break;
        }
    }

    // Copy any remaining elements from one of the sides
    while (sorted_index < size)
    {
        sorted[sorted_index++] = arr[right];
        right = (right + 1) % size; // Wrap around if needed
    }

    // Copy the sorted array back to the original array
    if (isAscending == isAscendingFirst)
    {
        for (int i = 0; i < size; i++)
        {
            arr[i] = sorted[size - i - 1];
        }
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            arr[i] = sorted[i];
        }
    }

    free(sorted);

    elbow_sort_timer.end = MPI_Wtime();
    elbow_sort_timer.elapsed += elbow_sort_timer.end - elbow_sort_timer.start;
}

bool getSortingDirection(int rank, int dimension)
{
    return ((rank & (1 << dimension)) == 0); // 0 for ascending, 1 for descending
}

int getNumChunks(int local_size)
{
    int num_chunks = 1;
    while (local_size / num_chunks > (int)(1.5e6))
    {
        num_chunks *= 2;
    }

    return num_chunks;
}

void localSort(int *local_array, int local_size, int rank)
{
    local_sort_timer.start = MPI_Wtime(); // Start timer

    // Sort local array using qsort as initial local sorting
    if (rank % 2 != 0)
        qsort(local_array, local_size, sizeof(int), compare_des);
    else
        qsort(local_array, local_size, sizeof(int), compare_asc);

    local_sort_timer.end = MPI_Wtime(); // Stop timer
}

void hypercube_exchange_cmpswp(int *send_buffer, int *recv_buffer, int local_size, int partner, bool keep_min)
{
    communication_timer.start = MPI_Wtime();

    int num_chunks = getNumChunks(local_size);
    int chunk_size = local_size / num_chunks; // Divide local_size into smaller parts

    for (int chunk = 0; chunk < num_chunks; chunk++)
    {

        int offset = chunk * chunk_size;

        // Non-blocking send and receive for the current chunk
        MPI_Request send_request, recv_request;
        MPI_Irecv(recv_buffer + offset, chunk_size, MPI_INT, partner, 0, MPI_COMM_WORLD, &recv_request);
        MPI_Isend(send_buffer + offset, chunk_size, MPI_INT, partner, 0, MPI_COMM_WORLD, &send_request);

        // Perform computation on already received previous chunks (if any)
        if (chunk > 0)
            compare_and_swap(&send_buffer[(chunk - 1) * chunk_size], &recv_buffer[(chunk - 1) * chunk_size], keep_min, chunk_size);

        // Wait for the current chunk communication to complete
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
    }

    // Process the last chunk
    compare_and_swap(&send_buffer[(num_chunks - 1) * chunk_size], &recv_buffer[(num_chunks - 1) * chunk_size], keep_min, chunk_size);

    communication_timer.end = MPI_Wtime();
    communication_timer.elapsed += communication_timer.end - communication_timer.start;
}

void bitonicSort(int *send_buffer, int *recv_buffer, int local_size, int rank, int numProcesses, int p)
{
    int max_hypercube_dimension = p; // log2(numProcesses)

    // if (rank == 0)
    //     printf("\nNum chunks: %d\n", getNumChunks(local_size));

    for (int dimension = 1; dimension <= max_hypercube_dimension; dimension++)
    {

        bool isAscending = getSortingDirection(rank, dimension);

        for (int distance = 1 << (dimension - 1); distance > 0; distance >>= 1)
        {

#ifdef DEBUG
            debug_print(send_buffer, local_size, rank, numProcesses, dimension, distance);
#endif

            // Determine partner process
            int partner = rank ^ distance;                   // XOR operation to find partner process (bitwise complement)
            bool keep_min = (rank < partner) == isAscending; // If ASC lower rank keeps min, if DES higher rank keeps min (to form BITONIC sequence)

            hypercube_exchange_cmpswp(send_buffer, recv_buffer, local_size, partner, keep_min);
        }

        elbowSort(send_buffer, local_size, isAscending);

        // MPI_Barrier(MPI_COMM_WORLD);
    }
}

void sequential_sort_verify(int *send_buffer, int *recv_buffer, int *sequential_array, long int total_size, int local_size, int rank, int numProcesses)
{

#ifdef VERIFY
    // Gather results to process 0 for verification
    int *global_array = NULL;

    if (rank == 0)
    {
        global_array = malloc(total_size * sizeof(int));
    }

    MPI_Gather(send_buffer, local_size, MPI_INT, global_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    // Verify sorting on process 0
    if (rank == 0)
    {

        sequential_timer.start = MPI_Wtime(); // Start timer

        // Sort the verification copy using standard qsort
        qsort(sequential_array, total_size, sizeof(int), compare_asc);

        sequential_timer.end = MPI_Wtime(); // Stop timer

        sequential_timer.elapsed = sequential_timer.end - sequential_timer.start; // Calculate elapsed time
        printf("Sequential time : %f\n", sequential_timer.elapsed);

#ifdef VERIFY
        // Compare sorted results
        int is_sorted = 1;
        for (int i = 0; i < total_size; i++)
        {
            if (global_array[i] != sequential_array[i])
            {
                printf("Error: Mismatch at index %d: %d != %d\n", i, global_array[i], sequential_array[i]);
                is_sorted = 0;
                break;
            }
        }

        printf("\n%s sorting %ld elements across %d processes (%d elements each)\n\n\n", is_sorted ? "SUCCESSFUL" : "FAILED", total_size, numProcesses, local_size);

        printf("Speedup: %f   ||   Efficiency: %f\n\n", sequential_timer.elapsed / total_timer.elapsed, (sequential_timer.elapsed / total_timer.elapsed) / numProcesses);
        free(global_array);
#endif

        free(sequential_array);
    }
}

void argument_check(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <q> <p>\n", argv[0]);
        fprintf(stderr, "    q: log2 of local array size\n");
        fprintf(stderr, "    p: log2 of number of processes\n");
        exit(1);
    }
}

void validate_parameters(int rank, int numProcesses, int p)
{
    if (numProcesses != (1 << p))
    {
        if (rank == 0)
        {
            fprintf(stderr, "Error: Number of processes must be 2^p\n");
        }
        MPI_Finalize();
        exit(1);
    }
}

void initialize_array(int *local_array, int local_size, int rank, int numProcesses)
{
    // Seed random number generator differently for each process
    srand(time(NULL) + rank); // Seed based on current time and rank

    // Generate random integers for local array
    for (int i = 0; i < local_size; i++)
    {
        local_array[i] = rand() % 1000; // Random integers between 0 and 9999
    }
}

int main(int argc, char **argv)
{

    // Verify for correct number of arguments
    argument_check(argc, argv);

    // Parse command-line arguments
    int q = atoi(argv[1]); // log2 of local array size
    int p = atoi(argv[2]); // log2 of number of processes

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank and total number of processes
    int rank, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    // Validate input parameters
    validate_parameters(rank, numProcesses, p);

    // Calculate local array size
    int local_size = 1 << q; // 2^q elements per process
    long int total_size = (long int)local_size * numProcesses;

    // Allocate memory for local array
    int *send_buffer = malloc(local_size * sizeof(int));
    int *recv_buffer = malloc(local_size * sizeof(int));

    initialize_array(send_buffer, local_size, rank, numProcesses);


#ifdef VERIFY
#define COMPARE_SECUENTIAL
#endif

    #ifdef COMPARE_SECUENTIAL
        int *sequenteal_array = NULL;

    if (rank == 0)
    {
        sequenteal_array = malloc(numProcesses * local_size * sizeof(int));
    }

    MPI_Gather(send_buffer, local_size, MPI_INT, sequenteal_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    total_timer.start = MPI_Wtime(); // Start timer

    localSort(send_buffer, local_size, rank);

    bitonicSort(send_buffer, recv_buffer, local_size, rank, numProcesses, p);

    total_timer.end = MPI_Wtime(); // Stop timer

    total_timer.elapsed = total_timer.end - total_timer.start;                // Calculate elapsed time
    local_sort_timer.elapsed = local_sort_timer.end - local_sort_timer.start; // Calculate elapsed time

#ifdef DEBUG
    debug_print(send_buffer, local_size, rank, numProcesses, 0, 0);
#endif

    double global_communication, global_initialSort, global_merge, global_all;
    MPI_Reduce(&communication_timer.elapsed, &global_communication, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elbow_sort_timer.elapsed, &global_merge, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sort_timer.elapsed, &global_initialSort, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_timer.elapsed, &global_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    global_initialSort /= numProcesses;
    global_merge /= numProcesses;
    global_communication /= numProcesses;
    global_all /= numProcesses;

    if (rank == 0)
    {
        printf("Sort: %f seconds\n", global_initialSort);
        printf("Merge: %f seconds\n", global_merge);
        printf("Communication: %f seconds\n", global_communication);
        printf("Total: %f seconds\n", global_all);
    }

#ifdef COMPARE_SECUENTIAL

    sequential_sort_verify(send_buffer, recv_buffer, sequenteal_array, total_size, local_size, rank, numProcesses);

#endif

    // Clean up
    free(send_buffer);
    free(recv_buffer);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

void debug_print(int *local_array, int local_size, int rank, int numProcesses, int dimension, int distance)
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Gather and print the global array for debugging (rank 0 only)
    int *global_array = NULL;

    if (rank == 0)
    {
        global_array = malloc(numProcesses * local_size * sizeof(int));
    }

    MPI_Gather(local_array, local_size, MPI_INT, global_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

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
        free(global_array);

        if (distance == 1 << dimension - 1)
        {
            printf("\n\n\n\n");

            for (int rank = 0; rank < numProcesses; rank++)
            {
                int partner = rank ^ distance; // XOR operation to find partner process (bitwise complement)
                bool direction = getSortingDirection(rank, dimension);

                // bool keep_min = (rank < partner) == direction;

                printf("\nDir of proccess %d: %s", rank, direction ? "ASC" : "DES");
            }
        }

        printf("\n\nAfter: Dimenison %d, Distance %d\n", dimension, distance);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
