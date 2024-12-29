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

#include "bitonic_sort.h"

int main(int argc, char **argv)
{
    int q, p;

    // Parse command-line arguments
    if (!parse_arguments(argc, argv, &q, &p))
    {
        return 1;
    }

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank and total number of processes
    int rank, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    // Validate input parameters
    if (!validate_parameters(numProcesses, p, rank))
    {
        MPI_Finalize();
        return 1;
    }

    // Calculate local array size
    int local_size = 1 << q; // 2^q elements per process
    long int total_size = (long int)local_size * numProcesses;

    // Allocate memory for local array
    int *send_buffer = malloc(local_size * sizeof(int));
    int *recv_buffer = malloc(local_size * sizeof(int));

    // Initialize local array with random data (each process generates its own data)
    initialize_data(send_buffer, local_size, rank);

#ifdef VERIFY_SEQUENTIAL
    // Initialize sequential array which gathers all data to process 0 
    int *sequential_array = initialize_sequential_array(rank, numProcesses, local_size, send_buffer);
#endif


    // ********** Start of the sorting process **********

    start_timer(&total_timer); 

    // Step 1: Local sort
    localSort(send_buffer, local_size, rank);

    // Step 2: Bitonic sort (communication + elbow sort)
    bitonicSort(send_buffer, recv_buffer, local_size, rank, p);

    stop_timer(&total_timer);

    // ********** End of the sorting process **********


    // Gather and print global statistics
    gather_and_print_statistics(rank, numProcesses, total_size);

#ifdef VERIFY_SEQUENTIAL
    // Verify sorting (SUCCESSFUL or FAILED)
    verify_sequential_sort(send_buffer, sequential_array, total_size, local_size, rank, numProcesses);
#endif

    // Clean up
    free(send_buffer);
    free(recv_buffer);

    MPI_Finalize();

    return 0;
}
