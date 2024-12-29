// Implementation File: bitonic_sort.c

#include "bitonic_sort.h"

void bitonicSort(int *send_buffer, int *recv_buffer, int local_size, int rank, int p)
{
    // Number of dimensions in the hypercube - log2(numProcesses)
    int max_hypercube_dimension = p;

    // Iterate through each dimension of the hypercube.
    for (int dimension = 1; dimension <= max_hypercube_dimension; dimension++)
    {
        // Determine the sorting direction (ascending or descending) for this dimension.
        bool isAscending = getSortingDirection(rank, dimension);

        // Iterate through the distances within this dimension (powers of 2) starting from the farthest.
        for (int distance = 1 << (dimension - 1); distance > 0; distance >>= 1)
        {

#ifdef DEBUG
            debug_print(send_buffer, local_size, rank, 1 << p, dimension, distance);
#endif

            // Determine the partner process for communication using XOR operation (bitwise complement).
            int partner = rank ^ distance;

            // Determine whether to keep the minimum or maximum of the compared elements.
            bool keep_min = (rank < partner) == isAscending; // If ascending and rank < partner, keep min, otherwise keep max.

            // Perform the communication and compare-and-swap operation with the partner.
            hypercube_exchange_cmpswp(send_buffer, recv_buffer, local_size, partner, keep_min);
        }

        // Perform elbow sort to locally sort the bitonic sequence after communication within the dimension.
        elbowSort(send_buffer, local_size, isAscending);
    }
}

void hypercube_exchange_cmpswp(int *send_buffer, int *recv_buffer, int local_size, int partner, bool keep_min)
{
    communication_timer.start = MPI_Wtime();

    int num_chunks = getNumChunks(local_size); // Get the number of chunks to divide the local array.
    int chunk_size = local_size / num_chunks;  // Size of each chunk.

    // Iterate through each chunk.
    for (int chunk = 0; chunk < num_chunks; chunk++)
    {
        int offset = chunk * chunk_size; // Calculate the offset for the current chunk.

        // Non-blocking send and receive for the current chunk.
        MPI_Request send_request, recv_request;
        MPI_Irecv(recv_buffer + offset, chunk_size, MPI_INT, partner, chunk, MPI_COMM_WORLD, &recv_request); // Receive from partner
        MPI_Isend(send_buffer + offset, chunk_size, MPI_INT, partner, chunk, MPI_COMM_WORLD, &send_request); // Send to partner

        // Perform computation (compare_and_swap) on previously received chunks (pipelining).
        if (chunk > 0)
            compare_and_swap(&send_buffer[(chunk - 1) * chunk_size], &recv_buffer[(chunk - 1) * chunk_size], keep_min, chunk_size);

        // Wait for the current chunk communication to complete.
        MPI_Wait(&send_request, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
    }

    // Process the last chunk (compare and swap) after communication is complete.
    compare_and_swap(&send_buffer[(num_chunks - 1) * chunk_size], &recv_buffer[(num_chunks - 1) * chunk_size], keep_min, chunk_size);

    communication_timer.end = MPI_Wtime();
}

int find_elbow(int *arr, int size, bool *isAscendingFirst)
{
    int i = 0;

    // Skip over any initial sequence of equal elements. This is important to determine the initial trend.
    while (i < size - 1 && arr[i] == arr[i + 1])
        i++;

    // If we haven't reached the end of the array, it means there are at least two distinct elements.
    if (i != size - 1)
    {
        // Determine the initial trend (ascending or descending) based on the first two distinct elements.
        *isAscendingFirst = arr[i] < arr[i + 1];

        // Iterate through the array to find the "elbow" (the point where the trend changes).
        for (int i = 0; i < size - 1; i++)
            // Check if the current pair of elements violates the initial trend.
            if ((*isAscendingFirst && arr[i] > arr[i + 1]) || (!*isAscendingFirst && arr[i] < arr[i + 1]))
                return i; // Return the index of the elbow
    }

    return -1; // If no elbow is found, the sequence is entirely increasing or decreasing or all elements are equal.
}

void elbowSort(int *arr, int size, bool isAscending)
{
    start_timer(&elbow_sort_timer);

    bool isAscendingFirst;
    int elbow = find_elbow(arr, size, &isAscendingFirst);

    // No elbow found
    if (elbow == -1)
    {
        // If the array is not in the correct order based on the initial trend and the desired isAscending order we reverse it
        if (isAscending != isAscendingFirst)
        {
            // Reverse the array to correct the order.
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
    int right = (elbow + 1) % size; // Start at the first element after the elbow (using modulo for circularity)

    int *sorted = malloc(size * sizeof(int)); // Allocate memory for the sorted array.
    int sorted_index = 0;                     // Index for the sorted array.

    // Merge the two parts of the bitonic sequence into the sorted array.
    while (sorted_index < size)
    {
        // Check the initial trend to determine which element to add.
        if (isAscendingFirst)
        {
            // Take the larger of the two elements.
            if (arr[left] >= arr[right])
            {
                sorted[sorted_index++] = arr[left];
                left = (left - 1 + size) % size; // Move left index backward (circularly).
            }

            else
            {
                sorted[sorted_index++] = arr[right];
                right = (right + 1) % size; // Move right index forward (circularly).
            }
        }

        else
        {
            // Take the smaller of the two elements.
            if (arr[left] <= arr[right])
            {
                sorted[sorted_index++] = arr[left];
                left = (left - 1 + size) % size; // Move left index backward (circularly).
            }

            else
            {
                sorted[sorted_index++] = arr[right];
                right = (right + 1) % size; // Move right index forward (circularly).
            }
        }
    }

    // Copy the sorted array back to the original array.
    if (isAscending == isAscendingFirst)
    {
        // Reverse the sorted array if it's not in the correct order
        for (int i = 0; i < size; i++)
        {
            arr[i] = sorted[size - i - 1];
        }
    }

    else
    {
        // Copy the sorted array directly if it's in the correct order.
        memcpy(arr, sorted, size * sizeof(int));
    }

    free(sorted); // Free the allocated memory.

    stop_timer(&elbow_sort_timer);
}
