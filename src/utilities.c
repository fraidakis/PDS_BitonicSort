// Implementation File: utils.c
#include "utilities.h"

// Define timer variables.
struct timer communication_timer;
struct timer elbow_sort_timer;
struct timer local_sort_timer;
struct timer total_timer;
struct timer sequential_timer;

bool parse_arguments(int argc, char **argv, int *q, int *p)
{
    // Check for correct number of arguments.
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <q> <p>\n", argv[0]);
        return false;
    }

    // Convert arguments to integers.
    *q = atoi(argv[1]);
    *p = atoi(argv[2]);

    return true;
}

bool validate_parameters(int numProcesses, int p, int rank)
{
    // Check if the number of processes running is same as processes requested.
    if (numProcesses != (1 << p))
    {
        if (rank == 0)
        {
            fprintf(stderr, "Error: Number of processes must be 2^p. Size = %d\n", numProcesses);
        }

        return false;
    }

    if (rank == 0)
    {
        printf("q: %d, p: %d\n\n", p, numProcesses);
    }

    return true;
}

void initialize_data(int *localArray, int size, int rank)
{
    // Seed the random number generator differently for each process.
    srand(time(NULL) + rank);

    // Fill the array with random numbers.
    for (int i = 0; i < size; i++)
    {
        localArray[i] = rand() % 1000;
    }
}

int *initialize_sequential_array(int rank, int numProcesses, int local_size, int *send_buffer)
{
    int *sequential_array = NULL;

    // Allocate memory on rank 0 only.
    if (rank == 0)
    {
        sequential_array = malloc(numProcesses * local_size * sizeof(int));
    }

    // Gather local data to rank 0.
    MPI_Gather(send_buffer, local_size, MPI_INT, sequential_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    return sequential_array;
}

void localSort(int *local_array, int local_size, int rank)
{
    start_timer(&local_sort_timer);

    // Sort in ascending order if rank is even, descending if odd
    if (rank & 1)
    {
        qsort(local_array, local_size, sizeof(int), compareDescending);
    }

    else
    {
        qsort(local_array, local_size, sizeof(int), compareAscending);
    }

    stop_timer(&local_sort_timer);
}

int getNumChunks(int local_size)
{
    const int threshold = (int)(1.5e6);

    // Calculate the required number of chunks directly
    int num_chunks = (local_size + threshold - 1) / threshold;

    // Ensure the number of chunks is a power of 2
    int power_of_2 = 1;

    while (power_of_2 < num_chunks)
    {
        power_of_2 *= 2;
    }

    return power_of_2;
}

bool getSortingDirection(int rank, int dimension)
{
    // Returns true (0) for ascending, false (1) for descending.
    return ((rank & (1 << dimension)) == 0);
}

void compare_and_swap(int *a, int *b, bool keep_min, int size)
{
    // Iterate through the arrays.
    for (int i = 0; i < size; i++)
    {
        if ((a[i] > b[i]) == keep_min) // If a[i] > b[i] and keep_min is true, swap
        {
            int temp = a[i];
            a[i] = b[i];
            b[i] = temp;
        }
    }
}

int compareAscending(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int compareDescending(const void *a, const void *b)
{
    return (*(int *)b - *(int *)a);
}

void start_timer(struct timer *t)
{
    t->start = MPI_Wtime();
}

void stop_timer(struct timer *t)
{
    t->end = MPI_Wtime();
    t->elapsed += t->end - t->start;
}

void gather_and_print_statistics(int rank, int numProcesses, long total_size)
{
    double global_communication, global_localSort, global_elbowMerge, global_all;
    
    // Reduce timer values across all processes.
    MPI_Reduce(&communication_timer.elapsed, &global_communication, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elbow_sort_timer.elapsed, &global_elbowMerge, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sort_timer.elapsed, &global_localSort, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_timer.elapsed, &global_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Average the timer values.
    global_localSort /= numProcesses;
    global_elbowMerge /= numProcesses;
    global_communication /= numProcesses;
    global_all /= numProcesses;

    // Print statistics on rank 0.
    if (rank == 0)
    {
        printf("Sorting %ld elements:\n", total_size);
        printf("LocalSort: %f seconds\n", global_localSort);
        printf("ElbowSort: %f seconds\n", global_elbowMerge);
        printf("Communication: %f seconds\n", global_communication);
        printf("Total: %f seconds\n\n\n", global_all);
    }
}

void verify_sequential_sort(int *send_buffer, int *sequential_array, long total_size, int local_size, int rank, int numProcesses)
{
    int *global_array = NULL;

    // Allocate memory for the global array on rank 0.
    if (rank == 0)
    {
        global_array = malloc(total_size * sizeof(int));
    }

    // Gather the local arrays to the global array on rank 0.
    MPI_Gather(send_buffer, local_size, MPI_INT, global_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform verification on rank 0.
    if (rank == 0)
    {
        start_timer(&sequential_timer); 
        // Sort the sequential array.
        qsort(sequential_array, total_size, sizeof(int), compareAscending);
        stop_timer(&sequential_timer); 


        printf("Sequential time : %f\n", sequential_timer.elapsed);

        // Compare the sorted global array with the sequentially sorted array.
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
        printf("%s sorting %ld elements across %d processes\n\n\n", is_sorted ? "SUCCESSFUL" : "FAILED", total_size, numProcesses);
        printf("Speedup: %f   ||   Efficiency: %f\n\n", sequential_timer.elapsed / total_timer.elapsed, (sequential_timer.elapsed / total_timer.elapsed) / numProcesses);
        free(global_array);
        free(sequential_array);
    }
}

void debug_print(int *local_array, int local_size, int rank, int numProcesses, int dimension, int distance)
{
    int *global_array = NULL;

    // Allocate memory for the global array on rank 0 only.
    if (rank == 0)
    {
        global_array = malloc(numProcesses * local_size * sizeof(int));
        if (global_array == NULL)
        {
            perror("malloc failed in debug_print"); // Handle potential malloc failure
            MPI_Abort(MPI_COMM_WORLD, 1);           // Abort MPI if allocation fails
        }
    }

    // Gather the local arrays to the global array on rank 0.
    MPI_Gather(local_array, local_size, MPI_INT, global_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the global array on rank 0.
    if (rank == 0)
    {
        int processor = 0;
        // Iterate through the global array and print elements.
        for (int i = 0; i < numProcesses * local_size; i++)
        {
            // Print a newline and process ID at the start of each process's data.
            if (i % local_size == 0)
            {
                printf("\nProcess %d: ", processor++);
            }
            
            printf("%d ", global_array[i]);
        }
        free(global_array); // Free the allocated memory.

        // Print dimension and distance information.
        printf("\n\nAfter: Dimension %d, Distance %d\n", dimension, distance);
    }
}
