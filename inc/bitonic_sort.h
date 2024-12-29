// Header File: bitonic_sort.h
#ifndef BITONIC_SORT_H
#define BITONIC_SORT_H

#include "utilities.h"

/**
 * Sorts a local array using the Bitonic Sort algorithm.
 * @param send_buffer - The local array to sort.
 * @param recv_buffer - A buffer for receiving data during communication.
 * @param local_size - Number of elements in the local array.
 * @param rank - Rank of the current MPI process.
 * @param p - Log2 of the number of total MPI processes.
 */
void bitonicSort(int *send_buffer, int *recv_buffer, int local_size, int rank, int p);

/**
 * Exchanges data and performs compare-and-swap operations in the hypercube network.
 * @param send_buffer - Local array to send.
 * @param recv_buffer - Buffer for receiving data.
 * @param local_size - Number of elements in each buffer.
 * @param partner - Rank of the partner process for exchange.
 * @param keep_min - If true, keeps the minimum values in 'send_buffer'; otherwise, keeps the maximum.
 */
void hypercube_exchange_cmpswp(int *send_buffer, int *recv_buffer, int local_size, int partner, bool keep_min);

/**
 * Merges a cyclic bitonic sequence using the "elbow" strategy.
 * @param arr - The array to merge.
 * @param size - Number of elements in the array.
 * @param isAscending - Sorting direction: true for ascending, false for descending.
 */
void elbowSort(int *arr, int size, bool isAscending);

#endif // BITONIC_SORT_H
