#include <stdlib.h>

/* Allocate memory for a rows*cols array of doubles.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 */
double* alloc_doublematrix(int cols, int rows)
{
    return (double*)calloc(rows*cols, sizeof(double));
} 

/* Allocate memory for a rows*cols array of chars. */
char* alloc_charmatrix(int cols, int rows)
{
    return malloc(rows*cols*sizeof(char));
} 

/* Free the memory of a matrix allocated with alloc_{double|char}matrix*/
void free_matrix(void *m)
{    
    free(m);
}
