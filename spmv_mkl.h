#include <stddef.h>
#include <mkl.h>

int my_dense(const unsigned int n, const double mat[], const double vec[], double result[]);
int my_csr(const unsigned int n, const sparse_matrix_t *m, const double vec[], double result[]);
int my_csc(const unsigned int n, const MKL_INT *cols_start, const MKL_INT *rows_indx, const double *values, const double vec[], double result[]);
int my_coo(const unsigned int n, const unsigned int nnz, const MKL_INT *rows_indx, const MKL_INT *cols_indx, const double *values, const double vec[], double result[]);
