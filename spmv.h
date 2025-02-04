
#include <stddef.h>
#include <gsl/gsl_spmatrix.h>

int my_dense(const unsigned int n, const double mat[], const double vec[], double result[]);
int my_csr(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[]);
int my_csc(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[]);
int my_coo(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[]);
