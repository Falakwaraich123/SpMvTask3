#ifdef _GSL_
#include "spmv.h"
#include <stddef.h>
#include <gsl/gsl_spmatrix.h>
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif

#ifdef _GSL_
int my_coo(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[])

{
  unsigned int i, size = m->nz;
 
  for (i=0; i < n; i++)
  	result[i]=0;
  for(i=0; i < size; i++)
	result[m->i[i]]	+= m->data[i] * vec[m->p[i]];
  
  return 0;
}
#endif

#ifdef _MKL_
int my_coo(const unsigned int n, const unsigned int nnz, const MKL_INT *rows_indx, const MKL_INT *cols_indx, const double *values, const double vec[], double result[])
{
  unsigned int i;
  for (i=0; i < n; i++)
  	result[i]=0;
  for(i=0; i < nnz; i++)
	result[rows_indx[i]] += values[i] * vec[cols_indx[i]];
  
  return 0;
}
#endif
