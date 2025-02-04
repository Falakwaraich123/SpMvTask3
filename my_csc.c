#ifdef _GSL_
#include "spmv.h"
#include <gsl/gsl_spmatrix.h>
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif

#ifdef _GSL_
int my_csc(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[])
{
  unsigned int i, j=0, k=0;
  double *Md = m->data;
  int *Mp = m->p, *Mi = m->i;
  for (i=0; i < n; i++)
  	result[i]=0;

 
  
  for(i=0; i < n; i++){
  	j = Mp[i+1];
	for(k = Mp[i]; k < j; k++){
			result[Mi[k]] += Md[k] * vec[i];
		}
	}

  return 0;
}
#endif
#ifdef _MKL_
int my_csc(const unsigned int n, const MKL_INT *cols_start, const MKL_INT *rows_indx, const double *values, const double vec[], double result[])
{
  unsigned int i, j=0, k=0;
  for (i=0; i < n; i++)
  	result[i]=0;

  for(i=0; i < n; i++){
  	j = cols_start[i+1];
		for(k=cols_start[i] ; k < j; k++){
			result[rows_indx[k]] += values[k] * vec[i];
		}
	}
  
  return 0;
}
#endif
