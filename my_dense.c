#ifdef _GSL_
  #include "spmv.h"
#endif
#ifdef _MKL_
  #include "spmv_mkl.h"
#endif
#include <math.h>
#include <stdio.h>

int my_dense(const unsigned int n, const double *restrict mat, const double *restrict vec, double *restrict result)
{
  unsigned int i, j;
  double tmp = 0.0;
  double *m = (double *)__builtin_assume_aligned(mat, 32);
  double *v = (double *)__builtin_assume_aligned(vec, 32);
  double *r = (double *)__builtin_assume_aligned(result, 32);
  double *M = m;
  for (i=0; i < n; i++){
    for (j=0; j < n ; j++){
	tmp += M[j] * v[j];
    }
    M += n;
	r[i] = tmp;
	tmp=0;
  }

  return 0;
}
