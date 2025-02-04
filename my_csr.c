#include <stddef.h>
#ifdef _GSL_
#include <gsl/gsl_spmatrix.h>
#include "spmv.h"
#endif
#ifdef _MKL_
#include <mkl.h>
#include "spmv_mkl.h"
#endif
//#include "csr.h"

#ifdef _GSL_
int my_csr(const unsigned int n, const gsl_spmatrix *m, const double vec[], double result[])
{
  // code your own solver
  unsigned int i, j=0, k=0;
  double tmp=0.0;
  //m->data = values
  //m->p = row_offsets
  //m->i = column_indices
  //Docs: https://www.gnu.org/software/gsl/doc/html/spmatrix.html#c.gsl_spmatrix_csr
  double *Md = m->data;
  int *Mp = m->p, *Mi = m->i;

  for(i=0; i < n; i++){
	j = Mp[i+1];
	for(k=Mp[i]; k < j; k++){
		tmp += Md[k] * vec[Mi[k]];
	}
	result[i] = tmp;
	tmp = 0;
  }

  return 0;
}
#endif
#ifdef _MKL_
int my_csr(const unsigned int n, const sparse_matrix_t *m, const double vec[], double result[])
{
  // code your own solver
  unsigned int i, j=0, k=0;
  double tmp=0.0;
  sparse_index_base_t indexing;
  MKL_INT nrows, ncols; 
  MKL_INT *rows_start, *rows_end, *cols_indx;
  double *values;

  //You CANNOT access to the elements of sparse_matrix_t by hand (are opaque)
  mkl_sparse_d_export_csr(*m, &indexing, &nrows, &ncols, &rows_start, &rows_end, &cols_indx, &values);
  for(i=0; i < n; i++){
	j = rows_start[i+1];
	for(k=rows_start[i]; k < j; k++){
		tmp += values[k] * vec[cols_indx[k]];
	}
	result[i] = tmp;
	tmp = 0;
  }
  
  return 0;
}
#endif
