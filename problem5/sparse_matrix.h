#ifndef SPARSE_MATRIX
#define SPARSE_MATRIX
#include <cstddef>

using namespace std;
typedef int row;
typedef int column;
struct sparse_matrix{
	int nnz ;//Number of non zero values
	int sizeofmatrix ;
	double *values ;
	int *col_index ;
	int *row_pointer ;
	};
#endif
