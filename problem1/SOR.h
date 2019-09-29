#ifndef WGT_JACOBI
#define WGT_JACOBI
#include "sparse_matrix.h"

void create_iterationsystem(sparse_matrix *A,double* x,double* b,double kfmode);
void SOR(sparse_matrix *A,double* x,double* b,double tolerance,double omega);


#endif
