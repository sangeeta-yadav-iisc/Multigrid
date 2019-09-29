#include "sparse_matrix.h"
#include "csr_mv.h"
#include "norm.h"
#include  <cmath>
#include  <fstream>
#include <iostream>
#include <cblas.h>
#define PI 3.1415926535897932
void printarray(double *x, int n)
{
    cout << endl;
    for (int i = 0; i < n; i++) {
	printf("%lf ", x[i]);
    }
    cout << endl;
}

using namespace std;
void
create_iterationsystem(sparse_matrix * A, double *x, double *b,
		       double kfmode)
{
    double h = 1 / (double) ((A->sizeofmatrix) + 1);
/****************Sparse matrix allocation***************************/
    A->values[0] = (2 + pow(h, 2.0)) / pow(h, 2.0);
    A->values[1] = -1 / pow(h, 2.0);

    for (int i = 2; i <= (A->nnz) - 3; i += 3) {
	A->values[i] = -1 / pow(h, 2);
	A->values[i + 1] = (2 + pow(h, 2.0)) / pow(h, 2.0);
	A->values[i + 2] = -1 / pow(h, 2);

    }

    A->values[A->nnz - 2] = -1 / pow(h, 2.0);
    A->values[A->nnz - 1] = (2 + pow(h, 2.0)) / pow(h, 2.0);
/****************Sparse matrix allcoation****************************/

/****************Column Index Allocation**************************/
    A->col_index[0] = 0;
    A->col_index[1] = 1;

    for (int i = 2; i <= (A->nnz) - 3; i += 3) {
	A->col_index[i] = (i - 2) / 3;
	A->col_index[i + 1] = A->col_index[i] + 1;
	A->col_index[i + 2] = A->col_index[i] + 2;
    }

    A->col_index[A->nnz - 2] = A->sizeofmatrix - 2;	//N-3
    A->col_index[A->nnz - 1] = A->sizeofmatrix - 1;	//N-2
/**********************Column Index Allocation***********************/

/********************Row Pointer Allocation**********************/
    A->row_pointer[0] = 0;
    A->row_pointer[A->sizeofmatrix] = (A->nnz) + 1;

    for (int i = 1; i < (A->sizeofmatrix); i++) {
	A->row_pointer[i] = (3 * i) - 1;
    }

/************************Row Pointer Allocation******************/
    for (int i = 0; i < (A->sizeofmatrix); i++) {
	x[i] = sin((i + 1) * kfmode * PI / ((A->sizeofmatrix) + 1));
	b[i] = 1.0;
    }
}

void
SOR(sparse_matrix * A, double *x, double *b, double tolerance,
    double omega)
{
    double h = 1 / (double) ((A->sizeofmatrix) + 1);
    double *y = new double[A->sizeofmatrix];
    int iterations = 0;
    double rel_error;
    ofstream datafile;

    datafile.open("datafile.dat");
    do {
	for (int row_index = 0; row_index < (A->sizeofmatrix); row_index++) {
	    y[row_index] = x[row_index];
	    x[row_index] =omega * ((-csr_mv(A, x, row_index) + b[row_index]) *
			 pow(h,2.0) / (2 + pow(h,2.0))) + (1 - omega)*x[row_index];
	}
	cblas_daxpy(A->sizeofmatrix, -1.0, x, 1, y, 1);
	rel_error = norm(y, A->sizeofmatrix) / norm(x, A->sizeofmatrix);
	iterations++;
    }
    while ((rel_error > tolerance) && (iterations < 100000));	//(rel_error > tolerance);
    cout << "Converged in " << iterations << " iterations" << endl;
    for (int i = 0; i < (A->sizeofmatrix); i++) {
	datafile << (i + 1) * h << "\t" << x[i] << endl;
    }

    FILE *gnuplot = popen("gnuplot", "w");
    fprintf(gnuplot,
	    "set key top right\n set xlabel \'x\' \n set ylabel \"u(x)\"\n");
    fprintf(gnuplot,
	    "set terminal wxt 0 \n plot \"datafile.dat\" using 1:2  with linespoints title \"Solution with %d intervals\"\n ",
	    A->sizeofmatrix + 1);
    fflush(gnuplot);
    datafile.close();

}
