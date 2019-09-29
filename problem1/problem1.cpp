#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include "csr_mv.h"
#include "norm.h"
#include "SOR.h"
#define PI 3.1415926535897932
using namespace std;
int N;				//Number of intervals in (0,1)
double tolerance = pow(10.0, -8.0);
double omega = 1.9;
double *x;
double *b;
int main(int argc, char *argv[])
{
    double kfmode = 5.0;

    if (argc == 2) {
	N = atoi(argv[1]);
    } else {
	cout << "Error:Enter the number of intervals" << endl;
	exit(1);
    }

    sparse_matrix *A;
    A = new sparse_matrix;
    A->nnz = (3 * N) - 5;
    A->sizeofmatrix = N - 1;
    A->values = new double[(3 * N) - 5];
    A->col_index = new int[(3 * N) - 5];
    A->row_pointer = new int[N];

    x = new double[N - 1];
    b = new double[N - 1];
    create_iterationsystem(A, x, b, kfmode);
    SOR(A, x, b, tolerance, omega);

    delete[]A->values;
    delete[]A->col_index;
    delete[]A->row_pointer;
    getchar();
    return 0;

}
