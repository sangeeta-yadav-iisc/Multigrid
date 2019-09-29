#include "csr_mv.h"
#include "sparse_matrix.h"
using namespace std;

double csr_mv(sparse_matrix *A,double *x, row i){

double y = 0.0;
for(int j= A->row_pointer[i];j<=(A->row_pointer[i+1]-1);j++){
if(i == A->col_index[j]){
continue;}
else{
y = y + A->values[j]*x[A->col_index[j]]; }
}
return y;

}



