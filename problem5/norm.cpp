#include "norm.h"
#include <cmath>
using namespace std;

double norm(double *array, int size){
double max=0.0;
for(int i=0;i<size;i++){
if(abs(array[i]) > max){
max =  abs(array[i]);
}
else continue;
}
return max;
}

