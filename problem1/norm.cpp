#include "norm.h"
#include <cmath>
using namespace std;

double norm(double *array, int size)
{
    double max = 0.0;
    for (int i = 0; i < size; i++) {
//sum += pow(array[i],2.0);}

//sum = pow(sum,0.5);
	if (abs(array[i]) > max) {
	    max = abs(array[i]);
	}
    }
    return max;
}
