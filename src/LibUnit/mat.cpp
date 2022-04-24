#include "rtklib.h"
#include <iostream>
#include <Dense>
//#include"Dense"
//#include <thirdparty/Eigen/Dense>
using namespace Eigen;
using namespace std;

extern double matdet(double *A, int n)
{
    MatrixXd B;
    B = Map<Matrix<double, Dynamic, Dynamic, ColMajor>>(A, n, n);

    double b = 0.0;
    b = B.determinant();

    return b;
}