#ifndef SOLVEMATRIX_H
#define SOLVEMATRIX_H
#include<stdlib.h>
class solveMatrix
{
public:
   void product(double **matrix,double *vector,int n,double *result);
   void solveMatrixEquation(double **matrix,double *vector,int n,double* result);
};
#endif