#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include<stdlib.h>
#include<math.h>
#define PI 3.1415926
typedef double shape;//the variable of shape function 
struct vector
{
	double x,y,z;
};
class ParametricInterpolation_1D
{
public:
	shape *ShapeFunction(double xi);
	shape *ShapeFunction_firstDeri(double xi);
	double Lagrange_interp(double xi);//lagrange interpolations, caculate the value of the approximately function at the point x
	double First_Derivatives(double xi);
	void initial(int nn_new,double *xi_new,double *x_new);
	ParametricInterpolation_1D(void);
private:
	int	 nn;//node number
	double xi[100];//node coordinate
	double x[100];//Image space
};
#endif