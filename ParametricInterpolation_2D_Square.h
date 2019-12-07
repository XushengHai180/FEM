#ifndef ParametricInterpolation_2D_Square_H
#define ParametricInterpolation_2D_Square_H

#include<stdlib.h>
#include<vector>
#include<iostream>

using namespace std;

class ParametricInterpolation_2D_Square
{
public:
	int NNPE_m, NNPE_n,NNPE;
	vector<double> xi,yita;//×ÔÈ»×ø±ê
	vector<vector<double>> value;
	ParametricInterpolation_2D_Square(int nnpe_m,int nnpe_n);
	double ShapeFunction(int m,int n,double xi0,double yita0);
	double ShapeFunction_firstDeri_xi(int m,int n,double xi0,double yita0);
	double ShapeFunction_firstDeri_yita(int m,int n,double xi0,double yita0);
	double Lagrange_interp(double xi0,double yita0);//lagrange interpolations, caculate the value of the approximately function at the point x
	double First_Derivatives_xi(double xi0,double yita0);
	double First_Derivatives_yita(double xi0,double yita0);
};

#endif