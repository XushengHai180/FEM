#include"boundaryCondition.h"
double** boundaryCondition::pointLoad(double t)//return location and the force(not pressure) at time point t;//the location ranges from 0 to 1
{
	double** p;
	p=(double**)malloc(sizeof(double*)*2);
	p[0]=(double*)malloc(sizeof(double)*3);
	p[1]=(double*)malloc(sizeof(double)*3);
	p[0][0]=1;
	p[0][1]=1;
	return p;
}
double** boundaryCondition::esseBC(double t)//写的不好，要修改！
{
	double** e;
	e=(double**)malloc(sizeof(double*)*2);
	e[0]=(double*)malloc(sizeof(double)*3);
	e[1]=(double*)malloc(sizeof(double)*3);
	e[0][0]=0;
	e[0][1]=0;
	return e;
}
