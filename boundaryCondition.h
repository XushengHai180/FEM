#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H
#include<stdlib.h>
class boundaryCondition
{
public:
	double** pointLoad(double t);//return location and the force(not pressure) at time point t;
	double** esseBC(double t);//return location and the disaplacement at time point t; 
};
#endif