#ifndef EXPLICIT_H
#define EXPLICIT_H
#include"entity.h"
#include"boundaryCondition.h"//Î´Á¬½Ó
#include"solveMatrix.h"
class explicitSolver
{
public:
	entity en;
	double t,dt;//the time to calculate and the time step
	double *k11,*k12,*k21,*k22,*k31,*k32,*k41,*k42;
	double *p,*essen;
	void initial(entity en_new,double t_new,double dt_new);
	void run();
	void refreshK();
};
#endif