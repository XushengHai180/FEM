#ifndef IMPLICIT_H
#define IMPLICIT_H
#include"entity.h"
#include"boundaryCondition.h"//Î´Á¬½Ó
#include"solveMatrix.h"
class implicitSolver
{
public:
	entity en;
	double tAll,dt;//the time to calculate and the time step
	double *p,*essen;
	double beta,gama;
	void initial(entity en_new,double t_new,double dt_new);
	void run();
	void refreshK();
};
#endif