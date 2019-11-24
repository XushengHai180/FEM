#ifndef ELEMENT_H
#define ELEMENT_H
#include<stdio.h>
#include<stdlib.h>
#include"intergration.h"
#include"interpolation.h"
class node
{
public:
	double location;//location 
	double a;//cross-sectional area
	double E,rho,u,strain;//moudul,density,displacement
	double Constitution(double strain);//return stress
	double Etan;
};
class element
{
public:
	node *node_obj;
	void initial(double a,double e,double rho,double sX,double eX);//initial node in this element
	int	nnpe,nig;//node number per element
	double intergrater_w[5],intergrater_xi[5],intergrater_w_1d[5],intergrater_xi_1d[5];
	double **m,**k;
	double** calculateStiffnessMatrix();
	double** calculateMassMatrix();
	double* ShapeFunction(double xi);
	double* ShapeFunction_firstDeri(double xi);
	double interpolater(double xi,int mode);//you can konow the walue at xi,  mode =1 is cross-sectional area; mode=2 is moudul Etan; =3 is density rho; mode=4 is strain
};
#endif