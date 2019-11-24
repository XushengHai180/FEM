#ifndef ENTITY_H
#define ENTITY_H
#include<stdio.h>
#include"element.h"
class entity
{
public:
	int ne,nnpe,nn;//element number, node number per element,node number
	double l;//length of the bar
	void mesh(int ne_new,int nnpe_new,double l_new,double a,double e,double rho);
	element *ele;
	double **m,**k,*u,*v,*a;
	void assemble(void),assembleKOnly(void);
};
#endif