#include<stdio.h>
#include"interpolation.h"
#include"intergration.h"
#include"element.h"
#include"entity.h"
#include"explicit.h"
#include"solveMatrix.h"
#include"implicit.h"
int main()//内存管理很垃圾
{
	entity me;
	explicitSolver ex;
	implicitSolver im;
	me.mesh(4,6,10,1,1,1);//ne,nnpe,l,a,e,rho
	me.assemble();
	ex.initial(me,40,0.01);
	ex.run();
	//im.initial(me,1.2,0.1);
	//im.run();
	system("pause");
	return 0;
}