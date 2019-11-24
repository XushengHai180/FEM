#include"intergration.h"
#include"interpolation.h"
#include<stdio.h>
void GaussQuadratureRules_1D::calculate(int nig)
{
	switch(nig)
	{
	case 1:
		w[0]=2;
		xi[0]=0;
		break;
	case 2:
		w[0]=1;
		w[1]=1;
		xi[0]=-1/sqrt(3.0);
		xi[1]=1/sqrt(3.0);
		break;
	case 3:
		w[0]=8.0/9;
		w[1]=5.0/9;
		w[2]=5.0/9;
		xi[0]=0;
		xi[1]=-sqrt(3.0/5);
		xi[2]=sqrt(3.0/5);
		break;
	case 4:
		w[0]=(18.0+sqrt(30.0))/36;
		w[1]=(18.0+sqrt(30.0))/36;
		w[2]=(18.0-sqrt(30.0))/36;
		w[3]=(18.0-sqrt(30.0))/36;
		xi[0]=sqrt((3.0-2*sqrt(6.0/5))/7);
		xi[1]=-sqrt((3.0-2*sqrt(6.0/5))/7);
		xi[2]=sqrt((3.0+2*sqrt(6.0/5))/7);
		xi[3]=-sqrt((3.0+2*sqrt(6.0/5))/7);
		break;
	case 5:
		w[0]=128.0/225;
		w[1]=(322+13*sqrt(70.0))/900;
		w[2]=(322+13*sqrt(70.0))/900;
		w[3]=(322-13*sqrt(70.0))/900;
		w[4]=(322-13*sqrt(70.0))/900;
		xi[0]=0;
		xi[1]=sqrt(5-2*sqrt(10.0/7))/3;
		xi[2]=-sqrt(5-2*sqrt(10.0/7))/3;
		xi[3]=sqrt(5+2*sqrt(10.0/7))/3;
		xi[4]=-sqrt(5+2*sqrt(10.0/7))/3;
		break;
	}
}