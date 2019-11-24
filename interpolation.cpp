#include"interpolation.h"
#include<stdio.h>
shape* ParametricInterpolation_1D::ShapeFunction(double xi)
{
   shape *sp;
   sp=(shape*)malloc(sizeof(shape)*nn);
   	for(int i=0;i<ParametricInterpolation_1D::nn;i++)//the shape function i
	{
		double denominator=1,numerator=1;//分母和分子
		for(int j=0;j<ParametricInterpolation_1D::nn;j++)
		{
			if(i!=j)
			{
				numerator*=(xi-ParametricInterpolation_1D::xi[j]);
				denominator*=(ParametricInterpolation_1D::xi[i]-ParametricInterpolation_1D::xi[j]);
			}
		}
		sp[i]=numerator/denominator;
	}
	return(sp);
}
shape* ParametricInterpolation_1D::ShapeFunction_firstDeri(double xi)
{
	shape *sp;
    sp=(shape*)malloc(sizeof(shape)*nn);
	for(int i=0;i<ParametricInterpolation_1D::nn;i++)//the shape function i
	{
		double denominator=1,numerator=0;
		for(int j=0;j<ParametricInterpolation_1D::nn;j++)//calculate the deno first（different from 0 deri）
		{
			if(i!=j)
			{
				denominator*=(ParametricInterpolation_1D::xi[i]-ParametricInterpolation_1D::xi[j]);
			}
		}
		for(int j=0;j<ParametricInterpolation_1D::nn;j++)
		{
			if(i==j) continue;
			double sum=1;//intermediate variable
			for(int k=0;k<ParametricInterpolation_1D::nn;k++)
			{
				if(i!=k&&j!=k)
				{
					sum*=(xi-ParametricInterpolation_1D::xi[k]);
				}
			}
			numerator+=sum;
		}
		sp[i]=numerator/denominator;
	}
	return sp;
}
double ParametricInterpolation_1D::Lagrange_interp(double xi)//作用是根据xi[i]来求拟合曲线
{
	double x=0;
	shape *sp;
	sp=ShapeFunction(xi);
	for(int i=0;i<ParametricInterpolation_1D::nn;i++)//shape function linear superposition
	{
		x+=ParametricInterpolation_1D::x[i]*sp[i];
	}
	return x;
}
double ParametricInterpolation_1D::First_Derivatives(double xi)
 {
	double x=0;
	shape *sp;
	sp=ShapeFunction_firstDeri(xi);
	for(int i=0;i<ParametricInterpolation_1D::nn;i++)//shape function linear superposition
	{
		x+=ParametricInterpolation_1D::x[i]*sp[i];
	}
	return x;
 }
ParametricInterpolation_1D::ParametricInterpolation_1D(void)//构造函数,initialization the parameter//没有乘系数
{
	ParametricInterpolation_1D::nn=11;
	for(int i=0;i<ParametricInterpolation_1D::nn;i++)
	{
			ParametricInterpolation_1D::xi[i]=i/5.0-1;
			ParametricInterpolation_1D::x[i]=xi[i];//没有乘系数
	}
}
void ParametricInterpolation_1D::initial(int nn_new,double *xi_new,double *x_new)
{
	nn=nn_new;
	for(int i=0;i<nn_new;i++)
	{
		xi[i]=xi_new[i];
		x[i]=x_new[i];
	}
}