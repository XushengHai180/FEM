#include"solveMatrix.h"
#include<stdio.h>
#include<stdlib.h>
void solveMatrix::product(double **matrix,double *vector,int n,double *result)
{
	//result=(double*)malloc(sizeof(double)*n);
	for(int i=0;i<n;i++)
	{
		result[i]=0;
		for(int j=0;j<n;j++)
		{
			result[i]+=matrix[i][j]*vector[j];
		}
	}
	return;
}
void solveMatrix::solveMatrixEquation(double **matrix,double *vector,int n,double* result)
{
	double error;
	for(int i=0;i<n;i++)
	{
			result[i]=0.1;
	}
	for(int k=0;k<1000;k++)
	{
		for(int i=0;i<n;i++)
		{
			double middle=0;
			for(int j=0;j<n;j++)
			{
				if(i!=j)
				{
					middle+=matrix[i][j]*result[j];
				}
			}
			result[i]=(vector[i]-middle)/matrix[i][i];
		}
	}
	return;
}
