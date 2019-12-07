#include"ParametricInterpolation_2D_Square.h"
double ShapeFunction_1D(int n,int nnpe,double xi0)
{
	vector<double> xi;
	for(int i=0;i<nnpe;i++) xi.push_back(-1+i*2.0/(nnpe-1));
	double denominator=1,numerator=1;//分母和分子
	for(int i=0;i<nnpe;i++)
	{
		if(i!=n)
		{
			denominator*=(xi[n]-xi[i]);
			numerator*=(xi0-xi[i]);
		}
	}
	return numerator/denominator;
}
double ShapeFunction_1D_firstDeri(int n,int nnpe,double xi0)
{
	vector<double> xi;
	for(int i=0;i<nnpe;i++) xi.push_back(-1+i*2.0/(nnpe-1));
	double denominator_xi=1,numerator_xi=0;//分母和分子
	for(int i=0;i<nnpe;i++)
	{
		if(i!=n)
		{
			denominator_xi*=(xi[n]-xi[i]);
		}
		if(i==n) continue;
		double numeratorSub=1;
		for(int j=0;j<nnpe;j++)
		{
			if(j!=i&&j!=n)
			{
				numeratorSub*=(xi0-xi[j]);
			}
		}
		numerator_xi+=numeratorSub;
	}
	return numerator_xi/denominator_xi;
}
ParametricInterpolation_2D_Square::ParametricInterpolation_2D_Square(int nnpe_m,int nnpe_n)
{
	NNPE_m=nnpe_m;
	NNPE_n=nnpe_n;
	NNPE=nnpe_m*nnpe_n;
	for(int i=0;i<nnpe_m;i++) yita.push_back(-1+i*2.0/(nnpe_m-1));
	for(int j=0;j<nnpe_n;j++) xi.push_back(-1+j*2.0/(nnpe_n-1));
	for(int i=0;i<nnpe_m;i++)
	{
		vector<double> path;
		for(int j=0;j<nnpe_n;j++)
			path.push_back(2);
		value.push_back(path);
	}
			

};
double ParametricInterpolation_2D_Square::ShapeFunction(int m,int n,double xi0,double yita0)
{
	return ShapeFunction_1D(m,NNPE_m,yita0)*ShapeFunction_1D(n,NNPE_n,yita0);
}
double ParametricInterpolation_2D_Square::ShapeFunction_firstDeri_xi(int m,int n,double xi0,double yita0)
{
	return ShapeFunction_1D_firstDeri(n,NNPE_n,xi0)*ShapeFunction_1D(m,NNPE_m,yita0);
}
double ParametricInterpolation_2D_Square::ShapeFunction_firstDeri_yita(int m,int n,double xi0,double yita0)
{
	return ShapeFunction_1D(n,NNPE_n,xi0)*ShapeFunction_1D_firstDeri(m,NNPE_m,yita0);
}
double ParametricInterpolation_2D_Square::Lagrange_interp(double xi0,double yita0)
{
	double ans=0;
	for(int i=0;i<NNPE_m;i++)
		for(int j=0;j<NNPE_n;j++)
			ans+=value[i][j]*ShapeFunction(i,j,xi0,yita0);
	return ans;
}
double ParametricInterpolation_2D_Square::First_Derivatives_xi(double xi0,double yita0)
{
	double ans=0;
	for(int i=0;i<NNPE_m;i++)
		for(int j=0;j<NNPE_n;j++)
			ans+=value[i][j]*ShapeFunction_firstDeri_xi(i,j,xi0,yita0);
	return ans;
}
double ParametricInterpolation_2D_Square::First_Derivatives_yita(double xi0,double yita0)
{
	double ans=0;
	for(int i=0;i<NNPE_m;i++)
		for(int j=0;j<NNPE_n;j++)
			ans+=value[i][j]*ShapeFunction_firstDeri_yita(i,j,xi0,yita0);
	return ans;
}