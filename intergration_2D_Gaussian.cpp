#include<intergration_2D_Gaussian.h>
GaussQuadratureRules_2D::GaussQuadratureRules_2D(int nig)
{
	if(nig!=1||nig!=4||nig!=9||nig!=16||nig!=25) printf("GaussQuadratureRules_2D--积分点数目非法");
	if(nig==1)
	{
		vector<double> path;
		path.push_back(0);
		path.push_back(0);
		path.push_back(2);
		xiYitaWeight.push_back(path);
	}
	else if(nig==4)
	{
		double xiValue=0.577350269189626;
		double xiS[4]={-1,-1,1,1},yitaS[4]={-1,1,-1,1};
		for(int i=0;i<4;i++)
		{
			vector<double> path;
			path.push_back(xiS[i]*xiValue);
			path.push_back(yitaS[i]*xiValue);
			path.push_back(1);
			xiYitaWeight.push_back(path);
		}
	}
	else if(nig==9)
	{
		double xiValue=0.774596669241483;
		double weightValue0=0.555555555555556;
		double weightValue1=0.888888888888889;
		double xi[3]={-1*xiValue,0,xiValue},yita[3]={-1*xiValue,0,xiValue},weight[3]={weightValue0,weightValue1,weightValue0};
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			{
				vector<double> path;
				path.push_back(xi[i]);
				path.push_back(yita[j]);
				path.push_back(weight[i]*weight[j]);
				xiYitaWeight.push_back(path);
			}
	}
	else if(nig==16)
	{
		double xiValue0=0.861136311594053;
		double xiValue1=0.339981043584856;
		double weightValue0=0.347854845137454;
		double weightValue1=0.652145154862546;
		double xi[4]={-1*xiValue0,-1*xiValue1,xiValue1,xiValue0},yita[4]={-1*xiValue0,-1*xiValue1,xiValue1,xiValue0},weight[4]={weightValue0,weightValue1,weightValue1,weightValue0};
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				vector<double> path;
				path.push_back(xi[i]);
				path.push_back(yita[j]);
				path.push_back(weight[i]*weight[j]);
				xiYitaWeight.push_back(path);
			}
	}
	else if(nig==25)
	{
	}
}