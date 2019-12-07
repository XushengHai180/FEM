#ifndef intergration_2D_Gaussian_H
#define intergration_2D_Gaussian_H

#include<vector>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

class GaussQuadratureRules_2D
{
public:
	vector<vector<double>> xiYitaWeight;
	GaussQuadratureRules_2D(int nig);
};
#endif