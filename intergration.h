#ifndef INTERGRATION_H
#define INTERGRATION_H
#include"interpolation.h"
#include"math.h"
#include"stdlib.h"
class GaussQuadratureRules_1D
{
public:
	void calculate(int nig);//initial w and xi
	double w[5],xi[5];
};
#endif