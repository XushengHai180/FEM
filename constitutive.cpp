#include"1dBar.h"
#include<stdlib.h>
typedef double   strain;
typedef double   shearmodulus;
typedef double   stress;
stress ConstitutiveRubber(strain eps)
{
	stress segma;
	shearmodulus G;
	segma=G*((eps+1)*(eps+1)-1/(eps+1));
}