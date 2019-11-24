#include"element.h"
double node::Constitution (double strain)//return etan
{
	double etan;
	//stress=E*((strain+1)*(strain+1)-1/(strain+1))/2;
	//etan=E+2*E/((1+strain)*(1+strain)*(1+strain));//rubber
	etan=1;
	//etan=3*E*((1+strain)*(1+strain)+1/((1+strain)*(1+strain)*(1+strain)*(1+strain)));//Ogden
	return etan;
}
void element::initial(double a,double e,double rho,double sX,double eX)
{
	node_obj=(node*)malloc(sizeof(node)*nnpe);
	for(int i=0;i<nnpe;i++)
	{
		node_obj[i].a=a;
		node_obj[i].E=e;
		node_obj[i].Etan=node_obj[i].Constitution(0);
		node_obj[i].rho=rho;
		node_obj[i].location=sX+(eX-sX)/(nnpe-1)*i;
	}
	if(nnpe<6) nig=nnpe;
	else nig=5;
	GaussQuadratureRules_1D c;
	c.calculate(nig);
	for(int i=0;i<nig;i++)
	{
		intergrater_w[i]=c.w[i];
		intergrater_xi[i]=c.xi[i];
	}
	m=(double**)malloc(sizeof(double*)*(nnpe+1));
	for(int i=0;i<nnpe+1;i++)
	{
		m[i]=(double*)malloc(sizeof(double)*(nnpe+1));
	}
	k=(double**)malloc(sizeof(double*)*(nnpe+1));
	for(int i=0;i<nnpe;i++)
	{
		k[i]=(double*)malloc(sizeof(double)*(nnpe+1));
	}
	calculateMassMatrix();
	calculateStiffnessMatrix();
	return;
}
shape* element::ShapeFunction(double xi)
{
	double *xiNew,*xNew;
	shape *sp;
	sp=(shape*)malloc(sizeof(shape)*nnpe);
	ParametricInterpolation_1D h;
	xiNew=(double*)malloc(sizeof(double)*(nnpe+1));
	xNew=(double*)malloc(sizeof(double)*(nnpe+1));
	for(int i=0;i<nnpe;i++)
	{
		xiNew[i]=node_obj[i].location;
		xNew[i]=node_obj[i].a;
	}
	h.initial(nnpe,xiNew,xNew);
	sp=h.ShapeFunction(xi);
	free(xiNew);
	free(xNew);
	return sp;
}
shape* element::ShapeFunction_firstDeri(double xi)
{
	double *xiNew,*xNew;
	shape *sp;
	sp=(shape*)malloc(sizeof(shape)*nnpe);
	ParametricInterpolation_1D h;
	xiNew=(double*)malloc(sizeof(double)*(nnpe+5));
	xNew=(double*)malloc(sizeof(double)*(nnpe+5));
	for(int i=0;i<nnpe;i++)
	{
		xiNew[i]=node_obj[i].location;
		xNew[i]=node_obj[i].a;
	}
	h.initial(nnpe,xiNew,xNew);
	sp=h.ShapeFunction_firstDeri(xi);
	free(xiNew);
	free(xNew);
	return sp;
}
double element::interpolater(double xi,int mode)//you can konow the walue at xi,  mode =1 is cross-sectional area; mode=2 is moudul Etan; =3 is density rho;
{
	double *xiNew,*xNew;
	ParametricInterpolation_1D in;
	xiNew=(double*)malloc(sizeof(double)*(nnpe+1));
	xNew=(double*)malloc(sizeof(double)*(nnpe+1));
	if(mode==1)
	{
		for(int i=0;i<nnpe;i++)
		{
			xiNew[i]=node_obj[i].location;
			xNew[i]=node_obj[i].a;
		}
		in.initial(nnpe,xiNew,xNew);
	}
	else if(mode==2)
	{
		for(int i=0;i<nnpe;i++)
		{
			xiNew[i]=node_obj[i].location;
			//printf("%d\n",i);
			xNew[i]=node_obj[i].Etan;
		}
		in.initial(nnpe,xiNew,xNew);
	}
	else if(mode==3)
	{
		for(int i=0;i<nnpe;i++)
		{
			xiNew[i]=node_obj[i].location;
			xNew[i]=node_obj[i].rho;
		}
		in.initial(nnpe,xiNew,xNew);
	}
	else if(mode==4)//strain is the gradient of u
	{
		for(int i=0;i<nnpe;i++)
		{
			xiNew[i]=node_obj[i].location;
			xNew[i]=node_obj[i].u;
		}
		in.initial(nnpe,xiNew,xNew);
		free(xiNew);
	    free(xNew);
		return(in.First_Derivatives(xi));
	}
	else
	{
		printf("请输入正确的mode");
		free(xiNew);
	    free(xNew);
		return -100;
	}
	free(xiNew);
    free(xNew);
	return(in.Lagrange_interp(xi));
}
double XiToX(double xMin,double xMax,double xi)//mapping from xi(-1,1) to x(location) 
{
	return(xMin+(xi+1)*(xMax-xMin)/2);
}
double XiToX_1d(double xMin,double xMax,double xi)
{
	return((xMax-xMin)/2);
}
double** element::calculateMassMatrix()
{
	shape *sp;
	for(int i=0;i<nnpe;i++)
	{
		for(int j=0;j<nnpe;j++)
		{
			m[i][j]=0;
			for(int k=0;k<nig;k++)
			{
				sp=ShapeFunction(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[k]));
				m[i][j]+=intergrater_w[k]*interpolater(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[k]),1)
					*interpolater(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[k]),3)
					*sp[i]*sp[j]
					*XiToX_1d(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[k]);
			}
		}
	}
	free(sp);
	return(m);
}
double** element::calculateStiffnessMatrix()
{
	shape *sp;
	for(int i=0;i<nnpe;i++)
	{
		for(int j=0;j<nnpe;j++)
		{
			k[i][j]=0;
			for(int l=0;l<nig;l++)
			{
				sp=ShapeFunction_firstDeri(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[l]));
				k[i][j]+=intergrater_w[l]*interpolater(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[l]),1)
					*interpolater(XiToX(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[l]),2)
						*sp[i]*sp[j]
						*XiToX_1d(node_obj[0].location,node_obj[nnpe-1].location,intergrater_xi[l]);
			}
		}
	} 
	free(sp);
	return(k);
}