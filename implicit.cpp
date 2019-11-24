#include"implicit.h"
void implicitSolver::initial(entity en_new,double t_new,double dt_new)//后续在这里加入边界条件参数
{
	en=en_new;
	tAll=t_new; dt=dt_new; 
	p=(double*)malloc(sizeof(double)*(1+en.nn));
	for(int i=0;i<en.nn;i++)
	{
		p[i]=0;
		if(i=en.nn-1)
		{
			p[i]=1;
		}
	}
	essen=(double*)malloc(sizeof(double)*(1+en.nn));
	for(int i=0;i<en.nn;i++)
	{
		if(i==0)
		{
			essen[i]=0;
		}
		else
		{
			essen[i]=-1;//代表当地位移未知
		}
	}
	beta=0.5;
	gama=0.5;
	return;
}
void implicitSolver::refreshK()
{
	for(int i=0;i<en.ne;i++)
	{
		for(int j=0;j<en.nnpe;j++)
		{
			en.ele[i].node_obj[j].u=en.u[(en.nnpe-1)*i+j];
		}
	}
	for(int i=0;i<en.ne;i++)
	{
		for(int j=0;j<en.nnpe;j++)
		{
			en.ele[i].node_obj[j].strain=en.ele[i].interpolater(en.ele[i].node_obj[j].location,4);
		}
	}
	for(int i=0;i<en.ne;i++)
	{
		for(int j=0;j<en.nnpe;j++)
		{
			en.ele[i].node_obj[j].Etan=en.ele[i].node_obj[j].Constitution(en.ele[i].node_obj[j].strain);
		}
	}
	for(int i=0;i<en.ne;i++)
	{
		en.ele[i].calculateStiffnessMatrix();
	}
	en.assembleKOnly();
}
void implicitSolver::run()
{
	FILE *fp;
	fp=fopen("result.txt","w");
	solveMatrix s;
	double *F,*Fsegama,**resiGradient,*detA,*resiFu,*middleMwithA,*aOld;
	F=(double*)malloc(sizeof(double)*(en.nn+1));
	Fsegama=(double*)malloc(sizeof(double)*(en.nn+1));
	resiGradient=(double**)malloc(sizeof(double*)*(en.nn+1));
	for(int i=0;i<en.nn+1;i++)
	{
		resiGradient[i]=(double*)malloc(sizeof(double)*(en.nn+1));
	}
	resiFu=(double*)malloc(sizeof(double)*(en.nn+1));
	detA=(double*)malloc(sizeof(double)*(en.nn+1));
	middleMwithA=(double*)malloc(sizeof(double)*(en.nn+1));
	aOld=(double*)malloc(sizeof(double)*(en.nn+1));
	double tStepN=tAll/dt;
	en.m[0][0]=100000000;
	for(int tStep=1;tStep<=tStepN;tStep++)
	{
		double t=tStep*dt;
		p[en.nn-1]=1;
		//p[en.nn-1]=sin(tStep*dt*3.1415926/4);
		s.product(en.k,en.u,en.nn,Fsegama);
		for(int i=0;i<en.nn;i++) F[i]=p[i]-Fsegama[i];
		for(int i=0;i<en.nn;i++) aOld[i]=en.a[i];
		int n=0;//debug
		while(1)
		{
			n++;
			if(n>10000)
			{
				printf("不收敛");
				break;
			}
			for(int i=0;i<en.nn;i++)
			{
				for(int j=0;j<en.nn;j++)
			   {
					resiGradient[i][j]=en.m[i][j]+2*beta*en.k[i][j];
			   }
			}
			s.product(en.m,en.a,en.nn,middleMwithA);
			double error=0;
			for(int i=0;i<en.nn;i++)
			{
				resiFu[i]=F[i]-middleMwithA[i];
				error+=abs(resiFu[i]);
			}
			if(error/en.nn<0.01)
			{
				    printf("%lf ",t);
					break;
			}
			s.solveMatrixEquation(resiGradient,resiFu,en.nn,detA);
			for(int i=0;i<en.nn;i++)
			{
				en.a[i]+=detA[i];
			}
		}
		for(int i=0;i<en.nn;i++)
		{
			en.u[i]+=dt*en.v[i]+(dt*dt/2)*((1-2*beta)*aOld[i]+2*beta*en.a[i]);
			en.v[i]+=dt*((1-gama)*aOld[i]+gama*en.a[i]);
		}
		if(dt*tStep==1.25||dt*tStep==2.5||dt*tStep==5||dt*tStep==7.5||dt*tStep==10||dt*tStep==12.5||dt*tStep==15||dt*tStep==17.5||dt*tStep==20||dt*tStep==22.5||dt*tStep==25||dt*tStep==27.5||dt*tStep==30)
		{
		fprintf(fp,"\nt=%lf\n",dt*tStep);
		fprintf(fp,"\nu=\n");
		for(int i=0;i<en.nn;i++)
		{
			fprintf(fp,"%lf\n",en.u[i]);
		}
		fprintf(fp,"\nv=\n");
		for(int i=0;i<en.nn;i++)
		{
			fprintf(fp,"%lf\n ",en.v[i]);
		}
		fprintf(fp,"\n");
		}
		refreshK();
	}
	free(F);
	free(Fsegama);
	for(int i=0;i<en.nn+1;i++) free(resiGradient[i]);
	free(resiGradient);
}