#include"explicit.h"
void explicitSolver::initial(entity en_new,double t_new,double dt_new)//后续在这里加入边界条件参数
{
	en=en_new;
	t=t_new; dt=dt_new; 
	p=(double*)malloc(sizeof(double)*(1+en.nn));
	k11=(double*)malloc(sizeof(double)*(1+en.nn));
	k12=(double*)malloc(sizeof(double)*(1+en.nn));
	k21=(double*)malloc(sizeof(double)*(1+en.nn));
	k22=(double*)malloc(sizeof(double)*(1+en.nn));
	k31=(double*)malloc(sizeof(double)*(1+en.nn));
	k32=(double*)malloc(sizeof(double)*(1+en.nn));
	k41=(double*)malloc(sizeof(double)*(1+en.nn));
	k42=(double*)malloc(sizeof(double)*(1+en.nn));
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
	return;
}
void explicitSolver::refreshK()
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
void explicitSolver::run()
{
	FILE *fp;
	fp=fopen("result.txt","w");
	solveMatrix s;
	double tStepN=t/dt;
	double *F,*Fsegama,*middleU;
	F=(double *)malloc(sizeof(double)*(en.nn+1));
	Fsegama=(double *)malloc(sizeof(double)*(en.nn+1));
	middleU=(double *)malloc(sizeof(double)*(en.nn+1));
	en.m[0][0]=1000000000;	
	for(int tStep=1;tStep<=tStepN;tStep++)
	{
		//p[en.nn-1]=cos(tStep*dt*3.1415926/40);
		p[en.nn-1]=1;
		s.product(en.k,en.u,en.nn,Fsegama);
		for(int i=0;i<en.nn;i++)
		{
			k11[i]=dt*en.v[i];			F[i]=(p[i]-Fsegama[i])*dt;
		}
		s.solveMatrixEquation(en.m,F,en.nn,k12);
		for(int i=0;i<en.nn;i++) 
		{
			middleU[i]=en.u[i]+0.5*k11[i];
		}
		s.product(en.k,middleU,en.nn,Fsegama);
		for(int i=0;i<en.nn;i++)
		{
			k21[i]=dt*(en.v[i]+0.5*k12[i]);
			F[i]=(p[i]-Fsegama[i])*dt;
		}
		s.solveMatrixEquation(en.m,F,en.nn,k22);
		for(int i=0;i<en.nn;i++)
		{
			middleU[i]=en.u[i]+0.5*k21[i];
		}
	    s.product(en.k,middleU,en.nn,Fsegama);
		for(int i=0;i<en.nn;i++)
		{
			k31[i]=dt*(en.v[i]+0.5*k22[i]);
			F[i]=(p[i]-Fsegama[i])*dt;
		}
		s.solveMatrixEquation(en.m,F,en.nn,k32);
		for(int i=0;i<en.nn;i++)
		{ 
			middleU[i]=en.u[i]+k31[i];
		}
		s.product(en.k,middleU,en.nn,Fsegama);
		for(int i=0;i<en.nn;i++)
		{
			k41[i]=dt*(en.v[i]+k32[i]);
			F[i]=(p[i]-Fsegama[i])*dt;
		}
		s.solveMatrixEquation(en.m,F,en.nn,k42);
		for(int i=0;i<en.nn;i++)
		{
			en.u[i]+=(k11[i]+2*k21[i]+2*k31[i]+k41[i])/6;
			en.v[i]+=(k12[i]+2*k22[i]+2*k32[i]+k42[i])/6;
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
	free(middleU);
}