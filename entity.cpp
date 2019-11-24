#include"entity.h"
#include"stdlib.h"
void entity::mesh(int ne_new,int nnpe_new,double l_new,double a,double e,double rho)
{
	ne=ne_new;
	nnpe=nnpe_new;
	l=l_new;
	nn=ne*(nnpe-1)+1;
	ele=(element*)malloc(sizeof(element)*(ne+1));
	for(int i=0;i<ne;i++)
	{
		ele[i].nnpe=nnpe;
		ele[i].initial(a,e,rho,i*l/ne,(i+1.0)*l/ne);//assignment to each element
	}
	k=(double**)malloc(sizeof(double*)*(nn+50));
	for(int i=0;i<(nn+50);i++)
	{
		k[i]=(double*)malloc(sizeof(double*)*(nn+50));
	}
	m=(double**)malloc(sizeof(double*)*(nn+50));//这里开辟了较大的空间是为了抵消某种未知的bug（以后有时间再弄清楚）
	for(int i=0;i<(nn+50);i++)
	{
		m[i]=(double*)malloc(sizeof(double*)*(nn+50));
	}
}
void entity::assembleKOnly(void)
{
	for(int i=0;i<nn;i++)
	{
		for(int j=0;j<nn;j++)
		{
			k[i][j]=0;
		}
	}
	int xx,yy;
	for(int i=0;i<ne;i++)
	{
		for(int j=0;j<nnpe;j++)
    	{
			for(int l=0;l<nnpe;l++)
			{
			    xx=i*(nnpe-1)+j;
				yy=i*(nnpe-1)+l;
				k[xx][yy]=k[xx][yy]+ele[i].k[j][l];
			}
	     }
	}
	return;
}
void entity::assemble(void)
{
	for(int i=0;i<nn;i++)
	{
		for(int j=0;j<nn;j++)
		{
			m[i][j]=0;
			k[i][j]=0;
		}
	}
	int xx,yy;
	for(int i=0;i<ne;i++)
	{
		for(int j=0;j<nnpe;j++)
    	{
			for(int l=0;l<nnpe;l++)
			{
			    xx=i*(nnpe-1)+j;
				yy=i*(nnpe-1)+l;
				m[xx][yy]=m[xx][yy]+ele[i].m[j][l];
				k[xx][yy]=k[xx][yy]+ele[i].k[j][l];
			}
	     }
	}
	u=(double*)malloc(sizeof(double)*(nn+1));
	v=(double*)malloc(sizeof(double)*(nn+1));
	a=(double*)malloc(sizeof(double)*(nn+1));
	for(int i=0;i<nn;i++)
	{
		u[i]=0;
		v[i]=0;
		a[i]=0;
	}
	return;
}