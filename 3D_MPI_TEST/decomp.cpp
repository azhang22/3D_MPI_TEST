#include <math.h>
#include <stdio.h>
#include "mygeneral.h"
void mpe_decomp1d(int n,int numprocs,int myid,int &s,int &e)
{  	
	int nlocal,deficit;
	nlocal=n/numprocs;
	s=myid*nlocal+1;
	deficit=n%numprocs;
	s=s+qu_min(myid,deficit);
	if(myid<deficit) nlocal=nlocal+1;
	e=s+nlocal-1;
	if(e>n || myid==numprocs-1)e=n;
}
int qu_min(int a,int b)
{
	if (a>=b) return b;
	else return a;
}