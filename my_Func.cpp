
#include "cpp_headers.h"
double *Left_Preconditioned(int M, int N, int nz, int *AI_CSR, int *RJ, double *Rval, double *b, int M_kind)
{
int i,j,i1,i2;
int ith_column;
// calculates x=(M^-1)b, using triangular nature of M
// b: input, x: output
// or we can write as: Mx=b
// M_kind; 1: Jacobian, 2: Gaus-Seidel other: no Left_Preconditioned

double *x;
x = (double *) malloc( N * sizeof(double) );
double dummy=0.0;
if (M_kind == 1){
	i=0;
	for(i=0;i<N;i++){		
		i1=AI_CSR[i];
		i2=AI_CSR[i+1];
		for(j=i1;j<i2;j++){
			if(RJ[j]==i){
				ith_column=j;
				j=i2;
			}
		}
		x[i]=(b[i])/(Rval[ith_column]);
	}
	return x;		
}
else if(M_kind == 2){
	i=0;
	for(i=0;i<N;i++){
		dummy=0.0;		
		i1=AI_CSR[i];
		i2=AI_CSR[i+1];
		for(j=i1;j<i2;j++){
			if(RJ[j]<i){
				dummy+=(Rval[j])*x[RJ[j]];
			}
			if(RJ[j]==i){
				ith_column=j;
				j=i2;
			}
		}
		x[i]=(b[i]-dummy)/(Rval[ith_column]);
	}
	return x;
}
else
{
	return b;// do nothing
}

}
void CSR_Create(int M,int N, int nz, int *I, int *J, double *val,int *AI_CSR, int *RI, int *RJ, double *Rval)
{
// M: number of rows
// N: number of columns
// nz: number of nonzero elements
// I, J, val are inputs in COO
// AI_CSR, RI, RJ, Rval are outputs in CSR
int i,j,k,temp;
k=0;
for(i=0;i<M;i++){
        for(j=0;j<nz;j++){
            if( I[j] == i){
		RI[k]=I[j];
		RJ[k]=J[j];
		Rval[k]=val[j];
		k++;}
	}
    }
    k=0;temp=-1;
    for(i=0;i<nz;i++){
	if(RI[i]>temp){
		AI_CSR[k]=i;
		temp=RI[i];
		k++;}
    }
AI_CSR[N]=nz;
}

void CSC_Create(int M,int N, int nz, int *I, int *J, double *val,int *AJ_CSC, int *CI, int *CJ, double *Cval)
{
// M: number of rows
// N: number of columns
// nz: number of nonzero elements
// I, J, val are inputs in COO
// AI_CSR, RI, RJ, Rval are outputs in CSR
int i,j,k,temp;
k=0;
    for(i=0;i<N;i++){
        for(j=0;j<nz;j++){
            if( J[j] == i){
		CI[k]=I[j];
		CJ[k]=J[j];
		Cval[k]=val[j];
		k++;}
	}
    }
    k=0;temp=-1;
    for(i=0;i<nz;i++){
	if(CJ[i]>temp){
		AJ_CSC[k]=i;
		temp=CJ[i];
		k++;}
    }

     AJ_CSC[M]=nz;
}
void Residual_Ax_b(int M,int N, int nz, int *AI_CSR, int *RJ, double *Rval, int *AJ_CSC, int *CI, double *Cval, double *x, double *b,double *r, char SG)
{
// calculating r = b - Ax
	double *A_times_x;
	A_times_x = (double *) malloc( M * sizeof(double) );
	Ax(M, N, nz, AI_CSR, RJ, Rval, AJ_CSC, CI, Cval, x, A_times_x, SG);
	int i;
	for(i=0;i<M;i++)
		r[i]=b[i];
	for(i=0;i<M;i++)
		r[i]-=A_times_x[i];
	free(A_times_x);
}

void Ax(int M,int N, int nz, int *AI_CSR, int *RJ, double *Rval, int *AJ_CSC, int *CI, double *Cval, double *x, double *y, char SG)
{
// y = A*x
// M: number of rows
// N: number of columns
// nz: number of nonzero elements
// A is stored in CSR, CSC (in case of symmetric matrix) Fromat
	int i,i1,i2;
	for(i=0;i<M;i++){
		i1=AI_CSR[i];
		i2=AI_CSR[i+1]-1;
		y[i] = dot_prod(Rval,x,i1,i2,RJ);
		if (SG=='S'){
			i1=AJ_CSC[i];
			i2=AJ_CSC[i+1]-1;
			y[i] += dot_prod(Cval,x,i1+1,i2,CI);}
	}
}

void Ax_tridiagonal_A(int m, double *alpha, double *beta , double *x, double *y)
{
// y = A*x
// m: size of alpha
// akpha: array, diagonal of A, from 0 to m-1
// beta: array, off-diagonal of A, from 0 to m-2
	int i;
	for(i=0;i<m;i++)
		y[i] = x[i]*alpha[i];
	for(i=0;i<m-1;i++)
		y[i] += x[i+1]*beta[i];
	for(i=1;i<m;i++)
		y[i] += x[i-1]*beta[i-1];
}

double dot_prod(double *val,double *x,int i1,int i2, int *J)
{
	int i;double dummy=0.0;
	for(i=i1;i<=i2;i++){
		dummy += val[i]*x[J[i]];}
	return dummy;
}
double norm_2(int n,double *x)
{
	int i;
	double norm=0.0;
	for(i=0;i<n;i++)
		norm+=x[i]*x[i];
	norm = sqrt(norm);
	return norm;
}

double scalar_prod(int n,double *x,double *y)
{
	int i;
	double dummy=0.0;
	for(i=0;i<n;i++)
		dummy+=x[i]*y[i];
	return dummy;
}
double scalar_prod_partial(int i1,int i2,double *x,double *y)
{
	int i;
	double dummy=0.0;
	for(i=i1;i<i2;i++)
		dummy+=x[i]*y[i];
	return dummy;
}
double norm_A(int M,int N, int nz, int *AI_CSR, int *RJ, double *Rval, int *AJ_CSC, int *CI, double *Cval, double *x, char SG)
{
// returning sqrt[(Ax,x)]
	
	double norm=0.0;
	double *A_times_x;
	A_times_x = (double *) malloc( M * sizeof(double) );
	Ax(M, N, nz, AI_CSR, RJ, Rval, AJ_CSC, CI, Cval, x, A_times_x, SG);
	norm = scalar_prod(M, A_times_x, x);
	free(A_times_x);
	norm = sqrt(norm);
	return norm;
}

