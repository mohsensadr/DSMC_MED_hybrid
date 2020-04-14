#include "cpp_headers.h"

#include <stdio.h>
#include <iostream>
#include<math.h>
using namespace std;

void rotate(double *x0,double *y0, double theta, double *x, double *y){
	*x = cos(theta)*(*x0)-sin(theta)*(*y0);
	*y = sin(theta)*(*x0)+cos(theta)*(*y0);
}

double geometric_sum(int n, double s){
// sum_{i=0}^{i=n} s^i
	if(n<0)
		return 0.0;
	else{
		if(fabs(s-1.0)<1e-15)
			return (n+1)*s;
		else
			return (1.0-pow(s,n+1.0))/(1.0-s);
	}
}
double volume_conical_frustum(double R1, double R2, double h){
// 	http://mathworld.wolfram.com/ConicalFrustum.html
return acos(-1.0)*h*(R1*R1 + R1*R2 + R2*R2)/3.0;
}
void Thompson_3d(double *a,double *b,double *c,double *x, double *RHS,int n)
{
// a_n is the main diagonal
// b_n-1 is lower-off diagonal
// c_n-1 is upper-off diagonal
    int i;
    double r;
    double Factor;
    //forward elimination
    for(i=1;i<n;i++)
    {
     	Factor=(1.0/a[i-1])*b[i-1];
        a[i]=a[i]-c[i-1]*Factor;
        RHS[i]=RHS[i]-Factor*RHS[i-1];
    }
    for(i=n-1;i>=0;i--)
    {
     	r=RHS[i];
        if(i!=n-1)
        {
            r=r-c[i]*x[i+1];
        }
	x[i]=r/a[i];
    }
}

void linear_fit(double n1, double r1, double n2, double r2, double n3, double r3, double *a, double *b){
      double sum_xi = r1+r2+r3;
      double sum_yi = n1+n2+n3;
      double sum_xixi = r1*r1+r2*r2+r3*r3;
      double sum_xiyi = n1*r1+n2*r2+n3*r3;
// n = a*r+b
      *b = ( sum_yi*sum_xixi-sum_xi*sum_xiyi )/( 3.0*sum_xixi-sum_xi*sum_xi );
      *a = ( 3.0*sum_xiyi-sum_xi*sum_yi )/( 3.0*sum_xixi-sum_xi*sum_xi );
}


double maxx(double *a,int n)
{
    double b;
    b=a[0];
    for(int i=0;i<n;i++)
    {
        if(b<a[i])
        {
            b=a[i];
        }
    }
    return b;
}

int max_int(int *a,int n)
{
    int b;
    b=a[0];
    for(int i=0;i<n;i++)
    {
        if(b<a[i])
        {
            b=a[i];
        }
    }
    return b;
}

double min(double a, double b)
{
  if(a>b)
    return b;
  else
    return a;
}
double max(double a, double b)
{
  if(a>b)
    return a;
  else
    return b;
}
int maxx_magnittude(double *a,int n)
{
    double b;int l=0;
    b=fabs(a[0]);
    for(int i=0;i<n;i++)
    {
        if(fabs(b)<fabs(a[i]))
        {
            b=fabs(a[i]);
            l=i;
        }
    }
    return l;
}


double standard_deviation(double data[], int n, double *mean)
{
    double  sum_deviation=0.0;
    mean[0] = 0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean[0]+=data[i];
    }
    mean[0]=mean[0]/(1.0*n);
    for(i=0; i<n;++i)
    sum_deviation+=(data[i]-mean[0])*(data[i]-mean[0]);
    return sqrt(sum_deviation/(1.0*n));
}
double mean(double *data, int n)
{
    double mean=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=data[i];
    }
    mean=mean/n;
    return mean;
}
double sum(double *data, int n)
{
    double dummy=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        dummy += data[i];
    }
    return dummy;
}
double Scalling(double *x,int i, int n)
{
    int b;
    b=maxx_magnittude(x,n);
    return fabs(x[i])/fabs(x[b]);
}

void Gauss_Jordan(double **A,double *RHS, int n, double *x)
{
//    double *x;
//    x=new double [n];
    double Factor,*z,b,c;
    int i,j,k=0, z_size;
    // Forward Elimination
    //cout<<"\nForwards Elimination\n";
    for(i=0;i<n-1;i++)
    {

        z_size=n;
	z = (double *) malloc( z_size * sizeof(double) );
	for(int ii=0; ii< z_size; ii++)
	  z[ii] = 0.0;
	k=0;
        for(j=i;j<n;j++)
        {
            z[k]=Scalling(A[j],i,n);
            k++;
        }
        k=maxx_magnittude(z,n);
	free(z);
        k=i+k;
        if(i!=k && k<n)
        {
            for(j=0;j<n;j++)
            {
                b=A[i][j];
                A[i][j]=A[k][j];
                A[k][j]=b;
            }
            c=RHS[i];
            RHS[i]=RHS[k];
            RHS[k]=c;
            //cout<< "\n** After Pivoting **"<<endl;
            //Printing_Matrix(A,RHS,n);
            //cout<<"\n\n";
        }

        for(j=i+1;j<n;j++)
        {
            Factor=(1.0/A[i][i])*A[j][i];
            for(k=0;k<n;k++)
            {
                A[j][k]=A[j][k]-Factor*A[i][k];
            }
            RHS[j]=RHS[j]-Factor*RHS[i];
        }

        //cout<< "*****  "<<i<<endl;
        //Printing_Matrix(A,RHS,n);
    }
//    cout<<"\n\n ** After Forward Elimination with pivoting and scalling **\n\n";
//    Printing_Matrix(A,RHS,n);
    //Backward Substitution
    double a;
    for(i=n-1;i>=0;i--)
    {
        a=0;
        for(j=i+1;j<n;j++)
        {
            a=a+A[i][j]*x[j];
        }
        x[i]=(RHS[i]-a)/(A[i][i]);
    }
   // return x;
}


int Factoriel(int n)
{
    double x=1;
    for(int i=1;i<n+1;i++)
    {
        x=x*i;
    }
    return x;
}

void curve_fit(double *x, double *y, int n, int *indices, double *a, double *b)
{
  double  sumx=0.0, sumx2=0.0, sumy=0.0, sumxy=0.0;
  int i, ii;
    for(ii=0;ii<=n-1;ii++)
    {
        i = indices[ii];
        sumx=sumx +x[i];
        sumx2=sumx2 +x[i]*x[i];
        sumy=sumy +y[i];
        sumxy=sumxy +x[i]*y[i];

    }
    *a=((sumx2*sumy -sumx*sumxy)*1.0/(n*sumx2-sumx*sumx)*1.0);
    *b=((n*sumxy-sumx*sumy)*1.0/(n*sumx2-sumx*sumx)*1.0);
}

// void curve_fit_U(struct CELLS *cells,  struct BOX *box, struct GAS *gas, double *slope, double *a)
// {
//   double *x,*y;
//   int n;
//   n = (*box).N[1];
//   x = (double *) malloc( n * sizeof(double) );
//   y = (double *) malloc( n * sizeof(double) );
//   int i;
//   int beg=0, end=n;
//   int many = end - beg;
//   for(i=beg; i<end-1; i++){
//     x[i] = 0.5* ( cells[i].dim[2]+cells[i].dim[3] );
//     y[i] =   cells[i].U[0];
//   }
//   double  sumx=0.0, sumx2=0.0, sumy=0.0, sumxy=0.0;
//     for(i=beg;i< end-1;i++)
//     {
//         sumx=sumx +x[i];
//         sumx2=sumx2 +x[i]*x[i];
//         sumy=sumy +y[i];
//         sumxy=sumxy +x[i]*y[i];
//
//     }
//     *a=((sumx2*sumy -sumx*sumxy)*1.0/(many*sumx2-sumx*sumx)*1.0);
//     *slope=((many*sumxy-sumx*sumy)*1.0/(many*sumx2-sumx*sumx)*1.0);
// }

// void mean_mu_cell(struct CELLS *cells,  struct BOX *box, struct GAS *gas, double *mean_mu_s, int new_step)
// {
//   int i, k=0;
//   double dummy = 0.0;
//     for(i=2; i<(*box).N[1]-3; i++)
//     {
//       dummy += cells[i].mu0_cal;
//       k++;
//     }
//     dummy = dummy/k;
//     *mean_mu_s = dummy;
// }



void CSR_Create(int M,int N, int nz, int *I, int *J, double *val,int *AI_CSR)
{
// M: number of rows
// N: number of columns
// nz: number of nonzero elements
// I, J, val are inputs in COO
// AI_CSR, RI, RJ, Rval are outputs in CSR
int i,j,k,temp;
int *RI = (int *) malloc( nz * sizeof(int) );
int *RJ = (int *) malloc( nz * sizeof(int) );
double *Rval = (double *) malloc( nz * sizeof(double) );
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
