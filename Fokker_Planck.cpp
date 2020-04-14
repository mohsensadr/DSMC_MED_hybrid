#include "cpp_headers.h"
/*
#include <math.h>
#include <random>
#include <iostream>
// #include <cmath>
#include <string>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>


#include <omp.h>
#include "omp.h"
#include <time.h>
*/

void random_generator(double *xi_1, double *xi_2, double *xi_3, double *xi_1f, double *xi_2f, double *xi_3f, struct CELLS *cells, struct BOX *box,  struct GAS *gas, int *index)
{
if ( (*gas).model == "FP_dense" || (*gas).model == "FP_ideal_gas" || (*gas).model == "FP_linear"){
	double start,end;
	//start = omp_get_wtime();

	int i, num, id;
	int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
	int num_inside;
	double sum_1, sum_2, sum_3;
	double std = 1.0;
	double mean = 0.0;

/*
// Testing
	double m1[3],m2[3],s1[3],s2[3];
	double st1, st2, f1, f2;
// testing NRVs
	st1 = omp_get_wtime();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<> dis1(mean,std);
	std::normal_distribution<> dis2(mean,std);
	std::normal_distribution<> dis3(mean,std);
	for(id=0; id<(*gas).N; id++){
	      	xi_1f[id] = 1.0*dis1(gen);
	      	xi_2f[id] = 1.0*dis2(gen);
	      	xi_3f[id] = 1.0*dis3(gen);
	}
	f1 = omp_get_wtime();
	s1[0] = standard_deviation(xi_1f, (*gas).N, &m1[0]);
	s1[1] = standard_deviation(xi_1f, (*gas).N, &m1[1]);
	s1[2] = standard_deviation(xi_1f, (*gas).N, &m1[2]);
// testing NRVs
	// Box–Muller transform
	// https://en.wikipedia.org/wiki/Box–Muller_transform
	st2 = omp_get_wtime();
	srand((*box).step+time(NULL));
	for(id=0; id<(*gas).N; id++)
	      	xi_1f[id] = generateGaussianNoise(mean, std);
	for(id=0; id<(*gas).N; id++)
	      	xi_2f[id] = generateGaussianNoise(mean, std);
	for(id=0; id<(*gas).N; id++)
	      	xi_3f[id] = generateGaussianNoise(mean, std);
	f2 = omp_get_wtime();
	s2[0] = standard_deviation(xi_1f, (*gas).N, &m2[0]);
	s2[1] = standard_deviation(xi_2f, (*gas).N, &m2[1]);
	s2[2] = standard_deviation(xi_3f, (*gas).N, &m2[2]);
	for(i=0; i<3; i++){
		printf("m1[%d] = %e and |s1[%d]-1| = %e and took %e\n", i, m1[i], i, fabs(s1[i]-1.0), f1-st1);
		printf("m2[%d] = %e and |s2[%d]-1| = %e and took %e\n", i, m2[i], i, fabs(s2[i]-1.0), f2-st2);
	}
	for(i=0; i<3; i++)
		printf("xi_1f[%i] = %e,  ",i, xi_1f[i]);
	printf("\n");
	for(i=0; i<3; i++)
		printf("xi_2f[%i] = %e,  ",i, xi_2f[i]);
	printf("\n");
	for(i=0; i<3; i++)
		printf("xi_3f[%i] = %e,  ",i, xi_3f[i]);
	printf("\n");
*/


	srand((*box).step+time(NULL));
	for(id=0; id<(*gas).N; id++)
	      	xi_1f[id] = generateGaussianNoise(mean, std);
	for(id=0; id<(*gas).N; id++)
	      	xi_2f[id] = generateGaussianNoise(mean, std);
	for(id=0; id<(*gas).N; id++)
	      	xi_3f[id] = generateGaussianNoise(mean, std);



	double *sum_11 = (double *) malloc(  num_cell * sizeof(double) );
	double *sum_22 = (double *) malloc(  num_cell * sizeof(double) );
	double *sum_33 = (double *) malloc(  num_cell * sizeof(double) );

	for(num=0; num<num_cell; num++){
		sum_11[num] = 0.0;
		sum_22[num] = 0.0;
		sum_33[num] = 0.0;
	}
	for(i=0; i<(*gas).N; i++){
		num 		= index[i];
		sum_11[num]    += xi_1f[i];
		sum_22[num]    += xi_2f[i];
		sum_33[num]    += xi_3f[i];
  	}
  	for(i=0; i<(*gas).N; i++){
		num 		= index[i];
		xi_1[i] = xi_1f[i] - sum_11[num]/(1.0*cells[num].num_inside);
		xi_2[i] = xi_2f[i] - sum_22[num]/(1.0*cells[num].num_inside);
		xi_3[i] = xi_3f[i] - sum_33[num]/(1.0*cells[num].num_inside);
  	}
  	for(num=0; num<num_cell; num++){
		sum_11[num] = 0.0;
		sum_22[num] = 0.0;
		sum_33[num] = 0.0;
  	}
  	for(i=0; i<(*gas).N; i++){
		num 		= index[i];
		sum_11[num]    += xi_1[i]*xi_1[i];
		sum_22[num]    += xi_2[i]*xi_2[i];
		sum_33[num]    += xi_3[i]*xi_3[i];
  	}
  	for(num=0; num<num_cell; num++){
		sum_11[num] = sqrt(sum_11[num]/(1.0*cells[num].num_inside));
		sum_22[num] = sqrt(sum_22[num]/(1.0*cells[num].num_inside));
		sum_33[num] = sqrt(sum_33[num]/(1.0*cells[num].num_inside));
	}
	for(i=0; i<(*gas).N; i++){
		num 		= index[i];
		if(cells[num].num_inside>1){
			xi_1f[i] = xi_1[i]/sum_11[num];
			xi_2f[i] = xi_2[i]/sum_22[num];
			xi_3f[i] = xi_3[i]/sum_33[num];
		}
	}
	free(sum_11); free(sum_22); free(sum_33);
	//end = omp_get_wtime();

	// Print time it took
	//if ( (*box).step % (*box).every == 0)
	//	printf("Random generator took %e s\n", end - start);
}

/*
//#pragma omp parallel for  private(num_inside, i, sum_1, sum_2, sum_3, id)
for(num=0; num<num_cell; num++){
  num_inside = cells[num].num_inside;
  double *xi_1 = (double *) malloc(  num_inside * sizeof(double) );
  double *xi_2 = (double *) malloc(  num_inside * sizeof(double) );
  double *xi_3 = (double *) malloc(  num_inside * sizeof(double) );
  double *xi_telda_1 = (double *) malloc(  num_inside * sizeof(double) );
  double *xi_telda_2 = (double *) malloc(  num_inside * sizeof(double) );
  double *xi_telda_3 = (double *) malloc(  num_inside * sizeof(double) );
  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
      	xi_1[i] = 1.0*dis1(gen);
      	xi_2[i] = 1.0*dis2(gen);
      	xi_3[i] = 1.0*dis3(gen);
  }
  sum_1 = 0.0;
  sum_2 = 0.0;
  sum_3 = 0.0;
  for(i=0; i<num_inside; i++){
      	sum_1 += xi_1[i];
      	sum_2 += xi_2[i];
      	sum_3 += xi_3[i];
  }
  for(i=0; i<num_inside; i++){
	      xi_telda_1[i] = xi_1[i] - sum_1/(1.0*num_inside);
      	xi_telda_2[i] = xi_2[i] - sum_2/(1.0*num_inside);
      	xi_telda_3[i] = xi_3[i] - sum_3/(1.0*num_inside);
  }
  sum_1 = 0.0;
  sum_2 = 0.0;
  sum_3 = 0.0;
  for(i=0; i<num_inside; i++){
	      sum_1 += xi_telda_1[i]*xi_telda_1[i];
      	sum_2 += xi_telda_2[i]*xi_telda_2[i];
      	sum_3 += xi_telda_3[i]*xi_telda_3[i];
  }
  sum_1 = sqrt(sum_1/(1.0*num_inside));
  sum_2 = sqrt(sum_2/(1.0*num_inside));
  sum_3 = sqrt(sum_3/(1.0*num_inside));

  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[ i ];
      	xi_1f[id] = xi_telda_1[i]/sum_1;
      	xi_2f[id] = xi_telda_2[i]/sum_2;
      	xi_3f[id] = xi_telda_3[i]/sum_3;
  }
  free(xi_1);free(xi_2);free(xi_3);
  free(xi_telda_1); free(xi_telda_2); free(xi_telda_3);
}
*/

}
void alpha_calculation( double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells, struct BOX *box){

  int n = 12;

//   double pi = acos(-1);
    double lambda_crr;
    double w, c_v;
  double mu_crr ;//= (*gas).mu_corr;
  double p_crr ;//= (*gas).n*(*gas).b;
  double k_crr,Y;
  int i,j,k, num;
  double a,b,c,d;
  double div_U = 0.0;
//  double MM[3*3], deviatoric[3*3];
  double *deviatoric = (double *) malloc( 6 * sizeof(double) );
  double *grad_U = (double *) malloc( 6 * sizeof(double) );
  double *grad_T = (double *) malloc( 3 * sizeof(double) );

  double *aa = (double *) malloc( n*n * sizeof(double) );

  double trace = 0.0;
  int ii,jj, kk;
  double pi = acos(-1.0);
 // double A[36];
      double **A;//*a;
    for( i=0; i<6; i++)
      grad_U[i] = 0.0;
    for( i=0; i<3; i++)
      grad_T[i] = 0.0;
    A=new double *[12];
    for(int i=0;i<12;i++)
    {
        A[i]=new double [12];
    }

    double *rhs = (double *) malloc( 12 * sizeof(double) );
    double *x = (double *) malloc( 12 * sizeof(double) );
//   double alpha[6];
//   double *x;
//   x = new double [6];
//int nrhs;
//int lda;
//int ipiv[n];
//int ldb;
//int info;
double coeff_mu;
double coeff_k;
int doit;
    for(j=0; j<(*box).N[1]; j++){
      for(i=0; i<(*box).N[0]; i++){
	num =  j*(*box).N[0] + i;
	doit = 1;
	//if( (*box).problem==inverted && (cells[num].cell_center[1]< (*box).ghost || cells[num].cell_center[1] > (*box).Len[1]-(*box).ghost) )
	//	doit = 0;
	if(cells[num].num_inside < 2)
		doit = 0;
	if(doit ==0){
		for(ii=0; ii<9; ii++)
			cells[num].alpha[ii] = 0.0;

		 cells[num].beta[0] = 0.0;
		 cells[num].beta[1] = 0.0;
		 cells[num].beta[2] = 0.0;
	}
	else{
	trace = cells[num].PIJ[0]+cells[num].PIJ[3]+cells[num].PIJ[5];

	for(ii=0; ii<6; ii++)
	    deviatoric[ii] = cells[num].PIJ[ii];
	deviatoric[0] =  deviatoric[0] - (1.0/3.0)*trace;
	deviatoric[3] =  deviatoric[3] - (1.0/3.0)*trace;
	deviatoric[5] =  deviatoric[5] - (1.0/3.0)*trace;

        a = cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5];// 111 + 221 + 331
	b = cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8];// 112 + 222 + 332
	c = cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9];// 113 + 223 + 333
	d = cells[num].M4[0] + cells[num].M4[3] + cells[num].M4[5] - cells[num].DM2*cells[num].DM2 ;
//			     + cells[num].U_space[0]*cells[num].Q[0] + cells[num].U_space[1]*cells[num].Q[1] + cells[num].U_space[2]*cells[num].Q[2]; // 11 + 22 + 33

	  for(ii=0; ii<n; ii++)
	    for(jj=0; jj<n; jj++){
	       A[ii][jj] = 0.0;
             aa[ii*n+jj] = 0.0;
         }

	 //A
	 A[0][0] = cells[num].PIJ[0];
	 A[0][1] = cells[num].PIJ[1];
	 A[0][2] = cells[num].PIJ[2];

	 A[0][9] = a;// 111 + 221 + 331

	 //B
	 A[1][0] = cells[num].PIJ[1];
	 A[1][1] = cells[num].PIJ[3];
	 A[1][2] = cells[num].PIJ[4];

	 A[1][9] = a;
 	 //A[1][9] = cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8];// 112 + 222 + 332

	 //C
	 A[2][0] = cells[num].PIJ[2];
	 A[2][1] = cells[num].PIJ[4];
	 A[2][2] = cells[num].PIJ[5];

	 A[2][9] = a;
  	 //A[2][9] = cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9];// 113 + 223 + 333

	 //D
	 A[3][3] =  cells[num].PIJ[0];
	 A[3][4] =  cells[num].PIJ[1];
	 A[3][5] =  cells[num].PIJ[2];

	 A[3][10] = b;
	// A[3][9] =  -( cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8] );// 112 + 222 + 332
	// A[3][10] =  cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5];    // 111 + 221 + 331

	 //E
	 A[4][3] = cells[num].PIJ[1];
	 A[4][4] = cells[num].PIJ[3];
	 A[4][5] = cells[num].PIJ[4];

	 A[4][10] = b;
	 //A[4][10] =  cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8];    // 112 + 222 + 332

	 //F
	 A[5][3] = cells[num].PIJ[2];
	 A[5][4] = cells[num].PIJ[4];
	 A[5][5] = cells[num].PIJ[5];

	 A[5][10] = b;
	 //A[5][10] = cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9];    // 113 + 223 + 333

	 //G
	 //A[6][0] = -cells[num].PIJ[2];
	 //A[6][1] = -cells[num].PIJ[4];
	 //A[6][2] = -cells[num].PIJ[5];
	 A[6][6] =  cells[num].PIJ[0];
	 A[6][7] =  cells[num].PIJ[1];
	 A[6][8] =  cells[num].PIJ[2];

	 A[6][11] = c;
	 //A[6][9]  =  -( cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9] );// 113 + 223 + 333
	 //A[6][11] =     cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5];    // 111 + 221 + 331

	 //H
	 //A[7][3] = -cells[num].PIJ[2];
	 //A[7][4] = -cells[num].PIJ[4];
	 //A[7][5] = -cells[num].PIJ[5];
	 A[7][6] =  cells[num].PIJ[1];
	 A[7][7] =  cells[num].PIJ[3];
	 A[7][8] =  cells[num].PIJ[4];

	 A[7][11] = c;
	 //A[7][10]  =  -( cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9] );// 113 + 223 + 333
	 //A[7][11] =     cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8];    // 112 + 222 + 332

	 //I
	 A[8][6] = cells[num].PIJ[2];
	 A[8][7] = cells[num].PIJ[4];
	 A[8][8] = cells[num].PIJ[5];

	 A[8][11] = c;
	 //A[8][11] = cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9];    // 113 + 223 + 333

	 //J
	 A[9][0] = a;
	 A[9][1] = b;
	 A[9][2] = c;
	 A[9][9] = d;

	//K
	 A[10][3] = a;
	 A[10][4] = b;
	 A[10][5] = c;
	 A[10][10] = d;

	 //L
	 A[11][6] = a;
	 A[11][7] = b;
	 A[11][8] = c;
	 A[11][11] = d;

//          mu_crr = 4.0*cells[num].n*(*gas).b/5.0;
          Y = increase_collision_rate(cells[num].n, gas);

          //mu_crr = pow( 1.0+2.0*cells[num].n*(*gas).b*Y/5.0 , 2.0 )/Y;
	  p_crr =  ( 1.0 + cells[num].n*(*gas).b*Y );
	  //k_crr =  pow( 1.0+3.0*cells[num].n*(*gas).b*Y/5.0, 2.0 )/Y;
         mu_crr = ( 1.0+2.0*cells[num].n*(*gas).b*Y/5.0 );
	 k_crr =  ( 1.0+3.0*cells[num].n*(*gas).b*Y/5.0 );


/*
//// no bulk term
	 rhs[0] = (mu_crr - 1.0)*deviatoric[0] + (p_crr-1.0)*trace/3.0;
	 rhs[1] = (mu_crr - 1.0)*deviatoric[1];
	 rhs[2] = (mu_crr - 1.0)*deviatoric[2];
	 rhs[3] = (mu_crr - 1.0)*deviatoric[1];
	 rhs[4] = (mu_crr - 1.0)*deviatoric[3] + (p_crr-1.0)*trace/3.0;
	 rhs[5] = (mu_crr - 1.0)*deviatoric[4];
	 rhs[6] = (mu_crr - 1.0)*deviatoric[2];
	 rhs[7] = (mu_crr - 1.0)*deviatoric[4];
	 rhs[8] = (mu_crr - 1.0)*deviatoric[5] + (p_crr-1.0)*trace/3.0;

// we need the 2 multiplication, because we have divided both sides of energy equation by rho/2
	 rhs[9]  = (k_crr-1.0)*a;
	 rhs[10] = (k_crr-1.0)*b;
	 rhs[11] = (k_crr-1.0)*c;
*/


//// with bulk term
    w = (4.0/9.0)*sqrt(pi)*Y*pow(cells[num].n,2.0)*pow((*gas).sigma,4.0)*sqrt((*gas).m*(*gas).kb*cells[num].T);
	//	w = 0.0;
   c_v = (3*(*gas).kb)/(2*(*gas).m);

	 //w = 0.0;
	 //cells[num].gamma_st = 0.0;
    for( ii=0; ii<6; ii++)
      grad_U[ii] = 0.0;

	/*
   if( fabs(cells[num].M_n[0])>1e-15 &&  fabs(cells[num].M_n[1])>1e-15 && fabs(cells[num].M_n[2])>1e-15 &&  fabs(cells[num].M_n[3])>1e-15){
     grad_U[0] = (cells[num].M1_f[1]/cells[num].M_n[1]-cells[num].M1_f[0]/cells[num].M_n[0])/(cells[num].dim[1]-cells[num].dim[0]) ;     // du/dx
     grad_U[1] = 0.5*(cells[num].M1_f[3]/cells[num].M_n[3]-cells[num].M1_f[2]/cells[num].M_n[2])/(cells[num].dim[3]-cells[num].dim[2])
                +0.5*(cells[num].M2_f[1]/cells[num].M_n[1]-cells[num].M2_f[0]/cells[num].M_n[0])/(cells[num].dim[1]-cells[num].dim[0]) ; // 1/2*(du/dy+dv/dx)
     grad_U[2] = 0.0
                +0.5*(cells[num].M3_f[1]/cells[num].M_n[1]-cells[num].M3_f[0]/cells[num].M_n[0])/(cells[num].dim[1]-cells[num].dim[0]) ; // 1/2*(du/dz+dw/dx)
     grad_U[3] = (cells[num].M2_f[3]/cells[num].M_n[3]-cells[num].M2_f[2]/cells[num].M_n[2])/(cells[num].dim[3]-cells[num].dim[2]);     // dv/dy
     grad_U[4] = 0.5*(cells[num].M3_f[3]/cells[num].M_n[3]-cells[num].M3_f[2]/cells[num].M_n[2])/(cells[num].dim[3]-cells[num].dim[2]); // 1/2*(dv/dz+dw/dy)
     grad_U[5] = 0.0;     // dw/dz
  }
//   div_U =  (cells[num].M2_f[3]-cells[num].M2_f[2])/(cells[num].dim[3]-cells[num].dim[2]);
   div_U = grad_U[0] + grad_U[3] + grad_U[5];
   grad_U[0] = grad_U[0] - div_U/3.0;
   grad_U[3] = grad_U[3] - div_U/3.0;
   grad_U[5] = grad_U[5] - div_U/3.0;
*/

	//div_U = 0.0;
	// following is correct. Just for testing purposes with Marcel is commented out

	 if(num==0){
		 div_U = (cells[num+1].U_space[1]-cells[num].U_space[1])/(*box).delta_dim[1];
	 }
	 else if(num==(*box).N[1]-1){
		 div_U = (cells[num].U_space[1]-cells[num-1].U_space[1])/(*box).delta_dim[1];
	 }
	 else{
		 div_U = (cells[num+1].U_space[1]-cells[num-1].U_space[1])/(*box).delta_dim[1]/2.0;
	 }


   for( ii=0; ii<3; ii++)
     grad_T[ii] = 0.0;
	cells[num].divU = div_U;
  /*
   if( fabs(cells[num].M_n[0])>1e-15 &&  fabs(cells[num].M_n[1])>1e-15)
      grad_T[0] = ((*gas).m/(3.0*(*gas).kb))*(cells[num].T_f[1]/cells[num].M_n[1]-cells[num].T_f[0]/cells[num].M_n[0])/(cells[num].dim[1]-cells[num].dim[0]);
   if( fabs(cells[num].M_n[2])>1e-15 &&  fabs(cells[num].M_n[3])>1e-15)
       grad_T[1] = ((*gas).m/(3.0*(*gas).kb))*(cells[num].T_f[3]/cells[num].M_n[3]-cells[num].T_f[2]/cells[num].M_n[2])/(cells[num].dim[3]-cells[num].dim[2]);
   grad_T[2] = 0.0;
   */
   coeff_mu  = 3.0*1.002/5.0*pow(cells[num].n*(*gas).b*Y,2.0)/(1+2.0*cells[num].n*(*gas).b*Y/5.0);
	 coeff_k = 1.002*pow(cells[num].n*(*gas).b*Y,2.0)/(1+3.0*cells[num].n*(*gas).b*Y/5.0)/5.0;

	 //div_U=0.0; cells[num].gamma_st=0.0; coeff_mu = 0.0; coeff_k = 0.0;
	 //cells[num].gamma_st=0.0; coeff_k = 0.0; k_crr = 1.0;

	 rhs[0] = (mu_crr - 1.0)*deviatoric[0] + (p_crr-1.0)*trace/3.0
		 - w*div_U/( (*gas).m*cells[num].n ) + coeff_mu*deviatoric[0]
		 - cells[num].gamma_st*cells[num].M4[0];
	 rhs[1] = (mu_crr - 1.0)*deviatoric[1] + coeff_mu*deviatoric[1]
		 - cells[num].gamma_st*cells[num].M4[1];
	 rhs[2] = (mu_crr - 1.0)*deviatoric[2] + coeff_mu*deviatoric[2]
		 - cells[num].gamma_st*cells[num].M4[2];
	 rhs[3] = (mu_crr - 1.0)*deviatoric[1] + coeff_mu*deviatoric[1]
		 - cells[num].gamma_st*cells[num].M4[1];
	 rhs[4] = (mu_crr - 1.0)*deviatoric[3] + (p_crr-1.0)*trace/3.0
		- w*div_U/( (*gas).m*cells[num].n ) + coeff_mu*deviatoric[3]
		 - cells[num].gamma_st*cells[num].M4[3];
	 rhs[5] = (mu_crr - 1.0)*deviatoric[4] + coeff_mu*deviatoric[4]
		 - cells[num].gamma_st*cells[num].M4[4];
	 rhs[6] = (mu_crr - 1.0)*deviatoric[2] + coeff_mu*deviatoric[2]
		 - cells[num].gamma_st*cells[num].M4[2];
	 rhs[7] = (mu_crr - 1.0)*deviatoric[4] + coeff_mu*deviatoric[4]
		 - cells[num].gamma_st*cells[num].M4[4];
	 rhs[8] = (mu_crr - 1.0)*deviatoric[5] + (p_crr-1.0)*trace/3.0
		 - w*div_U/( (*gas).m*cells[num].n )+ coeff_mu*deviatoric[5]
		 - cells[num].gamma_st*cells[num].M4[5];

// we need the 2 multiplication, because we have divided both sides of energy equation by rho/2
	 rhs[9]  = (k_crr-1.0)*a + coeff_k*a
		 - cells[num].gamma_st*(cells[num].M5[0])// + 2.0*cells[num].M4[0]*cells[num].U_space[0]+2.0*cells[num].M4[1]*cells[num].U_space[1]+2.0*cells[num].M4[2]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5] ) );
	 rhs[10] = (k_crr-1.0)*b + coeff_k*b
		 - cells[num].gamma_st*(cells[num].M5[1])// + 2.0*cells[num].M4[1]*cells[num].U_space[0]+2.0*cells[num].M4[3]*cells[num].U_space[1]+2.0*cells[num].M4[4]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8] ) );
	 rhs[11] = (k_crr-1.0)*c + coeff_k*c
		 - cells[num].gamma_st*(cells[num].M5[2])// + 2.0*cells[num].M4[2]*cells[num].U_space[0]+2.0*cells[num].M4[4]*cells[num].U_space[1]+2.0*cells[num].M4[5]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9] ) );

   /*
	 rhs[0] = (mu_crr - 1.0)*deviatoric[0] + (p_crr-1.0)*trace/3.0
		 - w*div_U/( (*gas).m*cells[num].n )- (6.0*w/5.0)*grad_U[0]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[0];
	 rhs[1] = (mu_crr - 1.0)*deviatoric[1] - (6.0*w/5.0)*grad_U[1]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[1];
	 rhs[2] = (mu_crr - 1.0)*deviatoric[2] - (6.0*w/5.0)*grad_U[2]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[2];
	 rhs[3] = (mu_crr - 1.0)*deviatoric[1] - (6.0*w/5.0)*grad_U[1]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[1];
	 rhs[4] = (mu_crr - 1.0)*deviatoric[3] + (p_crr-1.0)*trace/3.0
		- w*div_U/( (*gas).m*cells[num].n ) - (6.0*w/5.0)*grad_U[3]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[3];
	 rhs[5] = (mu_crr - 1.0)*deviatoric[4] - (6.0*w/5.0)*grad_U[4]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[4];
	 rhs[6] = (mu_crr - 1.0)*deviatoric[2] - (6.0*w/5.0)*grad_U[2]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[2];
	 rhs[7] = (mu_crr - 1.0)*deviatoric[4] - (6.0*w/5.0)*grad_U[4]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[4];
	 rhs[8] = (mu_crr - 1.0)*deviatoric[5] + (p_crr-1.0)*trace/3.0
		 - w*div_U/( (*gas).m*cells[num].n )- (6.0*w/5.0)*grad_U[5]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*cells[num].M4[5];

// we need the 2 multiplication, because we have divided both sides of energy equation by rho/2
	 rhs[9]  = (k_crr-1.0)*a - 2.0*c_v*w*grad_T[0]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*(cells[num].M5[0])// + 2.0*cells[num].M4[0]*cells[num].U_space[0]+2.0*cells[num].M4[1]*cells[num].U_space[1]+2.0*cells[num].M4[2]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5] ) );
	 rhs[10] = (k_crr-1.0)*b - 2.0*c_v*w*grad_T[1]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*(cells[num].M5[1])// + 2.0*cells[num].M4[1]*cells[num].U_space[0]+2.0*cells[num].M4[3]*cells[num].U_space[1]+2.0*cells[num].M4[4]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8] ) );
	 rhs[11] = (k_crr-1.0)*c - 2.0*c_v*w*grad_T[2]/( (*gas).m*cells[num].n )
		 - cells[num].gamma_st*(cells[num].M5[2])// + 2.0*cells[num].M4[2]*cells[num].U_space[0]+2.0*cells[num].M4[4]*cells[num].U_space[1]+2.0*cells[num].M4[5]*cells[num].U_space[2])
		 + cells[num].gamma_st*( cells[num].DM2*( cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9] ) );
  */
  // for(ii=0; ii<n; ii++)
  //   for(jj=0; jj<n; jj++){
  //     aa[ii*n+jj] = A[ii][jj];
  // }

/*
   nrhs = 1;
   lda=n;
   ipiv[n];
   ldb = n;
   info;

   dgesv(&n, &nrhs, aa, &lda, ipiv, rhs, &ldb, &info);

   //printf("info = %d\n", info);
   for(ii=0; ii<n; ii++)
     x[ii] = rhs[ii];
*/
   Gauss_Jordan(A, rhs, n, x);

   for(ii=0; ii<n; ii++){
     if(x[ii] != x[ii]){
        printf("Here, alpha_calculation!!, N = %d\n",cells[num].num_inside);
	printf("cells[%d].n = %e\n",num, cells[num].n);
	printf("cells[%d].weight = %lf\n", num, cells[num].weight) ;
	printf("cells[%d].indices_inside[0] = %d\n",num, cells[num].indices_inside[0]);
	for(jj=0; jj<9; jj++)
		printf("cell[%d].M3[%d] = %lf\n",num, jj, cells[num].M3[jj]);
	for(jj=0; jj<6; jj++)
		printf("cell[%d].PIJ[%d] = %lf\n",num,jj,cells[num].PIJ[jj]);
	for(jj=0; jj<9; jj++)
		printf("rhs[%d] = %lf\n",jj, rhs[jj]);
        exit(1);
     }
   }

   for(ii=0; ii<9; ii++)
	   cells[num].alpha[ii] = x[ii];
	 cells[num].beta[0] = x[9];
	 cells[num].beta[1] = x[10];
	 cells[num].beta[2] = x[11];
    }
      }
    }
   free(deviatoric); free(rhs); free(x);
   free(aa); free(grad_U), free(grad_T);

    for(ii=0;ii<n;ii++)
        delete [] A[ii];
   delete  [] A;
}


void c_gamma_calculation( double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cell, struct BOX *box){

	int n = 9;

  //int nrhs = 1;
  //int lda=n;
  //int ipiv[n];
  //int ldb = n;
  //int info;
  //double Lambda;

  int i,j,k, num,ii,jj, kk;
  double **A;//, **A0;//*a;
  double mu;
  double Y;
  double mu_crr, p_crr;
  double Pr;
  	 //double nubol;//= ( cell[num].n*(*gas).kb*cell[num].T )/mu;
	 double tau;
    A=new double *[9];
    //A0=new double *[9];
    for(int i=0;i<9;i++)
    {
        A[i]=new double [9];
	//A0[i]=new double [9];
    }

    double *rhs = (double *) malloc( 9 * sizeof(double) );
    //double *rhs0 = (double *) malloc( 9 * sizeof(double) );
    double *x = (double *) malloc( 9 * sizeof(double) );
    int doit = 0;
   for(k=0; k<(*box).N[2]; k++){
    for(j=0; j<(*box).N[1]; j++){
      for(i=0; i<(*box).N[0]; i++){
	num = k*(*box).N[0]*(*box).N[1] + j*(*box).N[0] + i;
	doit = 1;
	//if( (*box).problem==inverted && (cell[num].cell_center[1]< (*box).ghost || cell[num].cell_center[1] > (*box).Len[1]-(*box).ghost) )
	//	doit = 0;
	if(cell[num].num_inside < 2)
		doit = 0;
	if(doit ==0){
		for(ii=0; ii<6; ii++)
			cell[num].c[ii] = 0.0;

		 cell[num].gamma[0] = 0.0;
		 cell[num].gamma[1] = 0.0;
		 cell[num].gamma[2] = 0.0;
	}
	else{
		for(ii=0; ii<9; ii++)
	  		for(jj=0; jj<9; jj++)
	     			A[ii][jj] = 0.0;
		for(ii=0; ii<9; ii++)
			rhs[ii] = 0.0;
		for(ii=0; ii<9; ii++)
			x[ii] = 0.0;
	 A[0][0] = 2.0*cell[num].PIJ[0];//11
	 A[0][1] = 2.0*cell[num].PIJ[1];//12
	 A[0][2] = 2.0*cell[num].PIJ[2];//13
	 A[0][6] = 2.0*cell[num].Q[0];//1

	 A[1][0] = cell[num].PIJ[1];//12
	 A[1][1] = cell[num].PIJ[0] + cell[num].PIJ[3];//11+22
	 A[1][2] = cell[num].PIJ[4];//23
	 A[1][3] = cell[num].PIJ[1];//12
	 A[1][4] = cell[num].PIJ[2];//13
	 A[1][6] = cell[num].Q[1];//2
	 A[1][7] = cell[num].Q[0];//1

	 A[2][0] = cell[num].PIJ[2];//13
	 A[2][1] = cell[num].PIJ[4];//23
	 A[2][2] = cell[num].PIJ[0] + cell[num].PIJ[5];//11+33
	 A[2][4] = cell[num].PIJ[1];//12
	 A[2][5] = cell[num].PIJ[2];//13
	 A[2][6] = cell[num].Q[2];//3
	 A[2][8] = cell[num].Q[0];//1

	 A[3][1] = 2.0*cell[num].PIJ[1];//12
	 A[3][3] = 2.0*cell[num].PIJ[3];//22
	 A[3][4] = 2.0*cell[num].PIJ[4];//23
	 A[3][7] = 2.0*cell[num].Q[1];//2

	 A[4][1] = cell[num].PIJ[2];//13
	 A[4][2] = cell[num].PIJ[1];//12
	 A[4][3] = cell[num].PIJ[4];//23
	 A[4][4] = cell[num].PIJ[3]+cell[num].PIJ[5];//22+33
	 A[4][5] = cell[num].PIJ[4];//23
	 A[4][7] = cell[num].Q[2];//3
	 A[4][8] = cell[num].Q[1];//2

	 A[5][2] = 2.0*cell[num].PIJ[2];//13
	 A[5][4] = 2.0*cell[num].PIJ[4];//23
	 A[5][5] = 2.0*cell[num].PIJ[5];//33
	 A[5][8] = 2.0*cell[num].Q[2];//3

/////////////////////////////////////////////////////////////////////////////////////////

	 A[6][0] = cell[num].Q[0] + 2.0*cell[num].M3[0];//1+111
	 A[6][1] = cell[num].Q[1] + 4.0*cell[num].M3[1];//2+112
	 A[6][2] = cell[num].Q[2] + 4.0*cell[num].M3[2];//3+113
	 A[6][3] = 2.0*cell[num].M3[3];//122
	 A[6][4] = 4.0*cell[num].M3[4];//123
	 A[6][5] = 2.0*cell[num].M3[5];//133
	 A[6][6] = 2.0*cell[num].M4[0]-2.0*cell[num].PIJ[0]*cell[num].DM2;//11-11
	 A[6][7] = 2.0*cell[num].M4[1]-2.0*cell[num].PIJ[1]*cell[num].DM2;//12-12
	 A[6][8] = 2.0*cell[num].M4[2]-2.0*cell[num].PIJ[2]*cell[num].DM2;//13-13

	 A[7][0] = 2.0*cell[num].M3[1];//112
	 A[7][1] = cell[num].Q[0] + 4.0*cell[num].M3[3];//1+122
	 A[7][2] = 4.0*cell[num].M3[4];//123
	 A[7][3] = cell[num].Q[1] + 2.0*cell[num].M3[6];//2+222
	 A[7][4] = cell[num].Q[2] + 4.0*cell[num].M3[7];//3+223
	 A[7][5] = 2.0*cell[num].M3[8];//233
	 A[7][6] = 2.0*cell[num].M4[1]-2.0*cell[num].PIJ[1]*cell[num].DM2;//12-12
	 A[7][7] = 2.0*cell[num].M4[3]-2.0*cell[num].PIJ[3]*cell[num].DM2;//22-22
	 A[7][8] = 2.0*cell[num].M4[4]-2.0*cell[num].PIJ[4]*cell[num].DM2;//23-23

	 A[8][0] = 2.0*cell[num].M3[2];//113
	 A[8][1] = 4.0*cell[num].M3[4];//123
	 A[8][2] = cell[num].Q[0] + 4.0*cell[num].M3[5];//1+133
	 A[8][3] = 2.0*cell[num].M3[7];//223
	 A[8][4] = cell[num].Q[1] + 4.0*cell[num].M3[8];//2+233
	 A[8][5] = cell[num].Q[2] + 2.0*cell[num].M3[9];//3+333
	 A[8][6] = 2.0*cell[num].M4[2]-2.0*cell[num].PIJ[2]*cell[num].DM2;//13-13
	 A[8][7] = 2.0*cell[num].M4[4]-2.0*cell[num].PIJ[4]*cell[num].DM2;//23-23
	 A[8][8] = 2.0*cell[num].M4[5]-2.0*cell[num].PIJ[5]*cell[num].DM2;//33-33

	 A[6][6] += cell[num].DM4-(cell[num].DM2)*(cell[num].DM2);
	 A[7][7] += cell[num].DM4-(cell[num].DM2)*(cell[num].DM2);
	 A[8][8] += cell[num].DM4-(cell[num].DM2)*(cell[num].DM2);

	//for(jj=0; jj<9; jj++)
	//	for(kk=0; kk<9; kk++)
	//		A0[jj][kk] = A[jj][kk];

         Pr = 2.0/3.0;
         if( (*gas).model == "FP_dense"){
                 Y = increase_collision_rate(cell[num].n, gas);
        }
        else{
		Y = 1.0;
        }
        //Y = 1.0;
	scale_mu(cell[num].T, (*gas).T0, &mu, (*gas).mu, gas);
	tau = 2.0*mu/(cell[num].n*(*gas).kb*cell[num].T);

	rhs[0] = -2.0*cell[num].Lambda*cell[num].M4[0];
	rhs[1] = -2.0*cell[num].Lambda*cell[num].M4[1];
	rhs[2] = -2.0*cell[num].Lambda*cell[num].M4[2];
	rhs[3] = -2.0*cell[num].Lambda*cell[num].M4[3];
	rhs[4] = -2.0*cell[num].Lambda*cell[num].M4[4];
	rhs[5] = -2.0*cell[num].Lambda*cell[num].M4[5];



	 rhs[6] = Y*( 3.0/tau - 2.0*Pr/tau )*cell[num].Q[0];
	 rhs[7] = Y*( 3.0/tau - 2.0*Pr/tau )*cell[num].Q[1];
	 rhs[8] = Y*( 3.0/tau - 2.0*Pr/tau )*cell[num].Q[2];

       rhs[6] += cell[num].Lambda*( -3.0*cell[num].M5[0]+cell[num].DM2*cell[num].Q[0]+2.0*(cell[num].PIJ[0]*cell[num].Q[0]+cell[num].PIJ[1]*cell[num].Q[1]+cell[num].PIJ[2]*cell[num].Q[2]) );
	 rhs[7] += cell[num].Lambda*( -3.0*cell[num].M5[1]+cell[num].DM2*cell[num].Q[1]+2.0*(cell[num].PIJ[1]*cell[num].Q[0]+cell[num].PIJ[3]*cell[num].Q[1]+cell[num].PIJ[4]*cell[num].Q[2]) );
	 rhs[8] += cell[num].Lambda*( -3.0*cell[num].M5[2]+cell[num].DM2*cell[num].Q[2]+2.0*(cell[num].PIJ[2]*cell[num].Q[0]+cell[num].PIJ[4]*cell[num].Q[1]+cell[num].PIJ[5]*cell[num].Q[2]) );

/*
   for(ii=0; ii<n; ii++)
     for(jj=0; jj<n; jj++)
       aa[ii*n+jj] = A[ii][jj];
   dgesv(&n, &nrhs, aa, &lda, ipiv, rhs, &ldb, &info);
   for(ii=0; ii<n; ii++)
     x[ii] = rhs[ii];
*/

		//for(ii=0; ii<9; ii++)
		//	rhs0[ii] = rhs[ii];
   Gauss_Jordan(A, rhs, n, x);
   for(ii=0; ii<n; ii++){
     if(x[ii] != x[ii]){
        printf("Here, we have a nan, c_gamma!!, N = %d\n",cell[num].num_inside);
	printf("cells[%d].n = %e\n",num, cell[num].n);
	printf("cells[%d].weight = %lf\n", num, cell[num].weight) ;
	printf("cells[%d].indices_inside[0] = %d\n",num, cell[num].indices_inside[0]);
	for(jj=0; jj<3; jj++)
		printf("cell[%d].Q[%d] = %lf\n",num, jj, cell[num].Q[jj]);
	for(jj=0; jj<9; jj++)
		printf("cell[%d].M3[%d] = %lf\n",num, jj, cell[num].M3[jj]);
	for(jj=0; jj<6; jj++)
		printf("cell[%d].PIJ[%d] = %lf\n",num,jj,cell[num].PIJ[jj]);
	//for(jj=0; jj<9; jj++)
//		printf("rhs0[%d] = %lf\n",jj, rhs0[jj]);
	//for(jj=0; jj<9; jj++)
	//	for(kk=0; kk<9; kk++)
	//		printf("A0[%d][%d] = %lf\n",jj, kk, A0[jj][kk]);
  //      exit(1);
     }
   }
	 cell[num].c[0] = x[0];
	 cell[num].c[1] = x[1];
	 cell[num].c[2] = x[2];
	 cell[num].c[3] = x[3];
	 cell[num].c[4] = x[4];
	 cell[num].c[5] = x[5];

	 cell[num].gamma[0] = x[6];
	 cell[num].gamma[1] = x[7];
	 cell[num].gamma[2] = x[8];

	}
      }
    }
   }
   free(x); free(rhs);
    for(int i=0;i<9;i++)
        delete [] A[i];
   delete  [] A;
}


void Position_FP_new2(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  int ii;//, j, k;
   double U_mean[3];
int num;
double energy;
//#pragma omp parallel for  private(energy, num, U_mean)
   for(ii=0; ii<(*gas).N; ii++)
   {
      num = index[ii];
      U_mean[0] =  cells[num].U_space[0];
      U_mean[1] =  cells[num].U_space[1];
      U_mean[2] =  cells[num].U_space[2];
      energy =    (U1[ii]-  cells[num].U_space[0])*(U1[ii]-  cells[num].U_space[0])
		+ (U2[ii]-  cells[num].U_space[1])*(U2[ii]-  cells[num].U_space[1])
		+ (U3[ii]-  cells[num].U_space[2])*(U3[ii]-  cells[num].U_space[2]);
//      x1[ii] =  x1[ii] + (*gas).delta_t*(U1[ii]+ (U1[ii]-  U_mean[0])*cells[num].alpha[0] + (U2[ii]- U_mean[1])*cells[num].alpha[1] + (U3[ii]- U_mean[2])*cells[num].alpha[2]
//			+ cells[num].beta[0]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) );
      x2[ii] =  x2[ii] + (*gas).delta_t*U2[ii]+(*gas).delta_t*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
			+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) );
//       x3[ii] =  x3[ii] + (*gas).delta_t*(U3[ii]+ (U1[ii]-  U_mean[0])*cells[num].alpha[6] + (U2[ii]- U_mean[1])*cells[num].alpha[7] + (U3[ii]- U_mean[2])*cells[num].alpha[8] );
   }
}

/*
void Position_FP_new3(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt){
  int ii,i;//, j, k;
  int dist, id_temp1, id_temp2, num_temp1, num_temp2, num_temp, id_temp;
  double  x1_old, x2_old;
   double U_mean[3];
int num, num_1;
double energy;
double U1_eff, U2_eff;
int id_0[3], id_1[3];
int j,k,l;
double sign[3];
double dummy;
if( (*gas).model == "FP_dense" )
      dummy=1.0;
else
      dummy=0.0;

//#pragma omp parallel for  private(energy, num, U_mean)
   for(ii=0; ii<(*gas).N; ii++)
   {
      num = index[ii];

      U_mean[0] =  cells[num].U_space[0];
      U_mean[1] =  cells[num].U_space[1];
      U_mean[2] =  cells[num].U_space[2];
      energy =    (U1[ii]-  cells[num].U_space[0])*(U1[ii]-  cells[num].U_space[0])
		+ (U2[ii]-  cells[num].U_space[1])*(U2[ii]-  cells[num].U_space[1])
		+ (U3[ii]-  cells[num].U_space[2])*(U3[ii]-  cells[num].U_space[2]);
      if( (*box).direction[0] == 1 ){
            x1[ii] =  x1[ii] + dt*U1[ii]
                  +dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[0] + (U2[ii]- U_mean[1])*cells[num].alpha[1] + (U3[ii]- U_mean[2])*cells[num].alpha[2]
                  + cells[num].beta[0]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) )
                  + dummy*cells[num].gamma_st*( (U1[ii]-  cells[num].U_space[1])*energy - cells[num].M3[0]- cells[num].M3[3]- cells[num].M3[5] )*dt;

            x1[ii] = x1[ii] - dummy*dt*cells[num].phi[0];
      }
      if( (*box).direction[1] == 1 ){

//	  x2[ii] =  x2[ii] + dt*U2[ii]
//                +dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
//			+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) );

// this is the working implementation.
	   if( (*box).problem == inverted ){
		//if( x2[ii] >  (*box).ghost && x2[ii] < (*box).Len[1] -  (*box).ghost ){
		//if(cells[num].cell_center[1]> (*box).ghost && cells[num].cell_center[1] < (*box).Len[1]-(*box).ghost){
           		//x2[ii] =  x2[ii] + dt*U2[ii]
                	//	+dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
			//		+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) );
			x2[ii] =  x2[ii] + dt*U2[ii];

			//if( x2[ii] >  2.0*(*gas).sigma && x2[ii] < (*box).Len[1] -  2.0*(*gas).sigma ){
				x2[ii] =  x2[ii] +dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
					+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) );
		    		+ dummy*cells[num].gamma_st*( (U2[ii]-  cells[num].U_space[1])*energy - cells[num].M3[1]- cells[num].M3[6]- cells[num].M3[8] )*dt;
				x2[ii] = x2[ii] - dummy*dt*cells[num].phi[1];
			//}

		//}
	   }
	   else{
		x2[ii] =  x2[ii] + dt*U2[ii]
                	+dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
			+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) )
		    	+ dummy*cells[num].gamma_st*( (U2[ii]-  cells[num].U_space[1])*energy - cells[num].M3[1]- cells[num].M3[6]- cells[num].M3[8] )*dt;
		x2[ii] = x2[ii] - dummy*dt*cells[num].phi[1];
	   }

      }
   }
}
*/

void DFP_extra_streaming(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt){
  int ii,i;//, j, k;
  int dist, id_temp1, id_temp2, num_temp1, num_temp2, num_temp, id_temp;
  double  x1_old, x2_old;
   double U_mean[3];
int num, num_1;
double energy;
double U1_eff, U2_eff;
int id_0[3], id_1[3];
int j,k,l;
double sign[3];
double dummy;
if( (*gas).model == "FP_dense" )
      dummy=1.0;
else
      dummy=0.0;

//#pragma omp parallel for  private(energy, num, U_mean)
   for(ii=0; ii<(*gas).N; ii++)
   {
      num = index[ii];

      U_mean[0] =  cells[num].U_space[0];
      U_mean[1] =  cells[num].U_space[1];
      U_mean[2] =  cells[num].U_space[2];
      energy =    (U1[ii]-  cells[num].U_space[0])*(U1[ii]-  cells[num].U_space[0])
		+ (U2[ii]-  cells[num].U_space[1])*(U2[ii]-  cells[num].U_space[1])
		+ (U3[ii]-  cells[num].U_space[2])*(U3[ii]-  cells[num].U_space[2]);
      if( (*box).direction[0] == 1 ){
            x1[ii] =  x1[ii]
                  +dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[0] + (U2[ii]- U_mean[1])*cells[num].alpha[1] + (U3[ii]- U_mean[2])*cells[num].alpha[2]
                  + cells[num].beta[0]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) )
                  + dummy*cells[num].gamma_st*( (U1[ii]-  cells[num].U_space[1])*energy - cells[num].M3[0]- cells[num].M3[3]- cells[num].M3[5] )*dt;
            x1[ii] = x1[ii] - dummy*dt*cells[num].phi[0];
      }
      if( (*box).direction[1] == 1 ){
				x2[ii] =  x2[ii] +dummy*dt*( (U1[ii]-  U_mean[0])*cells[num].alpha[3] + (U2[ii]- U_mean[1])*cells[num].alpha[4] + (U3[ii]- U_mean[2])*cells[num].alpha[5]
									+ cells[num].beta[1]*( energy - 3.0*(*gas).kb*cells[num].T/(*gas).m ) )
		    					+ dummy*cells[num].gamma_st*( (U2[ii]-  cells[num].U_space[1])*energy - cells[num].M3[1]- cells[num].M3[6]- cells[num].M3[8] )*dt;
				x2[ii] = x2[ii] - dummy*dt*cells[num].phi[1];
      }
   }
}

void Position_FP_ideal_gas(double *U1,double *U2,double *U3, double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt){
  int ii;
//#pragma omp parallel for  private(energy, num, U_mean)
  if( (*box).direction[0] == 1 ){
      for(ii=0; ii<(*gas).N; ii++)
           x1[ii] =  x1[ii] + dt*U1[ii];
   }
   if( (*box).direction[1] == 1 ){
      for(ii=0; ii<(*gas).N; ii++)
	  if((*box).problem == inverted){
	 	 //if( x2[ii] >  (*box).ghost && x2[ii] < (*box).Len[1] -  (*box).ghost ){
            		x2[ii] =  x2[ii] + dt*U2[ii];
		//}
	  }
	  else
		x2[ii] =  x2[ii] + dt*U2[ii];
    }
    if( (*box).direction[2] == 1 ){
      for(ii=0; ii<(*gas).N; ii++)
             x3[ii] =  x3[ii] + dt*U3[ii];
    }
}

void Velocity_FP_opt(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3){
//omp_set_dynamic(0);     // Explicitly disable dynamic teams
//omp_set_num_threads( (*box).num_thr);
//printf("(*gas).T0 = %e\n", (*gas).T0);
  std::random_device rd;
  std::mt19937 gen(rd());
  double std = 1.0;
  double mean = 0.0;
  double tau;
  double mu;
double Y;
int i, num, id;
int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
// double *xi_1, *xi_2, *xi_3;
// double *xi_telda_1, *xi_telda_2, *xi_telda_3;
int num_inside;
// double sum_1, sum_2, sum_3;
double frac1 = 0.0, frac2 = 0.0;
double *Mp1, *Mp2, *Mp3;
double A1, A2, A3;
double energy;
double dummy[3];
double sqrt_fr1fr2;
  std::normal_distribution<> dis1(mean,std);
  std::normal_distribution<> dis2(mean,std);
  std::normal_distribution<> dis3(mean,std);
//#pragma omp parallel for  private(num_inside, xi_1, xi_2, xi_3, xi_telda_1, xi_telda_2, xi_telda_3, Mp1, Mp2, Mp3, sum_1, sum_2, sum_3, i, id, frac1, frac2, A1, A2, A3, dummy, tau, mu, energy)
for(num=0; num<num_cell; num++){
  num_inside = cells[num].num_inside;
  Mp1 = (double *) malloc(  num_inside * sizeof(double) );
  Mp2 = (double *) malloc(  num_inside * sizeof(double) );
  Mp3 = (double *) malloc(  num_inside * sizeof(double) );

  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
	Mp1[i] = M1[id] - cells[num].U_space[0];
	Mp2[i] = M2[id] - cells[num].U_space[1];
	Mp3[i] = M3[id] - cells[num].U_space[2];
  }
  frac1 = 0.0; frac2 = 0.0;
  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
 	scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
	tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );
        if( (*gas).model == "FP_dense"){
                 Y = increase_collision_rate(cells[num].n, gas);
        }
        else{
		Y = 1.0;
        }
	//tau = tau/Y;
        //double mu_crr = ( 1.0 + 4.0*cells[num].n*(*gas).b*Y/5.0 + 0.761*pow(cells[num].n*(*gas).b*Y,2.0) )/Y;

        //if( (*gas).model == "FP_dense")
	//	tau = tau/increase_collision_rate(cells[num].n, gas);
	//	tau = tau*mu_crr;
	A1 = exp( -(*gas).delta_t/tau );
	A3 = sqrt(  (1.0-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  );
	A2 = (1.0-exp(-(*gas).delta_t/tau))*tau;
	energy = Mp1[i]*Mp1[i] + Mp2[i]*Mp2[i] + Mp3[i]*Mp3[i];
	frac1 += energy;

        dummy[0] = cells[num].c[0]*Mp1[i] + cells[num].c[1]*Mp2[i] +  cells[num].c[2]*Mp3[i];
	dummy[0] += cells[num].gamma[0]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

	dummy[1] = cells[num].c[1]*Mp1[i] + cells[num].c[3]*Mp2[i] +  cells[num].c[4]*Mp3[i];
	dummy[1] += cells[num].gamma[1]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

	dummy[2] = cells[num].c[2]*Mp1[i] + cells[num].c[4]*Mp2[i] +  cells[num].c[5]*Mp3[i];
	dummy[2] += cells[num].gamma[2]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

      	Mp1[i] = A1 * Mp1[i]  + A3*xi_1[id];
        Mp2[i] = A1 * Mp2[i]  + A3*xi_2[id];
        Mp3[i] = A1 * Mp3[i]  + A3*xi_3[id];

// 	Mp1[i] = A1 * Mp1[i] + A2*dummy[0] + A3*xi_1[id];
//         Mp2[i] = A1 * Mp2[i] + A2*dummy[1] + A3*xi_2[id];
//         Mp3[i] = A1 * Mp3[i] + A2*dummy[2] + A3*xi_3[id];
  }

  for(i=0; i<num_inside; i++){
    frac2 += Mp1[i]*Mp1[i] + Mp2[i]*Mp2[i] + Mp3[i]*Mp3[i];}
  sqrt_fr1fr2 = sqrt(fabs(frac1/frac2));
  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
	M1[id] = Mp1[i]*sqrt_fr1fr2 + cells[num].U_space[0];
	M2[id] = Mp2[i]*sqrt_fr1fr2 + cells[num].U_space[1];
	M3[id] = Mp3[i]*sqrt_fr1fr2 + cells[num].U_space[2];
  }
  free(Mp1); free(Mp2); free(Mp3);
}

}

void Velocity_FP_parallel(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *Mp1, double *Mp2, double *Mp3){
  double std = 1.0;
  double mean = 0.0;
  double tau;
  double mu;
  double Y;
int num, id;
int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
double *frac1 = (double *) malloc(  num_cell * sizeof(double) );
double *frac2 = (double *) malloc(  num_cell * sizeof(double) );
double A1, A2, A3;
double energy;
double dummy[3];
double sqrt_fr1fr2;


for(num=0; num<num_cell; num++){
  frac1[num] = 0.0;
  frac2[num] = 0.0;
}

for(id=0; id<(*gas).N; id++){
        num = index[id];
	Mp1[id] = M1[id]-cells[num].U_space[0];
        Mp2[id] = M2[id]-cells[num].U_space[1];
        Mp3[id] = M3[id]-cells[num].U_space[2];
	frac1[num] += Mp1[id]*Mp1[id] + Mp2[id]*Mp2[id] + Mp3[id]*Mp3[id];
}

//#pragma omp parallel for  private(num, mu, tau, A1, A2, A3, energy, dummy)
  for(id=0; id< (*gas).N; id++){
        num = index[id];

	if( cells[num].num_inside>1 ){
	 	scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
		tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );

        	if( (*gas).model == "FP_dense"){
        	         Y = increase_collision_rate(cells[num].n, gas);
        	}
        	else{
			Y = 1.0;
        	}
		tau = tau/Y;

		A1 = exp( -(*gas).delta_t/tau );
		A3 = sqrt(  (1.0-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  );
		A2 = (1.0-exp(-(*gas).delta_t/tau))*tau;

		energy = Mp1[id]*Mp1[id] + Mp2[id]*Mp2[id] + Mp3[id]*Mp3[id];

      		dummy[0] = cells[num].c[0]*Mp1[id] + cells[num].c[1]*Mp2[id] +  cells[num].c[2]*Mp3[id];
		dummy[0] += cells[num].gamma[0]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);
      		dummy[0] += cells[num].Lambda*(Mp1[id]*energy-cells[num].Q[0]);

		dummy[1] = cells[num].c[1]*Mp1[id] + cells[num].c[3]*Mp2[id] +  cells[num].c[4]*Mp3[id];
		dummy[1] += cells[num].gamma[1]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);
      		dummy[1] += cells[num].Lambda*(Mp2[id]*energy-cells[num].Q[1]);

		dummy[2] = cells[num].c[2]*Mp1[id] + cells[num].c[4]*Mp2[id] +  cells[num].c[5]*Mp3[id];
		dummy[2] += cells[num].gamma[2]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);
      		dummy[2] += cells[num].Lambda*(Mp3[id]*energy-cells[num].Q[2]);

	  	Mp1[id] = A1 * Mp1[id] + A2*dummy[0] + A3*xi_1[id];
        	Mp2[id] = A1 * Mp2[id] + A2*dummy[1] + A3*xi_2[id];
        	Mp3[id] = A1 * Mp3[id] + A2*dummy[2] + A3*xi_3[id];
	}
	else{
		Mp1[id] = 0.0;
        	Mp2[id] = 0.0;
        	Mp3[id] = 0.0;
	}
        if(Mp2[id] != Mp2[id]){
       		printf("Here in cell %d we have a nan!, velocity ideal\n",num);
		printf("A1 = %lf\n", A1);
		printf("A2 = %lf\n", A2);
		printf("A3 = %lf\n", A3);
		printf("dummy[0] = %lf\n",dummy[0]);
		printf("dummy[1] = %lf\n",dummy[1]);
		printf("dummy[2] = %lf\n",dummy[2]);
		printf("xi_1[%d]=%lf\n", id, xi_1[id]);
		printf("xi_2[%d]=%lf\n", id, xi_2[id]);
		printf("xi_3[%d]=%lf\n", id, xi_3[id]);
        	 exit(2);
       }
  }

for(id=0; id<(*gas).N; id++){
        num = index[id];
	frac2[num] += Mp1[id]*Mp1[id] + Mp2[id]*Mp2[id] + Mp3[id]*Mp3[id];
}
 // #pragma omp parallel for  private(num, sqrt_fr1fr2)
int doit = 0;
  for(id=0; id<(*gas).N; id++){
        num = index[id];
	doit = 1;
	if( (*box).problem == inverted ){
		//if(cells[num].cell_center[1]< (*box).ghost || cells[num].cell_center[1] > (*box).Len[1]-(*box).ghost){
		//if(cells[num].cell_center[1]< 2.0*(*gas).sigma || cells[num].cell_center[1] > (*box).Len[1]-2.0*(*gas).sigma){
		//	doit = 0;
		//}
	}
	if(doit == 1){
		if( cells[num].num_inside>1 ){
			sqrt_fr1fr2 = sqrt(fabs(frac1[num]/frac2[num]));
			M1[id] = Mp1[id]*sqrt_fr1fr2 + cells[num].U_space[0];
			M2[id] = Mp2[id]*sqrt_fr1fr2 + cells[num].U_space[1];
			M3[id] = Mp3[id]*sqrt_fr1fr2 + cells[num].U_space[2];
		}
		else{
			M1[id] =  cells[num].U_space[0];
			M2[id] =  cells[num].U_space[1];
			M3[id] =  cells[num].U_space[2];
		}
  	}
  }
free(frac1); free(frac2);
}



void Velocity_FP_linear(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *Mp1, double *Mp2, double *Mp3){
//omp_set_dynamic(0);     // Explicitly disable dynamic teams
//omp_set_num_threads( (*box).num_thr);
//printf("(*gas).T0 = %e\n", (*gas).T0);
  std::random_device rd;
  std::mt19937 gen(rd());
  double std = 1.0;
  double mean = 0.0;
  double tau;
  double mu;
  double Y;
int num, id;
int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
// double *xi_1, *xi_2, *xi_3;
// double *xi_telda_1, *xi_telda_2, *xi_telda_3;

// double sum_1, sum_2, sum_3;
double *frac1 = (double *) malloc(  num_cell * sizeof(double) );
double *frac2 = (double *) malloc(  num_cell * sizeof(double) );
double A1, A2, A3;
double F;
double energy;
double ru, rl;
double sqrt_fr1fr2;
  std::normal_distribution<> dis1(mean,std);
  std::normal_distribution<> dis2(mean,std);
  std::normal_distribution<> dis3(mean,std);


for(num=0; num<num_cell; num++){
  frac1[num] = 0.0;
  frac2[num] = 0.0;
}

for(id=0; id<(*gas).N; id++){
        num = index[id];
	Mp1[id] = M1[id]-cells[num].U_space[0];
        Mp2[id] = M2[id]-cells[num].U_space[1];
        Mp3[id] = M3[id]-cells[num].U_space[2];
	      frac1[num] += Mp1[id]*Mp1[id] + Mp2[id]*Mp2[id] + Mp3[id]*Mp3[id];
}

//#pragma omp parallel for  private(num, mu, tau, A1, A2, A3, energy, dummy)
  for(id=0; id< (*gas).N; id++){
        num = index[id];

 	  scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
	  tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );

        if( (*gas).model == "FP_dense"){
                 Y = increase_collision_rate(cells[num].n, gas);
        }
        else{
		Y = 1.0;
        }
	  //Y = 1.0;
	  tau = tau/Y;
        A2 = (1.0-exp(-(*gas).delta_t/tau))*tau;
	  A1 = exp( -(*gas).delta_t/tau );
	  A3 = sqrt(  (1.0-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  );
        F = 0.0;
        if(num == 0){
           F = (cells[num+1].n-cells[num].n+cells[num+2].n/128.0-cells[num].n/128.0)/(*box).delta_dim[1];
        }
        else if(num == 1){
             F = (cells[num+1].n-cells[num-1].n+cells[num+2].n/128.0-cells[num].n/128.0)/(*box).delta_dim[1];
        }
        else if(num == (*box).N[1]-2){
           F = (cells[num+1].n-cells[num-1].n+cells[num].n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
        else if(num == (*box).N[1]-1){
           F = (cells[num].n-cells[num-1].n+cells[num].n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
        else{
           F = (cells[num+1].n-cells[num-1].n+cells[num+2].n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
//4.47478200e-21*
         F = 4.0*3.14*F*pow((*gas).sigma,3.0)*4.47478200e-21/(*gas).m;

	  Mp1[id] = A1 * Mp1[id] + A3*xi_1[id] ;//+ F*(*gas).delta_t;
        Mp2[id] = A1 * Mp2[id] + A3*xi_2[id] ;//*(*gas).delta_t;//*tau*(1.0-A1);
        Mp3[id] = A1 * Mp3[id] + A3*xi_3[id] ;//+ F*(*gas).delta_t;//F*tau*(1.0-A1);

       if(Mp2[id] != Mp2[id]){
         printf("Here, we have a nan!, velocity ideal\n");
         exit(2);
       }
  }

for(id=0; id<(*gas).N; id++){
        num = index[id];
	frac2[num] += Mp1[id]*Mp1[id] + Mp2[id]*Mp2[id] + Mp3[id]*Mp3[id];
}
 // #pragma omp parallel for  private(num, sqrt_fr1fr2)
  for(id=0; id<(*gas).N; id++){
        num = index[id];
	sqrt_fr1fr2 = sqrt(fabs(frac1[num]/frac2[num]));
//      F = (*gas).Fn*(cells[num].M_n[3] - cells[num].M_n[2])/( (*box).delta_dim[0]*(*box).delta_dim[2]*(*gas).delta_t );

	M1[id] = Mp1[id]*sqrt_fr1fr2 + cells[num].U_space[0];
	M2[id] = Mp2[id]*sqrt_fr1fr2 + cells[num].U_space[1];// - F*(*gas).delta_t;
	M3[id] = Mp3[id]*sqrt_fr1fr2 + cells[num].U_space[2];
  }
free(frac1); free(frac2);
}


void Velocity_FP(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index){
  // original with effect of dense flows included in tau but not any change in the prandtl number
  std::random_device rd;
  std::mt19937 gen(rd());
  double std = 1.0;
  double mean = 0.0;
	 double Y, rho;
  double tau;
  double mu;
//   double T_here;

//     tau = 2.0*mu/( (*gas).n*(*gas).kb*(*gas).T );
//   tau = (*gas).lambda/( sqrt((*gas).kb*(*gas).T/(*gas).m) );
double *xi_1, *xi_2, *xi_3;

double *f_xi_1, *f_xi_2, *f_xi_3;
double *s_xi_1, *s_xi_2, *s_xi_3;

double *xi_telda_1, *xi_telda_2, *xi_telda_3;
// double *xi_sq_1, *xi_sq_2, *xi_sq_3;
double *sum_1, *sum_2, *sum_3;
xi_1 = (double *) malloc(  (*gas).N * sizeof(double) );
xi_2 = (double *) malloc(  (*gas).N * sizeof(double) );
xi_3 = (double *) malloc(  (*gas).N * sizeof(double) );
f_xi_1 = (double *) malloc(  (*gas).N * sizeof(double) );
f_xi_2 = (double *) malloc(  (*gas).N * sizeof(double) );
f_xi_3 = (double *) malloc(  (*gas).N * sizeof(double) );
s_xi_1 = (double *) malloc(  (*gas).N * sizeof(double) );
s_xi_2 = (double *) malloc(  (*gas).N * sizeof(double) );
s_xi_3 = (double *) malloc(  (*gas).N * sizeof(double) );
  std::normal_distribution<> dis1(mean,std);
  std::normal_distribution<> dis2(mean,std);
  std::normal_distribution<> dis3(mean,std);
    int i, num;
    for(i=0; i<(*gas).N; i++){
      xi_1[i] = 1.0*dis1(gen);
      xi_2[i] = 1.0*dis2(gen);
      xi_3[i] = 1.0*dis3(gen);
    }
 double * f_M1, *f_M2, *f_M3;
double * s_M1, *s_M2, *s_M3;

f_M1 = (double *) malloc(  (*gas).N * sizeof(double) );
f_M2 = (double *) malloc(  (*gas).N * sizeof(double) );
f_M3 = (double *) malloc(  (*gas).N * sizeof(double) );

s_M1 = (double *) malloc(  (*gas).N * sizeof(double) );
s_M2 = (double *) malloc(  (*gas).N * sizeof(double) );
s_M3 = (double *) malloc(  (*gas).N * sizeof(double) );



int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
sum_1 = (double *) malloc(  num_cell * sizeof(double) );
sum_2 = (double *) malloc(  num_cell * sizeof(double) );
sum_3 = (double *) malloc(  num_cell * sizeof(double) );
xi_telda_1 = (double *) malloc(  (*gas).N * sizeof(double) );
xi_telda_2 = (double *) malloc(  (*gas).N * sizeof(double) );
xi_telda_3 = (double *) malloc(  (*gas).N * sizeof(double) );

//#pragma omp parallel for
for(num=0; num<num_cell; num++){
  sum_1[num] = 0.0;
  sum_2[num] = 0.0;
  sum_3[num] = 0.0;
}

for(i=0; i<(*gas).N; i++){
  num = index[i];
   sum_1[num] +=  xi_1[i];
   sum_2[num] +=  xi_2[i];
   sum_3[num] +=  xi_3[i];
}

//#pragma omp  parallel for  private(num)
for(i=0; i<(*gas).N; i++){
  num = index[i];
  xi_telda_1[i] = xi_1[i] - sum_1[num]/(1.0*cells[num].num_inside);
  xi_telda_2[i] = xi_2[i] - sum_2[num]/(1.0*cells[num].num_inside);
  xi_telda_3[i] = xi_3[i] - sum_3[num]/(1.0*cells[num].num_inside);
}
//#pragma omp  parallel for
for(num=0; num<num_cell; num++){
  sum_1[num] = 0.0;
  sum_2[num] = 0.0;
  sum_3[num] = 0.0;
}
for(i=0; i<(*gas).N; i++){
  num = index[i];
  sum_1[num] += xi_telda_1[i]*xi_telda_1[i];
  sum_2[num] += xi_telda_2[i]*xi_telda_2[i];
  sum_3[num] += xi_telda_3[i]*xi_telda_3[i];
}
//#pragma omp  parallel for
for(num=0; num<num_cell; num++){
  sum_1[num] = sqrt( sum_1[num]/(1.0*cells[num].num_inside) );
  sum_2[num] = sqrt( sum_2[num]/(1.0*cells[num].num_inside) );
  sum_3[num] = sqrt( sum_3[num]/(1.0*cells[num].num_inside) );
}
//#pragma omp parallel for  private(num)
  for(i=0; i<(*gas).N; i++){
    num = index[i];
    xi_1[i] = xi_telda_1[i]/sum_1[num];
    xi_2[i] = xi_telda_2[i]/sum_2[num];
    xi_3[i] = xi_telda_3[i]/sum_3[num];
  }

for(i=0; i<(*gas).N; i++){
	num = index[i];
	mu =  (*gas).mu * pow( cells[num].T/(*gas).T0, 0.5 );
	tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );
        M1[i] =   M1[i] - (1.0 - exp( -(*gas).delta_t/tau)) *( M1[i] - cells[num].U_space[0]) + sqrt(  (1-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  )*xi_1[i] ;
        M2[i] =   M2[i] - (1.0 - exp( -(*gas).delta_t/tau)) *( M2[i] - cells[num].U_space[1]) + sqrt(  (1-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  )*xi_2[i] ;
        M3[i] =   M3[i] - (1.0 - exp( -(*gas).delta_t/tau)) *( M3[i] - cells[num].U_space[2]) + sqrt(  (1-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  )*xi_3[i] ;
    }
    free(sum_1);free(sum_2);free(sum_3);
    free(xi_1);free(xi_2);free(xi_3);
    free(xi_telda_1); free(xi_telda_2); free(xi_telda_3);
}



/*
void Velocity_FP_opt_mod(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *xi2_1, double *xi2_2, double *xi2_3){
//omp_set_dynamic(0);     // Explicitly disable dynamic teams
//omp_set_num_threads( (*box).num_thr);
//printf("(*gas).T0 = %e\n", (*gas).T0);
  std::random_device rd;
  std::mt19937 gen(rd());
  double std = 1.0;
  double mean = 0.0;
  double tau;
  double mu;

int i, num, id;
int num_cell = (*box).N[2]*(*box).N[1]*(*box).N[0];
// double *xi_1, *xi_2, *xi_3;
// double *xi_telda_1, *xi_telda_2, *xi_telda_3;
int num_inside;
// double sum_1, sum_2, sum_3;
double frac1 = 0.0, frac2 = 0.0;
double *Mp1, *Mp2, *Mp3;
double A1, A2, A3;
double energy;
double dummy[3];
double sqrt_fr1fr2;
long double A,B,C, es;
long double first, second;
  std::normal_distribution<> dis1(mean,std);
  std::normal_distribution<> dis2(mean,std);
  std::normal_distribution<> dis3(mean,std);
//#pragma omp parallel for  private(num_inside, xi_1, xi_2, xi_3, xi_telda_1, xi_telda_2, xi_telda_3, Mp1, Mp2, Mp3, sum_1, sum_2, sum_3, i, id, frac1, frac2, A1, A2, A3, dummy, tau, mu, energy)
for(num=0; num<num_cell; num++){
  num_inside = cells[num].num_inside;
  Mp1 = (double *) malloc(  num_inside * sizeof(double) );
  Mp2 = (double *) malloc(  num_inside * sizeof(double) );
  Mp3 = (double *) malloc(  num_inside * sizeof(double) );

  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
	Mp1[i] = M1[id] - cells[num].U_space[0];
	Mp2[i] = M2[id] - cells[num].U_space[1];
	Mp3[i] = M3[id] - cells[num].U_space[2];
  }
  frac1 = 0.0; frac2 = 0.0;
  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
 	scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
	tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );

	A1 = exp( -(*gas).delta_t/tau );
	A3 = sqrt(  (1.0-exp(-2.0*(*gas).delta_t/tau))*(*gas).kb*cells[num].T/(1.0*(*gas).m)  );
	A2 = (1.0-exp(-(*gas).delta_t/tau))*tau;

	es = (3.0/2.0)*(*gas).kb*cells[num].T/(*gas).m;
	A = (2.0*es/3.0)*( 1.0-exp(-2.0*(*gas).delta_t/tau) );
        B = (2.0*es/3.0)*( 2.0*(*gas).delta_t/tau - (1.0-exp(-(*gas).delta_t/tau))*(3.0-exp(-(*gas).delta_t/tau)) );
	C = (2.0*es/3.0)*pow( 1.0-exp(-(*gas).delta_t/tau), 2.0 );

	energy = Mp1[i]*Mp1[i] + Mp2[i]*Mp2[i] + Mp3[i]*Mp3[i];
	frac1 += energy;

        dummy[0] = cells[num].c[0]*Mp1[i] + cells[num].c[1]*Mp2[i] +  cells[num].c[2]*Mp3[i];
	dummy[0] += cells[num].gamma[0]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

	dummy[1] = cells[num].c[1]*Mp1[i] + cells[num].c[3]*Mp2[i] +  cells[num].c[4]*Mp3[i];
	dummy[1] += cells[num].gamma[1]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

	dummy[2] = cells[num].c[2]*Mp1[i] + cells[num].c[4]*Mp2[i] +  cells[num].c[5]*Mp3[i];
	dummy[2] += cells[num].gamma[2]*(energy - 3.0*(*gas).kb*cells[num].T/(*gas).m);

	first = sqrt(C*C/B);
	second = sqrt(A - C*C/B);

	Mp1[i] = A1 * Mp1[i]  + first*xi_1[id] + second*xi2_1[id];
        Mp2[i] = A1 * Mp2[i]  + first*xi_2[id] + second*xi2_2[id];
        Mp3[i] = A1 * Mp3[i]  + first*xi_3[id] + second*xi2_3[id];

  }

  for(i=0; i<num_inside; i++){
    frac2 += Mp1[i]*Mp1[i] + Mp2[i]*Mp2[i] + Mp3[i]*Mp3[i];}
  sqrt_fr1fr2 = sqrt(fabs(frac1/frac2));
//   sqrt_fr1fr2 = 1.0;
  for(i=0; i<num_inside; i++){
        id = cells[num].indices_inside[i];
	M1[id] = Mp1[i]*sqrt_fr1fr2 + cells[num].U_space[0];
	M2[id] = Mp2[i]*sqrt_fr1fr2 + cells[num].U_space[1];
	M3[id] = Mp3[i]*sqrt_fr1fr2 + cells[num].U_space[2];
  }
  free(Mp1); free(Mp2); free(Mp3);
}

}*/

// void Velocity_FP_opt_mod(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *xi2_1, double *xi2_2, double *xi2_3){
//
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   double std = 1.0;
//   double mean = 0.0;
//   double tau;
//   double mu;
//
// int  num, id;
// double A1;
// long double A,B,C, es;
// long double first, second;
//   std::normal_distribution<> dis1(mean,std);
//   std::normal_distribution<> dis2(mean,std);
//   std::normal_distribution<> dis3(mean,std);
// //#pragma omp parallel for  private(num_inside, xi_1, xi_2, xi_3, xi_telda_1, xi_telda_2, xi_telda_3, Mp1, Mp2, Mp3, sum_1, sum_2, sum_3, i, id, frac1, frac2, A1, A2, A3, dummy, tau, mu, energy)
//
//   for(id=0; id<(*gas).N; id++){
//         num = index[id];
//  	scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
// 	tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );
//
// 	es = (3.0/2.0)*(*gas).kb*cells[num].T/(*gas).m;
// 	A = (2.0*es/3.0)*( 1.0-exp(-2.0*(*gas).delta_t/tau) );
//         B = (2.0*es/3.0)*( 2.0*(*gas).delta_t/tau - (1.0-exp(-(*gas).delta_t/tau))*(3.0-exp(-(*gas).delta_t/tau)) );
// 	C = (2.0*es/3.0)*pow( 1.0-exp(-(*gas).delta_t/tau), 2.0 );
//
// 	A1 = 1.0 - exp( -(*gas).delta_t/tau );
//
// 	first = sqrt(C*C/B);
// 	second = sqrt(A - C*C/B);
// 	M1[id] = M1[id]  - A1*(M1[id] - cells[num].U_space[0]) + first*xi_1[id] + second*xi2_1[id];
//         M2[id] = M2[id]  - A1*(M2[id] - cells[num].U_space[1]) + first*xi_2[id] + second*xi2_2[id];
//         M3[id] = M3[id]  - A1*(M3[id] - cells[num].U_space[2]) + first*xi_3[id] + second*xi2_3[id];
//   }
// }
