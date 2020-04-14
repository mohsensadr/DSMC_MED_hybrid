#include "cpp_headers.h"
/*
#include <math.h>
#include <random>
#include <iostream>
#include <cmath>
#include <string>
#include <string>
#include <iostream>
#include <sstream>

#include <vector>
*/
void measure_pressure(double *U2,double *x2, double *x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells){
	int numO, numN, num1, num2, i, num;
	double xnew;
	//printf("Here, measuring pressure\n" );
	for(i=0; i<(*gas).N; i++){
		numO = floor( x2_old[i]/(*box).delta_dim[1] );
		numN = floor( x2[i]/(*box).delta_dim[1] );
		if(numN>numO){
			num1 = numO; num2 = numN;
		}
		else if(numN<numO){
			num1 = numN; num2 = numO;
		}
		else{
			num1 = -1; num2 = -1;
		}
		for(num=num1; num<num2; num++){
			//printf("Then, measuring pressure\n");
			cells[num].p_m += (*gas).Fn*(*gas).m*fabs(U2[i]);
		}
	}
}

void compute_evap_coeff_evaporation(double *U2,double *x2,double *x2_old,struct GAS *gas, struct BOX *box, int *color,int *n_ratio,  struct CELLS *cells){
	int i, num;
	double x2_prob[4];
	double delta;
	double flux;

	if((*box).problem==evaporation){

		for(i=0; i<4; i++){
			x2_prob[i] =  (*box).x2_prob[i];
		}
		//x2_prob[0] = (*box).0.9*(*gas).Lc;
		//x2_prob[1] = 1.1*(*gas).Lc;
		//x2_prob[2] = (*box).Len[1]-1.1*(*gas).Lh;
		//x2_prob[3] = (*box).Len[1]-0.9*(*gas).Lh;

		// initialize colors
		if((*box).step == (*box).after){
			(*gas).count_SP_c=0.0;
			(*gas).count_reem_c=0.0;
			(*gas).count_evap_c=0.0;
			for(i=0; i<(*gas).N; i++){
				if(x2[i]<x2_prob[0] && x2[i]>x2_prob[3]){
					color[i] = 3;
				}
				else{
					color[i] = 0;
				}
			}
		}


		for(i=0; i<(*gas).N; i++){
			// coloring based on crossing probes
			if(x2[i]<x2_prob[1] && x2_old[i]>x2_prob[1] && x2[i]>x2_prob[0]){
				color[i] = 1;// just succeeded in getting to phase transition thickness
			}
			if(  x2[i]<x2_prob[0] && x2[i]>x2_prob[3] ){
				color[i] = 3;// from vapour made its way to liquid
			}
			if(x2[i]>x2_prob[2] && x2_old[i]<x2_prob[2] &&  x2[i]<x2_prob[3]){
				color[i] = 1;// just succeeded in getting to phase transition thickness
			}

			// counting what we have
			if( ( x2[i] > x2_prob[1] && x2_old[i] < x2_prob[1] ) || ( x2[i]<x2_prob[2] && x2_old[i]>x2_prob[2] )){
				if(color[i]==1){
					(*gas).count_SP_c = (*gas).count_SP_c+1.0;
					color[i] = 0;
				}
				else if(color[i]==3){
					(*gas).count_evap_c = (*gas).count_evap_c+1.0;
					color[i] = 0;
				}
			}
		}

		// print what we have done
		if((*box).step%(*box).every == 0){
			num = floor((*box).x2_prob[1]/(*box).delta_dim[1]);
			double total = (*gas).count_SP_c + (*gas).count_reem_c + (*gas).count_evap_c;
			if(total>1e-15 ){
				printf("count_SP_c = %e, count_evap_c = %e\n",  (*gas).count_SP_c, (*gas).count_evap_c);
			}
		}

		/*
		delta = (*gas).Lc*0.01;

		x2_prob[0] = 30*1e-10;
		x2_prob[1] = 50*1e-10;
		x2_prob[2] = 70*1e-10;
		x2_prob[3] = 90*1e-10;
		for(i=0; i<(*gas).N; i++){
			//flux = fabs( (*gas).Fn*(*gas).m*n_ratio[i]*U2[i] );
			flux = 1.0;
			// Check the 1st prob
			if(x2[i]>x2_prob[0] && x2_old[i]<x2_prob[0]){
				color[i] = 1000;
				if((*box).step > (*gas).times*(*box).after){
					(*box).JR[0] += flux;
					//(*box).JL[0] -= flux;
					(*box).NR[0] += 1.0;
				}
			}
			if(x2[i]<x2_prob[0] && x2_old[i]>x2_prob[0]){
				(*box).JL[0] += flux;
				//color[i] = 0;
				//if(color[i] == 1000){
				//	(*box).JR[0] -= flux;
				//	color[i] = 0;
				//}
				//(*box).JR[0] -= flux;
				(*box).NL[0] += 1.0;
			}
			// Check the 2nd prob
			if(x2[i]>x2_prob[1] && x2_old[i]<x2_prob[1]){
				if((*box).step > (*gas).times*(*box).after){
					if(color[i] == 1000){
						(*box).JR[1] += flux;
						//(*box).JL[1] -= flux;
						(*box).NR[1] += 1.0;
					}
				}
				color[i] = 0;
			}
			if(x2[i]<x2_prob[1] && x2_old[i]>x2_prob[1]){
				(*box).JL[1] += flux;
				//(*box).JR[1] -= flux;
				(*box).NL[1] += 1.0;
				//color[i] = 0;
			}
			// Check the 3rd prob
			if(x2[i]>x2_prob[2] && x2_old[i]<x2_prob[2]){
				(*box).JR[2] += flux;
				//(*box).JL[2] -= flux;
				(*box).NR[2] += 1.0;
				//color[i] = 0;
			}
			if(x2[i]<x2_prob[2] && x2_old[i]>x2_prob[2]){
				(*box).JL[2] += flux;
			//(*box).JR[2] -= flux;
			(*box).NL[2] += 1.0;
			//color[i] = 0;
			}
			// Check the 4th prob
			if(x2[i]>x2_prob[3] && x2_old[i]<x2_prob[3]){
				(*box).JR[3] += flux;
				//(*box).JL[3] -= flux;
				(*box).NR[3] += 1.0;
				//color[i] = 0;
			}
			if(x2[i]<x2_prob[3] && x2_old[i]>x2_prob[3]){
				(*box).JL[3] += flux;
				//(*box).JR[3] -= flux;
				(*box).NL[3] += 1.0;
				//color[i] = 0;
			}

		}
		*/
	}
	else if( (*box).problem == inverted){
		for(i=0; i<4; i++){
			x2_prob[i] =  (*box).x2_prob[i];
		}
		//x2_prob[0] = (*box).0.9*(*gas).Lc;
		//x2_prob[1] = 1.1*(*gas).Lc;
		//x2_prob[2] = (*box).Len[1]-1.1*(*gas).Lh;
		//x2_prob[3] = (*box).Len[1]-0.9*(*gas).Lh;

		// initialize colors
		if((*box).step == (*box).after){
			(*gas).count_SP_c=0.0;
			(*gas).count_SP_h=0.0;
			(*gas).count_reem_c=0.0;
			(*gas).count_reem_h=0.0;
			(*gas).count_evap_c=0.0;
			(*gas).count_evap_h=0.0;
			for(i=0; i<(*gas).N; i++){
				if(x2[i]<x2_prob[0]){
					color[i] = 3;
				}
				else if(x2[i]>x2_prob[3]){
					color[i] = 30;
				}
				else{
					color[i] = 0;
				}
			}
		}


		for(i=0; i<(*gas).N; i++){
			/*
			// coloring based on crossing probes
			if(x2[i]<x2_prob[1] && x2_old[i]>x2_prob[1]){
				color[i] = 1;// just succeeded in getting to phase transition thickness
			}
			if( (color[i] == 1 && x2[i]<x2_prob[0] && x2_old[i]>x2_prob[0] ) || (x2[i]<x2_prob[0] && x2_old[i]>x2_prob[1]) ){
				color[i] = 2;// from vapour made its way to liquid
			}
			if(x2[i]>x2_prob[2] && x2_old[i]<x2_prob[2]){
				color[i] = 10;// just succeeded in getting to phase transition thickness
			}
			if( (color[i] == 10 && x2[i]>x2_prob[3] && x2_old[i]<x2_prob[3] ) || (x2[i]>x2_prob[3] && x2_old[i]<x2_prob[2]) ){
				color[i] = 20;// from vapour made its way to liquid
			}
			*/



			/*
			// counting what we have
			if(x2[i] > x2_prob[1] && x2_old[i] < x2_prob[1]){
				if(color[i]==1){
					(*gas).count_SP_c = (*gas).count_SP_c+1.0;
					color[i] = 0;
				}
				else if(color[i]==2){
					(*gas).count_reem_c = (*gas).count_reem_c+1.0;
					color[i] = 0;
				}
				else if(color[i]==3){
					(*gas).count_evap_c = (*gas).count_evap_c+1.0;
					color[i] = 0;
				}
			}
			if(x2[i]<x2_prob[2] && x2_old[i]>x2_prob[2]){
				if(color[i]==10){
					(*gas).count_SP_h = (*gas).count_SP_h+1.0;
					color[i] = 0;
				}
				else if(color[i]==20){
					(*gas).count_reem_h = (*gas).count_reem_h+1.0;
					color[i] = 0;
				}
				else if(color[i]==30){
					(*gas).count_evap_h = (*gas).count_evap_h+1.0;
					color[i] = 0;
				}
			}
			*/
			if(x2[i] > x2_prob[1] && color[i] == 3){
				(*gas).count_evap_c = (*gas).count_evap_c+1.0;
				color[i] = 0;
			}
			else if(x2[i]<x2_prob[2] && color[i] == 30){
				(*gas).count_evap_h = (*gas).count_evap_h+1.0;
				color[i] = 0;
			}
			if(x2[i]<x2_prob[0]){
				color[i] = 3;
			}
			else if(x2[i]>x2_prob[3]){
				color[i] = 30;
			}
		}

		// print what we have done
		if((*box).step%(*box).every == 0){
			num = floor((*box).x2_prob[1]/(*box).delta_dim[1]);
			double total = (*gas).count_SP_c + (*gas).count_reem_c + (*gas).count_evap_c;
			if(total>1e-15 ){
				printf("mdot_c = %e\n",(*gas).count_evap_c/(*gas).t_post);
				//printf("sigma_c = %lf,  alpha_c = %lf\n",  ((*gas).count_evap_c+(*gas).count_reem_c)/total, (*gas).count_reem_c/(total-(*gas).count_evap_c));
				printf(" Tvc = %lf, nvc*sigma^3 = %e\n", cells[num].T, cells[num].n*pow( (*gas).sigma, 3.0 ));
			}
			total = (*gas).count_SP_h + (*gas).count_reem_h + (*gas).count_evap_h;
			num = floor((*box).x2_prob[2]/(*box).delta_dim[1]);
			if(total>1e-15 ){
				printf("mdot_h = %e\n",(*gas).count_evap_h/(*gas).t_post);
				//printf("sigma_h = %lf,  alpha_h = %lf\n",  ((*gas).count_evap_h+(*gas).count_reem_h)/total, (*gas).count_reem_h/(total-(*gas).count_evap_h));
				printf(" Tvh = %lf, nvh*sigma^3 = %e\n", cells[num].T, cells[num].n*pow( (*gas).sigma, 3.0 ));
			}
		}
	}
}


void compute_evap_coeff(double *U2,double *x2,double *x2_old,struct GAS *gas, struct BOX *box, int *color,int *n_ratio){
	int i, num;
	double x2_prob[4];
	double delta;
	double flux;
	delta = (*gas).Lc*0.01;
	x2_prob[0] = 8.5*(*gas).sigma;//(*gas).Lc_center + (*gas).Lc/2.0 - delta;
	x2_prob[1] = 14.0*(*gas).sigma;//(*gas).Lc_center + (*gas).Lc/2.0 + delta;
	x2_prob[2] = 31.0*(*gas).sigma;//(*gas).Lh_center - (*gas).Lh/2.0 - delta;
	x2_prob[3] = 36.5*(*gas).sigma;//(*gas).Lh_center - (*gas).Lh/2.0 + delta;
	for(i=0; i<(*gas).N; i++){
		flux = fabs( n_ratio[i]*U2[i] );
		// Check the 1st probe which is in liquid
		if(x2[i]>x2_prob[0] && x2_old[i]<x2_prob[0]){
			color[i] = 1000;
			(*box).JR[0] += flux;
		}
		if(x2[i]<x2_prob[0] && x2_old[i]>x2_prob[0]){
			if( color[i]==1000 ){
				(*box).Jref[0] += flux;
			}
			color[i] = 0;
			if( color[i] == 100)
				(*box).JL[0] += flux;
		}
		// Check the 2nd probe which is in vapour
		if(x2[i]>x2_prob[1] && x2_old[i]<x2_prob[1]){
			if(color[i]==100){
				(*box).Jref[1] += flux;
			}
			else if(color[i]==1000){
				(*box).JR[1] += flux;
			}
			color[i] = 0;
		}
		if(x2[i]<x2_prob[1] && x2_old[i]>x2_prob[1]){
			color[i] = 100;
			(*box).JL[1] += flux;
		}
		// Check the 3rd prob which is in vapour
		if(x2[i]>x2_prob[2] && x2_old[i]<x2_prob[2]){
			color[i] = 200;
			(*box).JR[2] += flux;
		}
		if(x2[i]<x2_prob[2] && x2_old[i]>x2_prob[2]){
			if(color[i]==200){
				(*box).Jref[2] += flux;
			}
			else if(color[i]==2000){
				(*box).JL[2] += flux;
			}
			color[i] = 0;
		}
		// Check the 4th prob which is in liquid
		if(x2[i]>x2_prob[3] && x2_old[i]<x2_prob[3]){
			if( color[i]==2000 ){
				(*box).Jref[3] += flux;
			}
			color[i] = 0;
			(*box).JR[3] += flux;
		}
		if(x2[i]<x2_prob[3] && x2_old[i]>x2_prob[3]){
			color[i] = 2000;
			(*box).JL[3] += flux;
		}
	}
}


void reset_inverted_ghost(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j;
int  Ntot;
int num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double stdh = sqrt((*gas).kb*(*gas).Th/(*gas).m);
double std;
double dt_remain,yw;
int count_reem_c=0, count_reem_h=0;
double pi = acos(-1.0);
double tol = 1e-15, toln = 0.9;
double epsilon = 0.01*(*gas).sigma;

// find the minimum point of density
//double Lh = 0.0;
double Lc = 0.0, n_c, n_h, L0, n;
double disp;
double volume, u2;
int stt = floor((*box).ghost/(*box).delta_dim[1])+1;
n_c = cells[stt].n;//(*gas).nc;//
n_h = cells[(*box).N[1]-1-stt].n;
int donec=0, doneh=0;


int Noutghost = 0;

	for(i=0; i<(*gas).N; i++){
		(*flag)[i] = 1;
		if( (*x2)[i] < (*box).ghost || (*x2)[i] >(*box).Len[1]-(*box).ghost ){
			(*flag)[i] = -1;
			Noutghost ++;
		}
	}

	volume = (*box).Len[0]*(*box).Len[2]*((*box).ghost);
	int N_new_ghost_c = floor(volume*n_c/(*gas).Fn+ dis_x(gen));
	int N_new_ghost_h = floor(volume*n_h/(*gas).Fn+ dis_x(gen));

	 Ntot = (*gas).N - Noutghost + N_new_ghost_c + N_new_ghost_h;

	 double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew;
	 int *colornew;
	 x1new = (double *) malloc( Ntot * sizeof(double) );
	 x2new = (double *) malloc( Ntot * sizeof(double) );
	 U1new = (double *) malloc( Ntot * sizeof(double) );
	 U2new = (double *) malloc( Ntot * sizeof(double) );
	 U3new = (double *) malloc( Ntot * sizeof(double) );
	 x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
	 x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
	 colornew = (int *) malloc( Ntot * sizeof(int) );
	 j=0;
	 for(i=0; i<(*gas).N; i++){
	 	if((*flag)[i] != -1){
			x1new[j] = (*x1)[i];
	 		x2new[j] = (*x2)[i];
	 		U1new[j] = (*U1)[i];
	 		U2new[j] = (*U2)[i];
	 		U3new[j] = (*U3)[i];
			colornew[j] = (*color)[i];
			x1_oldnew[j] = (*x1_old)[i];
	 		x2_oldnew[j] = (*x2_old)[i];
	 		j++;
	 	}
	 }

	 for(i=0; i<N_new_ghost_c; i++){
	 		x2new[j] = dis_x(gen)*((*box).ghost);
			x2_oldnew[j] = tol;
			x1new[j] = (*x1)[0];
			x1_oldnew[j] = (*x1)[0];

			U1new[j] = dis_norm(gen)*stdc;// + cells[stt].U_space[0];
			U2new[j] = dis_norm(gen)*stdc;// + cells[stt].U_space[1];
			U3new[j] = dis_norm(gen)*stdc;// + cells[stt].U_space[2];
			colornew[j] = 3;
	 		j++;
	 }

	 for(i=0; i<N_new_ghost_h; i++){
	 		x2new[j] = (*box).Len[1] - dis_x(gen)*((*box).ghost);
			x2_oldnew[j] = (*box).Len[1] - tol;
			x1new[j] = (*x1)[0];
			x1_oldnew[j] = (*x1)[0];

			U1new[j] = dis_norm(gen)*stdh;// + cells[(*box).N[1]-1-stt].U_space[0];
			U2new[j] = dis_norm(gen)*stdh;// + cells[(*box).N[1]-1-stt].U_space[1];
			U3new[j] = dis_norm(gen)*stdh;// + cells[(*box).N[1]-1-stt].U_space[2];
			colornew[j] = 30;
	 		j++;
	 }

	  free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);free(*color); free(*x1); free(*x1_old);
	 (*gas).N = Ntot;
	 *x1 = x1new;
	 *x2 = x2new;
	 *U1 = U1new;
	 *U2 = U2new;
	 *U3 = U3new;
	 *x1_old = x1_oldnew;
	 *x2_old = x2_oldnew;
	 *color = colornew;

	 free(*n_ratio);
	 free(*flag); free(*index);
	 free(*xi_1); free(*xi_2); free(*xi_3); free(*xi_1f); free(*xi_2f); free(*xi_3f);
	 free(*Mp1); free(*Mp2); free(*Mp3);

	 *n_ratio = (int *) malloc( (*gas).N * sizeof(int) );
	 *index = (int *) malloc( (*gas).N * sizeof(int) );
	 *flag = (int *) malloc( (*gas).N * sizeof(int) );
	 *xi_1f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_2f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_3f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_1 = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_2 = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_3 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp1 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp2 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp3 = (double *) malloc( (*gas).N * sizeof(double) );



	 for(i=0; i<(*gas).N; i++){
	 	num = floor((*x2)[i]/(*box).delta_dim[1]);
	 	(*index)[i] = num;
	 	(*n_ratio)[i] = cells[num].n_ratio;
	 	if(num<0 || num>=(*box).N[1]){
	 		printf("-------> num = %d\n",num);
	 	}
	 }

}







void shock_BC(double *U2, double *x2, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells){
//std::random_device rd;
//std::mt19937 gen(rd());
//std::uniform_real_distribution<> dis_x(-1.0, 1.0);
//std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,Nout=0,num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double tol = 1e-15;


	for(i=0; i<(*gas).N; i++){
		if( x2[i] < tol || x2[i] > (*box).Len[1]-tol){
				if( x2[i] < tol ){
					x2[i] = x2[i] + 2.0*fabs( x2[i]  ) + tol;
					T[i] = (*gas).Tv1;
					U2[i] = fabs(U2[i]);//cells[0].U_space[1];
				}
				else{
					x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] ) - tol;
					T[i] = (*gas).Tv3;
					U2[i] = -fabs(U2[i]);//cells[(*box).N[1]-1].U_space[1];
				}

		}
	}
}


void evap_BC(double *U1,double *U2, double *U3, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(-1.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,Nout=0,num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double tol = 1e-15, toln = 0.9;


	for(i=0; i<(*gas).N; i++){
		if( x2[i] < tol || x2[i] > (*box).Len[1]-tol){
					//if( x2[i] < tol){
						Nout += 1.0;//U2[i];
					//}
					x2[i] = (*gas).Lc_center+dis_x(gen)*(*box).thermc/2.0;
					num = floor(x2[i]/(*box).delta_dim[1]);
					// next line is the correct one
					//U2[i] = dis_norm(gen)*stdc + cells[num].U_space[1];

					// this is for testing with Marcel only
					U1[i] = dis_norm(gen)*stdc;// + cells[num].U_space[1];
					U2[i] = dis_norm(gen)*stdc;// + cells[num].U_space[1];
					U3[i] = dis_norm(gen)*stdc;// + cells[num].U_space[1];

		}
	}
	if( (*box).step > (*box).after ){
		(*gas).N_abs += Nout;
		if((*box).step%(*box).every == 0){
			printf("**N_abs= %e, time=%e, N_abs/time = %e\n", (*gas).N_abs, (*gas).t_post,(*gas).N_abs/ (*gas).t_post);
		}
	}
}

int reset_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j,NoutC=0, NoutH=0;
double Cratio, Hratio;
int N_new, Ntot;
int num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double stdh = sqrt((*gas).kb*(*gas).Th/(*gas).m);
double std;
double dt_remain,yw;
int count_reem_c=0, count_reem_h=0;
double pi = acos(-1.0);
double tol = 1e-15, toln = 0.9;
double epsilon = 0.01*(*box).delta_dim[1];

// find the minimum point of density
double Lh = 0.0;
double Lc = 0.0, n_c, n_h, L0, n;
double disp;
double volume, u2;
int stt = floor((*box).ghost/(*box).delta_dim[1])+1;
n_c = cells[stt].n_smooth;//(*gas).nc;//
n_h = cells[(*box).N[1]-1-stt].n_smooth;
int donec=0, doneh=0, idnum;

int post = 0;
for(num=stt; num<(*box).N[1]-1-stt; num++){
	if(cells[num].n_smooth<toln*n_c && donec==0){
		Lc = cells[num].cell_center[1] + (*box).delta_dim[1]*0.5;
		donec = 1;
	}
	if(cells[num].n_smooth<toln*n_h){
		Lh = (*box).Len[1] - cells[num].cell_center[1] - (*box).delta_dim[1]*0.5;
	}
}
disp = (Lc - Lh)/2.0;//(*gas).Lc;

/*
double delnU = 0.0;
delnU = cells[0].n*cells[0].U_space[1] - cells[(*box).N[1]-1].n*cells[(*box).N[1]-1].U_space[1];
if(delnU>0.0){
	n = n_h;
}
else{
	n = n_c;
}
disp = delnU*(*gas).delta_t/n;
*/
//printf("--->     disp = %e*s, Lc=%e*s\n", disp/(*gas).sigma, Lc/(*gas).sigma);
int Noutghost = 0;
// check if we need to reset the frame work
//double ghost = (*box).ghost;
if(fabs(disp) > epsilon ){
	post = 1;
	printf("--->     disp = %e*s, Lc=%e*s, Lh=%e*s\n", disp/(*gas).sigma, Lc/(*gas).sigma, Lh/(*gas).sigma);
	//printf("--->     disp = %e, Lc=%e*s Lh=%e*s\n", disp, Lc/(*gas).sigma, Lh/(*gas).sigma);
	//printf("     ----  >>>  disp=  %e\n",disp);
	for(i=0; i<(*gas).N; i++){
		(*flag)[i] = 1;
		num = floor((*x2)[i]/(*box).delta_dim[1]);
		//if( (*x2)[i] < (*box).ghost || (*x2)[i] >(*box).Len[1]-(*box).ghost ){
		//	(*flag)[i] = -1;
		//	Noutghost ++;
		//}
		//if( (*x2)[i] > (*box).ghost && (*x2)[i] <(*box).Len[1]-(*box).ghost ){
				(*x2)[i] = (*x2)[i] - disp;
				(*x2_old)[i] = (*x2_old)[i] - disp;
				//if( (*x2)[i] < tol + (*box).ghost ){
				if( (*x2)[i] < tol ){
					NoutC++;
					(*flag)[i] = -1;
				}
				//else if( (*x2)[i] > (*box).Len[1]-tol -(*box).ghost ){
				else if( (*x2)[i] > (*box).Len[1]-tol){
					NoutH++;
					(*flag)[i] = -1;
				}
			//}

	}
	if(disp>0.0){
		idnum = (*box).N[1]-1-stt;
		L0 = (*box).Len[1];// - (*box).ghost;
		std = stdh;
		n = n_h;//(*gas).nh;
		//u2 = cells[0].U_space[1];//-n_c*cells[0].U_space[1]/n_h;
	}
	else{
		idnum = stt;
		L0 = 0.0;//(*box).ghost;
		std = stdc;
		n = n_c;//(*gas).nc;
		//u2 = cells[(*box).N[1]-1].U_space[1];//-n_h*cells[(*box).N[1]-1].U_space[1]/n_c;
	}
	 volume = (*box).Len[0]*(*box).Len[2]*fabs(disp);
	 //num = floor(L0-(disp+(*box).delta_dim[1]));
	 //if(num<0 || num > (*box).N[1]-1){
	 //	 printf("num = %d\n", num);
	 //}
	 //n = cells[num].n;

	 N_new = floor(volume*n/(*gas).Fn+ dis_x(gen));

	 //volume = (*box).Len[0]*(*box).Len[2]*((*box).ghost);
	 //int N_new_ghost_c = floor(volume*n_c/(*gas).Fn+ dis_x(gen));
	 //int N_new_ghost_h = floor(volume*n_h/(*gas).Fn+ dis_x(gen));
	 if(N_new<1)
	 	printf("N_new = %d\n", N_new);
	 //Ntot = N_new + (*gas).N - NoutC - NoutH - Noutghost + N_new_ghost_c + N_new_ghost_h;
	 Ntot = N_new + (*gas).N - NoutC - NoutH;

	 double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew;
	 int *colornew;
	 x1new = (double *) malloc( Ntot * sizeof(double) );
	 x2new = (double *) malloc( Ntot * sizeof(double) );
	 U1new = (double *) malloc( Ntot * sizeof(double) );
	 U2new = (double *) malloc( Ntot * sizeof(double) );
	 U3new = (double *) malloc( Ntot * sizeof(double) );
	 x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
	 x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
	 colornew = (int *) malloc( Ntot * sizeof(int) );
	 j=0;
	 for(i=0; i<(*gas).N; i++){
	 	if((*flag)[i] != -1){
			x1new[j] = (*x1)[i];
	 		x2new[j] = (*x2)[i];
	 		U1new[j] = (*U1)[i];
	 		U2new[j] = (*U2)[i];
	 		U3new[j] = (*U3)[i];
			colornew[j] = (*color)[i];
			x1_oldnew[j] = (*x1_old)[i];
	 		x2_oldnew[j] = (*x2_old)[i];
	 		j++;
	 	}
	 }

	 for(i=0; i<N_new; i++){
	 		x2new[j] = L0-dis_x(gen)*disp;
			x2_oldnew[j] = L0;
			x1new[j] = (*x1)[0];
			x1_oldnew[j] = (*x1)[0];

			//num = floor(x2new[j]/(*box).delta_dim[1]);

			U1new[j] = dis_norm(gen)*std + cells[idnum].U_space[0];
			U2new[j] = dis_norm(gen)*std + cells[idnum].U_space[1];
			U3new[j] = dis_norm(gen)*std + cells[idnum].U_space[2];

			if(disp>0.0)
				colornew[j] = 30;
			else
				colornew[j] = 3;
	 		j++;
	 }
	 /*
	 for(i=0; i<N_new_ghost_c; i++){
	 		x2new[j] = dis_x(gen)*((*box).ghost);
			x2_oldnew[j] = tol;
			x1new[j] = (*x1)[0];
			x1_oldnew[j] = (*x1)[0];

			num = floor(x2new[j]/(*box).delta_dim[1]);

			U1new[j] = dis_norm(gen)*stdc;// + cells[num].U_space[0];
			U3new[j] = dis_norm(gen)*stdc;// + cells[num].U_space[1];
			U2new[j] = dis_norm(gen)*stdc;// + cells[num].U_space[2];
			colornew[j] = 3;
	 		j++;
	 }

	 for(i=0; i<N_new_ghost_h; i++){
	 		x2new[j] = (*box).Len[1] - dis_x(gen)*((*box).ghost);
			x2_oldnew[j] = (*box).Len[1] - tol;
			x1new[j] = (*x1)[0];
			x1_oldnew[j] = (*x1)[0];

			num = floor(x2new[j]/(*box).delta_dim[1]);

			U1new[j] = dis_norm(gen)*stdh;// + cells[num].U_space[0];
			U3new[j] = dis_norm(gen)*stdh;// + cells[num].U_space[1];
			U2new[j] = dis_norm(gen)*stdh;// + cells[num].U_space[2];
			colornew[j] = 30;
	 		j++;
	 }
	 */
	  free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);free(*color); free(*x1); free(*x1_old);
	 (*gas).N = Ntot;
	 *x1 = x1new;
	 *x2 = x2new;
	 *U1 = U1new;
	 *U2 = U2new;
	 *U3 = U3new;
	 *x1_old = x1_oldnew;
	 *x2_old = x2_oldnew;
	 *color = colornew;

	 // Resize the arrays than will be initilized in other functions
	 // first free them
	 // then allocate memory
	 //free(*x1);free(*x1_old);
	 free(*n_ratio);
	 free(*flag); free(*index);
	 free(*xi_1); free(*xi_2); free(*xi_3); free(*xi_1f); free(*xi_2f); free(*xi_3f);
	 free(*Mp1); free(*Mp2); free(*Mp3);

	 *n_ratio = (int *) malloc( (*gas).N * sizeof(int) );
	 *index = (int *) malloc( (*gas).N * sizeof(int) );
	 *flag = (int *) malloc( (*gas).N * sizeof(int) );
	 *xi_1f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_2f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_3f = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_1 = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_2 = (double *) malloc( (*gas).N * sizeof(double) );
	 *xi_3 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp1 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp2 = (double *) malloc( (*gas).N * sizeof(double) );
	 *Mp3 = (double *) malloc( (*gas).N * sizeof(double) );


	 //*x1 = (double *) malloc( Ntot * sizeof(double) );
	 //*x1_old = (double *) malloc( Ntot * sizeof(double) );

	 // set n_ratio of particles to the cell's n_ratios
	 for(i=0; i<(*gas).N; i++){
	 	num = floor((*x2)[i]/(*box).delta_dim[1]);
	 	(*index)[i] = num;
	 	(*n_ratio)[i] = cells[num].n_ratio;
	 	if(num<0 || num>=(*box).N[1]){
	 		printf("-------> num = %d\n",num);
	 	}
	 }
}
	return post;
}




void ghost_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j,NoutC=0, NoutH=0;
double epsilon = 1e-15;
double Cratio, Hratio;
int NnewC,NnewH,N_new, Ntot;
int num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double stdh = sqrt((*gas).kb*(*gas).Th/(*gas).m);
double dt_remain,yw;
int count_reem_c=0, count_reem_h=0;
double pi = acos(-1.0);
double tol = (*box).ghost;
// particles that leave the domain needed to be removed
for(i=0; i<(*gas).N; i++){
	if( (*x2)[i] < tol ){
			NoutC++;
			(*flag)[i] = -1;
	}
	else if( (*x2)[i] > (*box).Len[1]-tol ){
			NoutH++;
			(*flag)[i] = -1;
	}
	else{
		(*flag)[i] = 1;
	}
}
// I should remove those with negative flag.

// Now, how many new particles enter domain
double beta, dummy;
//nvc = 0.00287287/pow(3.405e-10,3.0);
//nvh = 0.00412534/pow(3.405e-10,3.0);

double volume = (*box).Len[0]*(*box).Len[2]*tol;

//beta = 1.0/(stdc*sqrt(2.0));
//beta = 1.0/(stdc);
//dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nvc/(beta*2.0*sqrt(pi));
//Cratio = 0.9;
//Hratio = 1.1;
//NnewC =  floor(Cratio*NoutC);
//NnewH =  floor(Hratio*NoutH);
NnewC = floor(volume*(*gas).nc/(*gas).Fn+ dis_x(gen));//-count_reem_c;

//beta = 1.0/(stdh*sqrt(2.0));
//beta = 1.0/(stdh);
//dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nvh/(beta*2.0*sqrt(pi));
NnewH = floor(volume*(*gas).nh/(*gas).Fn+ dis_x(gen));


N_new = NnewC + NnewH;
/*
if( NnewC<2){
	printf("\n*** NnewC<1, = %d\n\n",NnewC);
}
if( NnewH<2){
	printf("\n*** NnewH<1, = %d\n\n",NnewH);
}
*/
Ntot = N_new + (*gas).N - NoutC - NoutH;

double *x2new, *U1new, *U2new, *U3new, *x2_oldnew;
x2new = (double *) malloc( Ntot * sizeof(double) );
U1new = (double *) malloc( Ntot * sizeof(double) );
U2new = (double *) malloc( Ntot * sizeof(double) );
U3new = (double *) malloc( Ntot * sizeof(double) );
x2_oldnew = (double *) malloc( Ntot * sizeof(double) );

j=0;
for(i=0; i<(*gas).N; i++){
	if((*flag)[i] != -1){
		x2new[j] = (*x2)[i];
		U1new[j] = (*U1)[i];
		U2new[j] = (*U2)[i];
		U3new[j] = (*U3)[i];
		x2_oldnew[j] = (*x2_old)[i];
		j++;
	}
}

for(i=0; i<NnewC; i++){
	  U1new[j] = dis_norm(gen)*stdc;
		U3new[j] = dis_norm(gen)*stdc;
		U2new[j] = dis_norm(gen)*stdc;
		x2new[j] = dis_x(gen)*tol;
		j++;
}

for(i=0; i<NnewH; i++){
	  U1new[j] = dis_norm(gen)*stdh;
		U3new[j] = dis_norm(gen)*stdh;
		U2new[j] = dis_norm(gen)*stdh;
		x2new[j] = (*box).Len[1] - tol*dis_x(gen) ;
		j++;
}

 free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
(*gas).N = Ntot;
*x2 = x2new;
*U1 = U1new;
*U2 = U2new;
*U3 = U3new;
*x2_old = x2_oldnew;

// Resize the arrays than will be initilized in other functions
// first free them
// then allocate memory
free(*x1);free(*x1_old); free(*n_ratio);
free(*flag); free(*index);
free(*xi_1); free(*xi_2); free(*xi_3); free(*xi_1f); free(*xi_2f); free(*xi_3f);
free(*Mp1); free(*Mp2); free(*Mp3);

*n_ratio = (int *) malloc( (*gas).N * sizeof(int) );
*index = (int *) malloc( (*gas).N * sizeof(int) );
*flag = (int *) malloc( (*gas).N * sizeof(int) );
*xi_1f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_2f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_3f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_1 = (double *) malloc( (*gas).N * sizeof(double) );
*xi_2 = (double *) malloc( (*gas).N * sizeof(double) );
*xi_3 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp1 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp2 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp3 = (double *) malloc( (*gas).N * sizeof(double) );


*x1 = (double *) malloc( Ntot * sizeof(double) );
*x1_old = (double *) malloc( Ntot * sizeof(double) );

// set n_ratio of particles to the cell's n_ratios
for(i=0; i<(*gas).N; i++){
	num = floor((*x2)[i]/(*box).delta_dim[1]);
	(*index)[i] = num;
	(*n_ratio)[i] = cells[num].n_ratio;
	if(num<0 || num>=(*box).N[1]){
		printf("-------> num = %d\n",num);
	}
}

}

void zoomed_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j,NoutC=0, NoutH=0;
double epsilon = 1e-15;
double Cratio, Hratio;
int NnewC,NnewH,N_new, Ntot;
int num;
double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
double stdh = sqrt((*gas).kb*(*gas).Th/(*gas).m);
double dt_remain,yw;
int count_reem_c=0, count_reem_h=0;
double pi = acos(-1.0);
// particles that leave the domain needed to be removed
for(i=0; i<(*gas).N; i++){
	if( (*x2)[i] < epsilon ){
		if( dis_x(gen) < 1.0-(*gas).sigc){
			if( dis_x(gen) < (*gas).alphac){
				yw = 0.0;
				dt_remain = fabs ( ((*x2)[i]-yw)/(*U2)[i] );

				(*U1)[i] = stdc * dis_norm( gen );
				(*U3)[i] = stdc * dis_norm( gen );
				(*U2)[i] = sqrt(2.0)*stdc*sqrt(-log(dis_x(gen)));
				(*x2)[i] = yw + epsilon+(*U2)[i]*dt_remain;
			}
			else{
				yw = 0.0;
				dt_remain = fabs ( ((*x2)[i]-yw)/(*U2)[i] );

				(*U2)[i] = -(*U2)[i];
				(*x2)[i] = yw + epsilon+(*U2)[i]*dt_remain;
			}
			//(*U2)[i] = - (*U2)[i];
			//(*x2)[i] = yw + epsilon+(*U2)[i]*dt_remain;
			(*flag)[i] = 1;
			count_reem_c ++;
		}
		else{
			NoutC++;
			(*flag)[i] = -1;
		}
	}
	else if( (*x2)[i] > (*box).Len[1]-epsilon ){
		if( dis_x(gen) < 1.0-(*gas).sigh ){
			if( dis_x(gen) < (*gas).alphah){
				yw = (*box).Len[1];
				dt_remain = fabs ( ((*x2)[i]-yw)/(*U2)[i] );

				(*U1)[i] = stdh * dis_norm( gen );
				(*U3)[i] = stdh * dis_norm( gen );
				(*U2)[i] = -sqrt(2.0)*stdh*sqrt(-log(dis_x(gen)));
				//(*x2)[i] = (*box).Len[1]-epsilon+(*U2)[i]*(*gas).delta_t*dis_x(gen);
				(*x2)[i] = yw-epsilon+(*U2)[i]*dt_remain;
			}
			else{
				yw = (*box).Len[1];
				dt_remain = fabs ( ((*x2)[i]-yw)/(*U2)[i] );

				(*U2)[i] = -(*U2)[i];
				(*x2)[i] = yw - epsilon+(*U2)[i]*dt_remain;
			}
			(*flag)[i] = 1;
			count_reem_h ++;
		}
		else{
			NoutH++;
			(*flag)[i] = -1;
		}
	}
	else{
		(*flag)[i] = 1;
	}
}
// I should remove those with negative flag.

// Now, how many new particles enter domain
double beta, dummy;
//nvc = 0.00287287/pow(3.405e-10,3.0);
//nvh = 0.00412534/pow(3.405e-10,3.0);


beta = 1.0/(stdc*sqrt(2.0));
//beta = 1.0/(stdc);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nvc/(beta*2.0*sqrt(pi));
//Cratio = 0.9;
//Hratio = 1.1;
//NnewC =  floor(Cratio*NoutC);
//NnewH =  floor(Hratio*NoutH);
NnewC = floor((*gas).sigc*dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;

beta = 1.0/(stdh*sqrt(2.0));
//beta = 1.0/(stdh);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nvh/(beta*2.0*sqrt(pi));
NnewH = floor((*gas).sigh*dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;

N_new = NnewC + NnewH;
if( NnewC<2){
	printf("\n*** NnewC<1, = %d\n\n",NnewC);
}
if( NnewH<2){
	printf("\n*** NnewH<1, = %d\n\n",NnewH);
}
Ntot = N_new + (*gas).N - NoutC - NoutH;

double *x2new, *U1new, *U2new, *U3new, *x2_oldnew;
x2new = (double *) malloc( Ntot * sizeof(double) );
U1new = (double *) malloc( Ntot * sizeof(double) );
U2new = (double *) malloc( Ntot * sizeof(double) );
U3new = (double *) malloc( Ntot * sizeof(double) );
x2_oldnew = (double *) malloc( Ntot * sizeof(double) );

j=0;
for(i=0; i<(*gas).N; i++){
	if((*flag)[i] != -1){
		x2new[j] = (*x2)[i];
		U1new[j] = (*U1)[i];
		U2new[j] = (*U2)[i];
		U3new[j] = (*U3)[i];
		x2_oldnew[j] = (*x2_old)[i];
		j++;
	}
}

for(i=0; i<NnewC; i++){
	  U1new[j] = dis_norm(gen)*stdc;
		U3new[j] = dis_norm(gen)*stdc;
		U2new[j] = sqrt(2.0)*stdc*sqrt(-log(dis_x(gen)));
		x2new[j] = U2new[j]*(*gas).delta_t*dis_x(gen);
		j++;
}

for(i=0; i<NnewH; i++){
	  U1new[j] = dis_norm(gen)*stdh;
		U3new[j] = dis_norm(gen)*stdh;
		U2new[j] = -sqrt(2.0)*stdh*sqrt(-log(dis_x(gen)));
		x2new[j] = (*box).Len[1] + U2new[j]*(*gas).delta_t*dis_x(gen);
		j++;
}

/*
double VMP, SC, QA, FS1, FS2, U, pratio,FS1m;
VMP = sqrt(2.0*(*gas).kb*(*gas).Tc/(*gas).m);
SC = 0.0;//(*box).U_wall_5/VMP;
QA = 3.0;
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));

for(i=0; i<NnewC; i++){
	  U1new[j] = dis_norm(gen)*stdc;
		U3new[j] = dis_norm(gen)*stdc;
		U = (10.0*dis_x(gen)*FS1-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		}
		U2new[j] = (U+SC)*VMP;
		x2new[j] = U2new[j]*(*gas).delta_t*dis_x(gen);
		j++;
}

VMP = sqrt(2.0*(*gas).kb*(*gas).Th/(*gas).m);
SC = 0.0;//(*box).U_wall_5/VMP;
QA = 3.0;
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));
for(i=0; i<NnewH; i++){
	  U1new[j] = dis_norm(gen)*stdh;
		U3new[j] = dis_norm(gen)*stdh;
		U = (10.0*dis_x(gen)*FS1m-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1m-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		}
		U2new[j] = (U+SC)*VMP;
		x2new[j] = (*box).Len[1] + U2new[j]*(*gas).delta_t*dis_x(gen);
		j++;
}
*/

// generating new particles that should enter due to evaporation

// combining all
 free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
(*gas).N = Ntot;
*x2 = x2new;
*U1 = U1new;
*U2 = U2new;
*U3 = U3new;
*x2_old = x2_oldnew;

// Resize the arrays than will be initilized in other functions
// first free them
// then allocate memory
free(*x1);free(*x1_old); free(*n_ratio);
free(*flag); free(*index);
free(*xi_1); free(*xi_2); free(*xi_3); free(*xi_1f); free(*xi_2f); free(*xi_3f);
free(*Mp1); free(*Mp2); free(*Mp3);

*n_ratio = (int *) malloc( (*gas).N * sizeof(int) );
*index = (int *) malloc( (*gas).N * sizeof(int) );
*flag = (int *) malloc( (*gas).N * sizeof(int) );
*xi_1f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_2f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_3f = (double *) malloc( (*gas).N * sizeof(double) );
*xi_1 = (double *) malloc( (*gas).N * sizeof(double) );
*xi_2 = (double *) malloc( (*gas).N * sizeof(double) );
*xi_3 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp1 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp2 = (double *) malloc( (*gas).N * sizeof(double) );
*Mp3 = (double *) malloc( (*gas).N * sizeof(double) );


*x1 = (double *) malloc( Ntot * sizeof(double) );
*x1_old = (double *) malloc( Ntot * sizeof(double) );

// set n_ratio of particles to the cell's n_ratios
for(i=0; i<(*gas).N; i++){
	num = floor((*x2)[i]/(*box).delta_dim[1]);
	(*index)[i] = num;
	(*n_ratio)[i] = cells[num].n_ratio;
	if(num<0 || num>=(*box).N[1]){
		printf("-------> num = %d\n",num);
	}
}

}
void inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag, int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color){

	int num, Ncells, id;
	int *numbering;
	double r;
	int *id_split_part;
	int N_new, N_old, N_light;
	double *x1n,*x2n, *U1n, *U2n, *U3n, *x1_oldn,*x2_oldn;
	double *x1new,*x2new, *U1new, *U2new, *U3new,*x1_oldnew, *x2_oldnew;
	double U[3], *U1mean, *U2mean, *U3mean;
	double *sig;
	int *weight;
	double h;
	int *numb_new, *numb_rem;
	double energy;
	int count, i, j;
	double ratio;
	int *colornew,*colorn;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis_r(0.0, 1.0);
	std::normal_distribution<> dis_norm(0.0,1.0);

	Ncells = (*box).N[2]*(*box).N[1]*(*box).N[0];
	numb_new = (int *) malloc(  Ncells* sizeof(int) );
	numb_rem = (int *) malloc(  Ncells* sizeof(int) );
        numbering = (int *) malloc(  Ncells* sizeof(int) );
	weight = (int *) malloc(  Ncells* sizeof(int) );
	sig = (double *) malloc(  Ncells* sizeof(double) );
	U1mean = (double *) malloc(  Ncells* sizeof(double) );
	U2mean = (double *) malloc(  Ncells* sizeof(double) );
	U3mean = (double *) malloc(  Ncells* sizeof(double) );
	for(id=0; id<(*gas).N; id++)
		(*flag)[id] = 0;
	for (num=0; num < Ncells; num++){
		numb_new[num] = 0;
		numb_rem[num] = 0;
		numbering[num] = 0;
		sig[num] = 0.0;
		weight[num] = 0;
		U1mean[num] = 0.0;
		U2mean[num] = 0.0;
		U3mean[num] = 0.0;
	}
	N_new = 0;
	N_light = 0;
	N_old = (*gas).N;
	ratio = 1.0-1.0/((*gas).n_ratio*1.0);
	for(id=0; id<(*gas).N; id++){
		num = (*index)[id];
		if((*n_ratio)[id] > cells[num].n_ratio){
			numb_new[num] = numb_new[num] + 1;
			N_new++;
			weight[num] = weight[num] + (*n_ratio)[id];
			(*flag)[id] = 1;
		}
		else if((*n_ratio)[id] < cells[num].n_ratio){
			// light particle in heavy cell
			N_light ++;
			r = dis_r(gen);
			if(r < ratio ){
				// delete
				(*flag)[id] = -1;
				N_old = N_old - 1;
			}
			else{
				// accept with the new weight
				(*n_ratio)[id] = cells[num].n_ratio;
				numb_rem[num] = numb_rem[num]+1;
				weight[num] = weight[num] + (*n_ratio)[id];
			}
		}
		else{
			// do nothing
			numb_rem[num] = numb_rem[num]+1;
			weight[num] = weight[num] + (*n_ratio)[id];
		}
	}
	if(N_new>0){
		for (num=0; num<Ncells; num++){
			if((*box).step != 1)
				free(cells[num].id_split_part);
			cells[num].id_split_part = (int *) malloc( numb_new[num] * sizeof(int) );
    		}

		for(id=0; id<(*gas).N; id++){
			if((*flag)[id] == 1){
				num = (*index)[id];
				cells[num].id_split_part[numbering[num]] = id;
				numbering[num]++;
			}
		}

		for(id=0; id<(*gas).N; id++){
			if((*flag)[id] != -1){
				num = (*index)[id];
				U1mean[num] = (*U1)[id]*(*n_ratio)[id];
				U2mean[num] = (*U2)[id]*(*n_ratio)[id];
				U3mean[num] = (*U3)[id]*(*n_ratio)[id];
			}
		}
		for(num=0; num<Ncells; num++){
			U1mean[num] = U1mean[num]/(1.0*weight[num]);
			U2mean[num] = U2mean[num]/(1.0*weight[num]);
			U3mean[num] = U3mean[num]/(1.0*weight[num]);
		}
		for(id=0; id<(*gas).N; id++){
			if((*flag)[id] != -1){
				num = (*index)[id];
				energy = ((*U1)[id]-U1mean[num])*((*U1)[id]-U1mean[num])
				        +((*U2)[id]-U2mean[num])*((*U2)[id]-U2mean[num])
					+((*U3)[id]-U3mean[num])*((*U3)[id]-U3mean[num]);
				sig[num] = sig[num]+energy*(*n_ratio)[id];
			}
		}
		for(num=0; num<Ncells; num++){
			sig[num] = sqrt( sig[num]/(1.0*weight[num]) );
		}

		N_new = ((*gas).n_ratio-1) * N_new;
		x1n = (double *) malloc( N_new * sizeof(double) );
		x2n = (double *) malloc( N_new * sizeof(double) );
		U1n = (double *) malloc( N_new * sizeof(double) );
		U2n = (double *) malloc( N_new * sizeof(double) );
		U3n = (double *) malloc( N_new * sizeof(double) );
		x1_oldn = (double *) malloc( N_new * sizeof(double) );
		x2_oldn = (double *) malloc( N_new * sizeof(double) );
		colorn = (int *) malloc( N_new * sizeof(int) );
		// generate new particles

		//double mean[3], st[3];
		count = 0;
		for(num=0; num<Ncells; num++){
			for(i=0; i<numb_new[num]; i++){
				id = cells[num].id_split_part[i];
				h = 1.06*sig[num]*pow(1.0*(numb_rem[num]+(*gas).n_ratio),-1.0/5.0);
				for(j=0; j<(*gas).n_ratio-1; j++){
					U1n[count] = (*U1)[id]+h*dis_norm(gen);
					U2n[count] = (*U2)[id]+h*dis_norm(gen);
					U3n[count] = (*U3)[id]+h*dis_norm(gen);
					x1n[count] = (*x1)[id];
					x1_oldn[count] = (*x1_old)[id];
					x2n[count] = (*x2)[id];
					x2_oldn[count] = (*x2_old)[id];
					colorn[count] = (*color)[id];
					count++;
				}
				//st[0] = standard_deviation(&U1n[count-((*gas).n_ratio-1)], (*gas).n_ratio-1, &mean[0]);
				//st[1] = standard_deviation(&U2n[count-((*gas).n_ratio-1)], (*gas).n_ratio-1, &mean[1]);
				//st[2] = standard_deviation(&U3n[count-((*gas).n_ratio-1)], (*gas).n_ratio-1, &mean[2]);
			}
		}
	}
	if(N_new>0 || N_old!=(*gas).N){
		int Ntot;
		Ntot = N_new+N_old;
		x1new = (double *) malloc( Ntot * sizeof(double) );
		x2new = (double *) malloc( Ntot * sizeof(double) );
		U1new = (double *) malloc( Ntot * sizeof(double) );
		U2new = (double *) malloc( Ntot * sizeof(double) );
		U3new = (double *) malloc( Ntot * sizeof(double) );
		x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
		x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
		colornew = (int *) malloc( Ntot * sizeof(int) );
		j=0;
		for(i=0; i<(*gas).N; i++){
			if((*flag)[i] != -1){
				x1new[j] = (*x1)[i];
				x2new[j] = (*x2)[i];
				U1new[j] = (*U1)[i];
				U2new[j] = (*U2)[i];
				U3new[j] = (*U3)[i];
				x1_oldnew[j] = (*x1_old)[i];
				x2_oldnew[j] = (*x2_old)[i];
				colornew[j] = (*color)[i];
				j++;
			}
		}
		for(i=0; i<N_new; i++){
			x1new[j] = x1n[i];
			x2new[j] = x2n[i];
			U1new[j] = U1n[i];
			U2new[j] = U2n[i];
			U3new[j] = U3n[i];
			x1_oldnew[j] = x1_oldn[i];
			x2_oldnew[j] = x2_oldn[i];
			colornew[j]  = colorn[i];
			j++;
		}
		free(*x1);free(*x2);free(*U1);free(*U2);free(*U3);free(*x1_old);free(*x2_old); free(*color);
		if(N_new>0){
			free(x1n);
			free(x2n);
			free(U1n);
			free(U2n);
			free(U3n);
			free(x1_oldn);
			free(x2_oldn);
			free(colorn);
		}
		(*gas).N = Ntot;
		*x1 = x1new;
		*x2 = x2new;
		*U1 = U1new;
		*U2 = U2new;
		*U3 = U3new;
		*x1_old = x1_oldnew;
		*x2_old = x2_oldnew;
		*color = colornew;

		// Resize the arrays than will be initilized in other functions
		// first free them
		// then allocate memory
		free(*flag); free(*index); free(*n_ratio);
		free(*xi_1); free(*xi_2); free(*xi_3); free(*xi_1f); free(*xi_2f); free(*xi_3f);
		free(*Mp1); free(*Mp2); free(*Mp3);

		*index = (int *) malloc( (*gas).N * sizeof(int) );
		*flag = (int *) malloc( (*gas).N * sizeof(int) );
		*n_ratio = (int *) malloc( (*gas).N * sizeof(int) );
		*xi_1f = (double *) malloc( (*gas).N * sizeof(double) );
		*xi_2f = (double *) malloc( (*gas).N * sizeof(double) );
		*xi_3f = (double *) malloc( (*gas).N * sizeof(double) );
		*xi_1 = (double *) malloc( (*gas).N * sizeof(double) );
		*xi_2 = (double *) malloc( (*gas).N * sizeof(double) );
		*xi_3 = (double *) malloc( (*gas).N * sizeof(double) );
 		*Mp1 = (double *) malloc( (*gas).N * sizeof(double) );
 		*Mp2 = (double *) malloc( (*gas).N * sizeof(double) );
 		*Mp3 = (double *) malloc( (*gas).N * sizeof(double) );

		// set n_ratio of particles to the cell's n_ratios
		for(id=0; id<(*gas).N; id++){
			num = floor((*x2)[id]/(*box).delta_dim[1]);
			(*index)[id] = num;
			(*n_ratio)[id] = cells[num].n_ratio;
			if(num<0 || num>=(*box).N[1]){
				printf("-------> num = %d\n",num);
			}
		}
	}
	// free local memory
	free(numb_new); free(numb_rem); free(numbering);
	free(weight); free(sig); free(U1mean);  free(U2mean); free(U3mean);
}
void apply_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *flag){
	double start, end;

          /////////////////////////////////////////////////////////////////////////////////
      //##########################       Boundary Condition     ##########################
          /////////////////////////////////////////////////////////////////////////////////
          //start = omp_get_wtime();
          if( (*box).problem == closed_box )
            BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
            //BC_periodic(U1, U2, U3, x1, x2,x3, &gas, &box);
          else if ( (*box).problem == couette_flow)
            //BC_thermal_wall_3D(U1, U2, U3, x1, x2, x3, gas, box, index, cells, flag, x1_old, x2_old, x3_old);
            BC_thermal_wall( U1, U2, U3, x1, x2, gas, box, index, cells, flag, x1_old, x2_old);
	  else if (  (*box).problem == vacuum){
		//BC_specular_reflection(U1, U2, U3, x1, x2,x3, &gas, &box);
		//BC_thermal_wall( U1, U2, U3, x1, x2, &gas, &box, index, cells, flag, x1_old, x2_old);
		BC_evaporation(U1,U2,U3,x1,x2,x3,gas,box,index,cells,flag,x1_old,x2_old,x3_old);
	  }
	  else if( (*box).problem == inverted){
		BC_inverted(U1,U2,U3,x1,x2,x3,gas,box,index,cells,x1_old,x2_old,x3_old);
	  }
	  else if((*box).problem == evaporation){
			//working BC_evaporation2(U1,U2,U3,x1,x2,x3,gas,box,index,cells,x1_old,x2_old,x3_old);
			//evap_BC(U1, U2, U3, x2, gas, box, cells);
			BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
	  }
		else if ( (*box).problem == shock){
			BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
			//shock_BC(U2, x2, T, gas, box, cells);
		}
		else if(    (*box).problem==wall   ){
			BC_wall(U1, U2, U3, x2, gas, box, index, cells, flag, x2_old);

		}
		else if((*box).problem == injection){
			BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
		}
	  //end = omp_get_wtime();
	  //if ( (*box).step % (*box).every == 0)
		//printf("BC took %e s\n", end - start);

}

void stream_particles(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index)
{
	double start, end;
      /////////////////////////////////////////////////////////////////////////////////
      //##########################            move the particles (streaming)       ##########################
      /////////////////////////////////////////////////////////////////////////////////
        //start = omp_get_wtime();
         if ( (*gas).model == "FP_dense" && (*box).problem != dcones){
		//Position_FP_ideal_gas(U1, U2, U3,x1,x2,x3, gas, box, cells, index, (*gas).delta_t);
            //Position_FP_new3(U1, U2, U3,x1,x2, gas, box, cells, index, (*gas).delta_t);
	    //
								Position_FP_ideal_gas(U1, U2, U3,x1,x2,x3, gas, box, cells, index, (*gas).delta_t);
								//
         }
         else if ( (*gas).model == "MD" && (*box).problem != dcones){
         }
	 else if( (*box).problem == flatnose || (*box).problem == dcones){
		stream_axisym(U1, U2, U3,x1,x2, gas, box, cells, index, (*gas).delta_t);
	 }
         else{
         	Position_FP_ideal_gas(U1, U2, U3,x1,x2,x3, gas, box, cells, index, (*gas).delta_t);
		//update_face_info(U1, U2, U3, x1, x2, x1_old, x2_old, gas, box, cells, index);
        }
        //end = omp_get_wtime();
        //if ( (*box).step % (*box).every == 0)
        //	printf("Streaming took %e s\n", end - start);

}
void BC_for_CBA(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *flag)
{
	if( (*gas).model == "CBA_VHS" || (*gas).model == "CBA_HS" ){
		double start, end;
      		//##########################   ONLY for CBA    Boundary Condition     ##########################
          	//start = omp_get_wtime();
          	if( (*box).problem == closed_box)
            		BC_specular_reflection(U1, U2, U3, x1, x2, x3, gas, box);
          	else if ( (*box).problem == couette_flow)
            		BC_thermal_wall_3D(U1, U2, U3, x1, x2, x3, gas, box, index, cells, flag, x1_old, x2_old, x3_old);
          	//end = omp_get_wtime();

          	//if ( (*box).step % (*box).every == 0)
      			//printf("BC for CBA done in %e sec\n", end - start);
      	}
}
void velocity_update(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, struct Collision *CBA, double *xi_1f, double *xi_2f, double *xi_3f, double *Mp1, double *Mp2, double *Mp3, double *F1, double *F2, double *F3, double *rho, double *T, int* color){
	int i;
	double start, end;
	//start = omp_get_wtime();

	// Collision
   	if( (*box).step >= (*box).init_step){
          	if ( (*gas).model == "DSMC_VHS" ){
            		Collision_DSMC_VHS(U1,U2,U3,x1,x2, gas, box, cells, CBA);
      	   	}
          	else if ( (*gas).model == "CBA_HS" )
            		Collision_CBA_HS(U1,U2,U3,x1,x2, gas, box, cells, CBA);
						else if ( (*gas).model == "SPH" )
								SPH_simple_hybrid(x2, U2, rho, T, gas, box, cells, index, color);//SPH_simple(x2, U2, rho, T, gas, box, cells, index);
							//	SPH_velocity_temperature_update(x2, U2, rho, T, gas, box, cells, index);
          	else if( (*gas).model == "CBA_VHS" )
            		Collision_CBA_VHS(U1,U2,U3,x1,x2, gas, box, cells, CBA);
          	else if( (*gas).model == "ESMC"){
			if((*box).problem == dcones)
				Collision_ESMC_dcones(U1,U2,U3,x1,x2, gas, box, cells, index);
			else{
      				//Collision_ESMC_new(U1,U2,U3,x1,x2, gas, box, cells, index);
							//Collision_ESMC_new4(U1,U2,U3,x1,x2, gas, box, cells, index);
							Collision_ESMC_new3(U1,U2,U3,x1,x2, gas, box, cells, index);
        	 	   //	Collision_ESMC_1D(U1,U2,U3,x1,x2,x3, gas, box, cells, CBA, index);
			}
          	}
          	else if ( (*gas).model == "FP_dense" || (*gas).model == "FP_ideal_gas"){
            		Velocity_FP_parallel(U1,U2,U3, gas, cells, box, index, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3);
 			//	Velocity_FP_linear(U1,U2,U3, &gas, cells, &box, index, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3);
          	}
          	else if ( (*gas).model == "FP_linear")
            		Velocity_FP_linear(U1,U2,U3, gas, cells, box, index, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3);
						else if( (*gas).model == "Hybrid" ){
							SPH_simple_hybrid(x2, U2, rho, T, gas, box, cells, index, color);
							Collision_DSMC_VHS_hybrid(U1, U2, U3, x1, x2, gas, box, cells, CBA, color);
						}
          	else if ( (*gas).model == "MD" ){
            		//    Collision_DSMC_VHS(U1,U2,U3,x1,x2, &gas, &box, cells, &CBA);
			//    MD(U1, U2, U3, x1, x2, x3, &gas, &box, cells, index,F1, F2, F3, F1_old, F2_old, F3_old);
	    		MD_new(U1, U2, U3, x1, x2, x3, gas, box, cells, index,F1, F2, F3);
            		// check new values
            		for(i=0; i<(*gas).N; i++){
                  		if(U1[i] != U1[i]){
                        		printf("U1[%d] = %e in step = %d\n",i, U1[i] ,(*box).step);
                        		getchar();
                  		}
                  		if(U2[i] != U2[i]){
                        		printf("U2[%d] = %e in step = %d\n",i, U2[i] , (*box).step);
                        		getchar();
                  		}
                  		if(U3[i] != U3[i]){
                        		printf("U3[%d] = %e in step = %d\n",i, U3[i] , (*box).step);
                        		getchar();
                  		}
            		}
          	}
          	else if ( (*gas).model == "none"){
          	}
          	else{
            		printf("Error: Given model name does not have correspondence\n");
            		exit(1);
          	}
	}
	//end = omp_get_wtime();
	// Print time it took
	//if ( (*box).step % (*box).every == 0)
	//	printf("Collision + thermostat took %e s\n", end - start);

	// Long Range Interaction
	//start = omp_get_wtime();
	double dummyF;
	double pi = acos(-1);
	int num, id;
	/* long range interaction:
		1 = direct
		2 = density expansion
		3 = screened-Poisson
		0 or anything>3 = no attraction
	*/
	if((*gas).long_range == 1){
		if((*box).problem == dcones)
			dcones_vlasov_integral(U1, U2, gas, box, cells, index);
		else
			vlasov_taking_the_integral(U2, gas, box, cells, index, x1, x2, x3);
		if ((*box).step % (*box).every == 0)
			std::cout << "attraction with direct method!\n";
	}
	else if( (*gas).long_range == 2){
		// This is based on gradient of denisty (from Henning's derivation)
      		for(id=0; id< (*gas).N; id++){
	     		num = index[id];
             	if(num == 0)
             		dummyF = cells[num].dn2;//
								//dummyF = (cells[num+1].n-cells[num].n)/(2.0*(*box).delta_dim[1]);
             	else if(num == (*box).N[1]-1)
               		dummyF = cells[num].dn2;//
									//dummyF = (cells[num].n-cells[num-1].n)/(2.0*(*box).delta_dim[1]);
             	else{
								//dummyF = (cells[num+1].n-cells[num-1].n)/(2.0*(*box).delta_dim[1]);
             		dummyF = cells[num].dn2;//+gas.eqdist*gas.eqdist/2.0*cells[num].dn3;
								//dummyF = -(cells[num+1].n_smooth-cells[num-1].n_smooth)/(2.0*box.delta_dim[1]);
	      	}
              dummyF = dummyF*4.0*pi*pow((*gas).sigma*pow(2.0,1.0/6.0),3.0)*(*gas).phi/(*gas).m*2.0;
	      	if(num != 0 && num != (*box).N[1]-1)
						U2[id] = U2[id] - dummyF*(*gas).delta_t;
      	   	}
	   	if ((*box).step % (*box).every == 0)
	   		std::cout << "attraction with density expansion!\n";
		}
		else if((*gas).long_range == 3){
			if((*box).problem == injection){
				//compute_SP_2D(gas, box, cells,x1,x2,U1,U2,index);
			}
			else{
				poisson_solver( U2, x1, x2, x3, gas, box, cells, index);
			}
		if ((*box).step % (*box).every == 0)
				std::cout << "attraction with screened-Poisson!\n";
		}
	else{

		if( (*box).step % (*box).every == 0)
			std::cout << "no attraction!\n";
	}
//        end = omp_get_wtime();
//	if ( (*box).step % (*box).every == 0)
//		printf("Attraction took %e s\n", end - start);

}
void thermostat(double *U1,double *U2,double *U3, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  int i, num, id;
  double lamb;
  double therm_len = 4.0*(*gas).sigma;
  //double therm_len = 1.2*(*gas).Lc;
	if((*box).problem == injection){
		/*
		for(i=0; i<(*gas).N; i++){
				num = index[i];
				if(cells[num].T > 0.0){
					lamb = sqrt((*gas).T/cells[num].T);
					U1[i] = U1[i]*lamb;
					U2[i] = U2[i]*lamb;
					U3[i] = U3[i]*lamb;
			}
		}
		*/
	}
  else if ( (*box).problem == vacuum){
		/*
	for(i=0; i<(*gas).N; i++){
		if( x2[i]>(*box).Len[1]/2.0 - 2.0*(*gas).sigma && x2[i]<(*box).Len[1]/2.0 + 2.0*(*gas).sigma ){
			num = index[i];
			if(cells[num].T > 0.0){
				lamb = sqrt((*gas).T/cells[num].T);
				U1[i] = U1[i]*lamb;
				U2[i] = U2[i]*lamb;
				U3[i] = U3[i]*lamb;
			}
		}
	}
	*/
  }

  else if ( (*box).problem == evaporation ){
	/*
	if((*box).step % 5000 == 0){
		printf("!!!    gas.T=%lf\n", (*gas).T);
		for(i=0; i<(*gas).N; i++){
			if( fabs( x2[i]-(*box).Len[1]/2.0 ) <  (*box).ghost/2.0 ){
				//printf(" (*box).ghost = %e\n", (*box).ghost );
				num = index[i];
				if(cells[num].T > 0.0){
					lamb = sqrt((*gas).T/cells[num].T);
					U1[i] = lamb*(U1[i]-cells[num].U_space[0]) + cells[num].U_space[0];
					U2[i] = lamb*(U2[i]-cells[num].U_space[1]) + cells[num].U_space[1];
					U3[i] = lamb*(U3[i]-cells[num].U_space[2]) + cells[num].U_space[2];
					//U1[i] = U1[i]*lamb;
					//U2[i] = U2[i]*lamb;
					//U3[i] = U3[i]*lamb;
				}
			}
		}

	}
	*/
  }

//&& (*box).step % 500 == 0
  else if( (*box).problem ==inverted ){
 	std::random_device rd;
  	std::mt19937 gen(rd());
	int stt = floor((*box).ghost/(*box).delta_dim[1])+1;
	double stdc = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
	double stdh = sqrt((*gas).kb*(*gas).Th/(*gas).m);
  	std::normal_distribution<> dis_u1(0.0,1.0);
  	std::normal_distribution<> dis_v1(0.0,1.0);
  	std::normal_distribution<> dis_w1(0.0,1.0);


	//int ncp = floor( 2.0*(*box).ghost/(*box).delta_dim[1] ) + 1;
	//int nhm = floor( ((*box).Len[1] - 2.0*(*box).ghost)/(*box).delta_dim[1] ) - 1;
	//(*gas).Tc = cells[ncp].T;
	//(*gas).Th = cells[nhm].T;
	//therm_len = (*gas).Lc;
	//therm_len = (*box).therm;
	//if((*box).step % 1000 == 0){
	    for(num=0; num< (*box).N[1]; num++){

		if( cells[num].cell_center[1] < (*box).thermc ){

			/*
		  //if( cells[num].cell_center[1] < (*box).ghost){
			if( cells[num].cell_center[1] < -1.0){
				lamb = sqrt((*gas).Tc/cells[num].T);
				for(id=0; id<cells[num].num_inside; id++){
					i = cells[num].indices_inside[id];

					//U1[i] = lamb*(U1[i]-cells[num].U_space[0]);// + cells[num].U_space[0];
					//U2[i] = lamb*(U2[i]-cells[num].U_space[1]);// + cells[num].U_space[1];
					//U3[i] = lamb*(U3[i]-cells[num].U_space[2]);// + cells[num].U_space[2];

					U1[i] = stdc*dis_u1(gen);//+ cells[stt].U_space[0];
					U2[i] = stdc*dis_v1(gen);//+ cells[stt].U_space[1];
					U3[i] = stdc*dis_w1(gen);//+ cells[stt].U_space[2];

					//if( (*box).step > 500 ){
					//	U1[i] += cells[stt].U_space[0];
					//	U2[i] += cells[stt].U_space[1];
				//		U3[i] += cells[stt].U_space[2];
				//	}
				//
				}

			}
			*/

		//	else{
		//		if((*box).step % (*box).reset_every) == 0){
			  if(cells[num].T > 0.0){
					lamb = sqrt((*gas).Tc/cells[num].T);
					for(id=0; id<cells[num].num_inside; id++){
						i = cells[num].indices_inside[id];
						U1[i] = lamb*(U1[i]-cells[num].U_space[0]) + cells[num].U_space[0];
						U2[i] = lamb*(U2[i]-cells[num].U_space[1]) + cells[num].U_space[1];
						U3[i] = lamb*(U3[i]-cells[num].U_space[2]) + cells[num].U_space[2];
					}
			 	}
	//		}
	//	}


		}
		else if( cells[num].cell_center[1] >(*box).Len[1]- (*box).thermh   ){

			 /*
			 //if( cells[num].cell_center[1] >(*box).Len[1]-(*box).ghost ){
  		 if( cells[num].cell_center[1] >(*box).Len[1]+1.0 ){
				 lamb = sqrt((*gas).Th/cells[num].T);
			 	for(id=0; id<cells[num].num_inside; id++){
					i = cells[num].indices_inside[id];

					//U1[i] = lamb*(U1[i]-cells[num].U_space[0]);// + cells[num].U_space[0];
					//U2[i] = lamb*(U2[i]-cells[num].U_space[1]);// + cells[num].U_space[1];
					//U3[i] = lamb*(U3[i]-cells[num].U_space[2]);// + cells[num].U_space[2];

					U1[i] = stdh*dis_u1(gen);//+ cells[(*box).N[1]-stt-1].U_space[0];
					U2[i] = stdh*dis_v1(gen);//+ cells[(*box).N[1]-stt-1].U_space[1];
					U3[i] = stdh*dis_w1(gen);//+ cells[(*box).N[1]-stt-1].U_space[2];

					//if( (*box).step > 500 ){
					//	U1[i] += cells[(*box).N[1]-stt-1].U_space[0];
					//	U2[i] += cells[(*box).N[1]-stt-1].U_space[1];
					//	U3[i] += cells[(*box).N[1]-stt-1].U_space[2];
					//}
					//
				}

			 }
			 */

	//		else{
	//			if((*box).step % (*box).reset_every) == 0){
			   if(cells[num].T > 0.0){
					 lamb = sqrt((*gas).Th/cells[num].T);
					 for(id=0; id<cells[num].num_inside; id++){
						 i = cells[num].indices_inside[id];
						 U1[i] = lamb*(U1[i]-cells[num].U_space[0]) + cells[num].U_space[0];
						 U2[i] = lamb*(U2[i]-cells[num].U_space[1]) + cells[num].U_space[1];
						 U3[i] = lamb*(U3[i]-cells[num].U_space[2]) + cells[num].U_space[2];
					 }
			   }
	//		}
	//	}

		}
	   }
	//}
  }
/*
  else if( (*box).problem ==inverted){
	for(i=0; i<(*gas).N; i++){
		if( x2[i]>(*gas).Lc - 8.0*(*gas).sigma && x2[i]<(*gas).Lc + 8.0*(*gas).sigma ){
		//if( x2[i]<(*gas).Lc + 8.0*(*gas).sigma ){
			num = index[i];
			if(cells[num].T > 0.0){
				lamb = sqrt((*gas).Tc/cells[num].T);
				U1[i] = U1[i]*lamb;
				U2[i] = U2[i]*lamb;
				U3[i] = U3[i]*lamb;
			}
		}
		else if( x2[i] > 1.5*(*gas).Lc+(*gas).Lv+0.5*(*gas).Lh - 8.0*(*gas).sigma && x2[i]<1.5*(*gas).Lc+(*gas).Lv+0.5*(*gas).Lh + 8.0*(*gas).sigma ){
		//else if( x2[i] > 1.5*(*gas).Lc+(*gas).Lv+0.5*(*gas).Lh - 8.0*(*gas).sigma ){
			num = index[i];
			if(cells[num].T > 0.0){
				lamb = sqrt((*gas).Th/cells[num].T);
				U1[i] = U1[i]*lamb;
				U2[i] = U2[i]*lamb;
				U3[i] = U3[i]*lamb;
			}
		}
	}
  }
*/
}

void stream_axisym(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt){
/*
it's a rotation around the axis -x (in r,theta,x cylindrical coord.) in the direction of -theta
u' = v*cos(-theta)-w*sin(-theta) =  v*cos(theta)+w*sin(theta)
v' = v*sin(-theta)+w*cos(-theta) = -v*sin(theta)+w*cos(theta)
where theta is
cos(theta) = (YI+VI*dt)/Y
sin(theta) = (WI*dt)/Y
YI is initial radius. ZI is zero since everytime position is mapped back to z=0
VI and WI are initial velocities of particle in caresian coord.
output is rotated position Y, velocities U1 and U2
*/
	int ii;
	double XI, YI;
	double VI, WI;
	double X,Y;
	double DX,DY,DZ;
//#pragma omp parallel for  private(energy, num, U_mean)
	for(ii=0; ii<(*gas).N; ii++){
		XI = x1[ii];
		YI = x2[ii];
		VI = U2[ii];
		WI = U3[ii];
		DX = U1[ii]*dt;
		DY = VI*dt;
		DZ = WI*dt;
		X = XI + DX;
		Y = sqrt( (YI+DY)*(YI+DY)+DZ*DZ );
		U2[ii] = (VI*(YI+DY)+WI*DZ)/Y;
		U3[ii] = (WI*(YI+DY)-VI*DZ)/Y;
		x1[ii] = X;
		x2[ii] = Y;
	}
}

double increase_collision_rate(double n, struct GAS *gas)
{
  double  Y;
  double eta;
  eta = n*(*gas).b/4.0;
  Y = 1.0 + 0.625*n*(*gas).b + 0.2869*pow( n*(*gas).b,2.0 ) + 0.115*pow( n*(*gas).b,3.0 );
  //Y = (1.0-eta/2.0)/( pow(1.0-eta,3.0) );
  return Y;
	//return 1.0;
}


void cell_poly_velocity(double *U1,double *U2,double *U3, struct CELLS *cells,  struct BOX *box, int step, struct GAS *gas, int *index)
{
  int i, num;
	int m,n,l,id;
	double energy;
  int num_size = (*box).N[2]*(*box).N[1]*(*box).N[0];
  for(i=0; i<(*gas).N; i++){
    num = index[i];

    if(cells[num].in==1)
	    cells[num].sum_weight += 1.0;

    cells[num].sum_M[0] +=  U1[i];
    cells[num].sum_M[1] +=  U2[i];
    cells[num].sum_M[2] +=  U3[i];

    cells[num].sum_MM[0] += U1[i]*U1[i];
    cells[num].sum_MM[1] += U1[i]*U2[i];
    cells[num].sum_MM[2] += U1[i]*U3[i];
    cells[num].sum_MM[3] += U2[i]*U2[i];
    cells[num].sum_MM[4] += U2[i]*U3[i];
    cells[num].sum_MM[5] += U3[i]*U3[i];

		id = 0;
		for(m=0; m<4; m++){
			for(n=0; n<4-m; n++){
				l = 3 - m - n;
				cells[num].sum_MMM[id] += pow(U1[i],m)*pow(U2[i],n)*pow(U3[i],l);
				id++;
			}
		}
		/*
    cells[num].sum_MMM[0] +=  U1[i]*U1[i]*U1[i];//111
    cells[num].sum_MMM[1] +=  U2[i]*U2[i]*U1[i];//221
    cells[num].sum_MMM[2] +=  U3[i]*U3[i]*U1[i];//331

    cells[num].sum_MMM[3] +=  U1[i]*U1[i]*U2[i];//112
    cells[num].sum_MMM[4] +=  U2[i]*U2[i]*U2[i];//222
    cells[num].sum_MMM[5] +=  U3[i]*U3[i]*U2[i];//332
		*/
		energy = U1[i]*U1[i] + U2[i]*U2[i] + U3[i]*U3[i];

		cells[num].sum_MMMM[0] +=  U1[i]*U1[i]*energy;//11
		cells[num].sum_MMMM[1] +=  U1[i]*U2[i]*energy;//12
		cells[num].sum_MMMM[2] +=  U1[i]*U3[i]*energy;//13

		cells[num].sum_MMMM[3] +=  U2[i]*U2[i]*energy;//22
		cells[num].sum_MMMM[4] +=  U2[i]*U3[i]*energy;//23
		cells[num].sum_MMMM[5] +=  U3[i]*U3[i]*energy;//33
  }
}

void scale_mu(double T_want, double T, double *mu_want, double mu, struct GAS *gas)
{
  *mu_want = mu*pow(T_want/T,(*gas).visp);
}

void update_index_number_inside(double *x1, double *x2,  double *U1, double *U2, double *U3, double *Mp1, double *Mp2, double *Mp3, int *index, struct CELLS *cells,  struct BOX *box, struct GAS *gas){
int i,j,k, num, id;
double energy, vol;
  int *numbering = (int *) malloc(  (*box).N[2]*(*box).N[1]*(*box).N[0]* sizeof(int) );
  for (num=0; num < (*box).N[0]*(*box).N[1]*(*box).N[2]; num++){
	  cells[num].num_inside = 0;
	  numbering[num] = 0;}


  for (i=0; i < (*gas).N; i++){
		if( (*box).direction[0] == 1 )
    	j = floor(x1[i]/(*box).delta_dim[0]);
		else
			j=0;
    k = floor(x2[i]/(*box).delta_dim[1]);
    //l = floor(x3[i]/(*box).delta_dim[2]);
    if(j<0)
      printf("OOPS, particle is out from left, x1=%e x1<0\n",x1[i]);
    if(j>(*box).N[0]-1)
      printf("OOPS, particle is out from right, x1=%e x1>L\n",x1[i]);
    if(k<0)
      printf("OOPS, particle is out from down, x2=%e x2<0\n",x2[i]);
    if(k>(*box).N[1]-1)
      printf("OOPS, particle is out from up, x2=%e x2>L\n",x2[i]);

    num = k*(*box).N[0] + j;
    index[i] = num;
    cells[num].num_inside = cells[num].num_inside+1;
  }

  for (j=0; j < (*box).N[0]; j++)
    for (k=0; k < (*box).N[1]; k++){
	  num =  k*(*box).N[0] +j;
	  free(cells[num].indices_inside);
	  cells[num].indices_inside = (int *) malloc( cells[num].num_inside * sizeof(int) );
	}

  for (i=0; i < (*gas).N; i++){
		if(	(*box).direction[0]==1 )
    	j = floor(x1[i]/(*box).delta_dim[0]);
		else
			j=0;
    k = floor(x2[i]/(*box).delta_dim[1]);
    num = k*(*box).N[0] + j;
    cells[num].indices_inside[ numbering[num] ] = i;
    numbering[num]++;
  }


  for(id=0; id<(*gas).N; id++){
      num = index[id];
      Mp1[id] = U1[id]-cells[num].U_space[0];
      Mp2[id] = U2[id]-cells[num].U_space[1];
      Mp3[id] = U3[id]-cells[num].U_space[2];
  }

}
void cell_update_post(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *Mp1,double *Mp2,double *Mp3)
{
  double vol;
  int num, i, j, k;
  int *numbering = (int *) malloc(  (*box).N[2]*(*box).N[1]*(*box).N[0]* sizeof(int) );
  for (num=0; num < (*box).N[0]*(*box).N[1]*(*box).N[2]; num++){
	  cells[num].num_inside = 0;
	  numbering[num] = 0;}
  for (i=0; i < (*gas).N; i++){
    j = floor(x1[i]/(*box).delta_dim[0]);
    k = floor(x2[i]/(*box).delta_dim[1]);
    if(j<0)
      printf("OOPS, particle is out from left, x1=%e x1<0\n",x1[i]);
    if(j>(*box).N[0]-1)
      printf("OOPS, particle is out from right, x1=%e x1>L\n",x1[i]);
    if(k<0)
      printf("OOPS, particle is out from down, x2=%e x2<0\n",x2[i]);
    if(k>(*box).N[1]-1)
      printf("OOPS, particle is out from up, x2=%e x2>L\n",x2[i]);
    num =  k*(*box).N[0] + j;
    index[i] = num;
    cells[num].num_inside = cells[num].num_inside+1;
  }

  for (j=0; j < (*box).N[0]; j++)
    for (k=0; k < (*box).N[1]; k++){
	  num =  k*(*box).N[0] +j;
	  free(cells[num].indices_inside);
	  cells[num].indices_inside = (int *) malloc( cells[num].num_inside * sizeof(int) );
	}

  for (i=0; i < (*gas).N; i++){
    j = floor(x1[i]/(*box).delta_dim[0]);
    k = floor(x2[i]/(*box).delta_dim[1]);
    num = k*(*box).N[0] + j;
    cells[num].indices_inside[ numbering[num] ] = i;
    numbering[num]++;
  }

  for(num=0; num< (*box).N[2]*(*box).N[1]*(*box).N[0]; num++){
    	cells[num].U_space[0] = 0.0;
	cells[num].U_space[1] = 0.0;
	cells[num].U_space[2] = 0.0;
  }
  int id;
//   #pragma omp parallel for  private(i, id, num)
  for(num=0; num< (*box).N[2]*(*box).N[1]*(*box).N[0]; num++){
    for(i=0; i<cells[num].num_inside; i++){
      id = cells[num].indices_inside[i];
      cells[num].U_space[0] = cells[num].U_space[0]+U1[id];
      cells[num].U_space[1] = cells[num].U_space[1]+U2[id];
      cells[num].U_space[2] = cells[num].U_space[2]+U3[id];
    }
  }
  for(num=0; num<(*box).N[2]*(*box).N[1]*(*box).N[0]; num++){
    cells[num].U_space[0] = cells[num].U_space[0]/(1.0*cells[num].num_inside);
    cells[num].U_space[1] = cells[num].U_space[1]/(1.0*cells[num].num_inside);
    cells[num].U_space[2] = cells[num].U_space[2]/(1.0*cells[num].num_inside);
    vol = (cells[num].dim[1]-cells[num].dim[0])*(cells[num].dim[3]-cells[num].dim[2])*(cells[num].dim[5]-cells[num].dim[4]);
    cells[num].n = (*gas).Fn*(1.0*cells[num].num_inside)/vol;
  }
 free(numbering);
}

int dcone_index(double x, double y, struct CELLS *cells,  struct BOX *box, int *j, int *k){
	double dx[3];
//	int j,k;
	int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
	dx[0] = (*box).Lx[0]/(1.0*(*box).Nx[0]);
	dx[1] = (*box).Lx[1]/(1.0*(*box).Nx[1]);
	dx[2] = (*box).Lx[2]/(1.0*(*box).Nx[2]);
	if(x<(*box).Lx[0]){
	    	*j = floor(x/dx[0]);
    		*k = floor( log(1.0-(1.0-(*box).s)*y/(*box).dy0) / log((*box).s) );
	}
	else if(x<(*box).Lx[0]+(*box).Lx[1]){
		*j = (*box).Nx[0]+floor((x-(*box).Lx[0])/dx[1]);
		*k = floor( log(1.0-(1.0-(*box).s)*(y-(x-(*box).Lx[0])*tan((*box).alpha[0]))/(*box).dy0) / log((*box).s) );
	}
	else{
		*j = (*box).Nx[0]+(*box).Nx[1]+floor((x-(*box).Lx[0]-(*box).Lx[1])/dx[2]);
		*k = floor( log(1.0-(1.0-(*box).s)*(y-(*box).Lx[1]*tan((*box).alpha[0])-(x-(*box).Lx[0]-(*box).Lx[1])*tan((*box).alpha[1]))/(*box).dy0) / log((*box).s) );
	}
	return  (*k)*Nx + (*j);
}
int flatnose_index(double x, double y, struct CELLS *cells,  struct BOX *box, int *j, int *k){
	double dx[2],dy[2];
//	int j,k;
	int Nx = (*box).N1[0]+(*box).N1[1];
	dx[0] = (*box).L1[0]/(1.0*(*box).N1[0]);
	dx[1] = (*box).L1[1]/(1.0*(*box).N1[1]);

	dy[0] = (*box).L2[0]/(1.0*(*box).N2[0]);
	dy[1] = (*box).L2[1]/(1.0*(*box).N2[1]);
	if(x<(*box).L1[0])
		*j = floor( log(1.0-(1.0-(*box).s1)*x/(*box).dx0) / log((*box).s1) );
	else
		//*j = (*box).N1[0] + floor((x-(*box).L1[0])/dx[1]);
		*j = (*box).N1[0] +floor( log(1.0-(-1.0+(*box).s1)*(x-(*box).L1[0])/(*box).dx1) / log(2.0-(*box).s1) );
  	if(y<(*box).L2[0])
	    //	*k = floor(y/dy[0]);
		*k = floor( log(1.0-(-1.0+(*box).s2)*(y)/(*box).dy0) / log(2.0-(*box).s2) );
	else
		*k = (*box).N2[0] + floor( log(1.0-(1.0-(*box).s2)*(y-(*box).L2[0])/(*box).dy1) / log((*box).s2) );
	return  (*k)*Nx + (*j);
}
void cell_update(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double ** T, double * rho, int *color)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  int id;
  double vol;
  int num, i, j, k, num2, kk,jj;
  //double   *E_kin;
  double energy;
  //E_kin = (double *) malloc(  (*box).N[2]*(*box).N[1]*(*box).N[0]* sizeof(double) );
//   find_inside_cells(x1, x2, x3, (*gas).N, box, cells, index);
 // int number_insides[(*box).N[0]*(*box).N[1]*(*box).N[2]];




  double Mp1, Mp2, Mp3;
  int sgn;
  double r;
  int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
  int Ny;
	double *dx = (double *) malloc(  3* sizeof(double) );

  dx[0] = (*box).Lx[0]/(1.0*(*box).Nx[0]);
  dx[1] = (*box).Lx[1]/(1.0*(*box).Nx[1]);
  dx[2] = (*box).Lx[2]/(1.0*(*box).Nx[2]);
int Ncells;
if((*box).problem == dcones)
	Ncells = (*box).Ny*Nx;
if((*box).problem == flatnose){
	Nx = (*box).N1[0]+(*box).N1[1];
	Ny = (*box).N2[0]+(*box).N2[1];
	Ncells = Ny*Nx;
}
else
	Ncells = (*box).N[2]*(*box).N[1]*(*box).N[0];

    int *numbering = (int *) malloc(  Ncells* sizeof(int) );
  for (num=0; num < Ncells; num++){
		cells[num].num_inside = 0;
		numbering[num] = 0;
		for(i=0; i<4; i++){
			cells[num].n_n[i] = 0.0;
		}
  }

  if((*box).problem == dcones){
	for (i=0; i < (*gas).N; i++){
		if(x1[i]<(*box).Lx[0]){
		    	j = floor(x1[i]/dx[0]);
    			k = floor( log(1-(1-(*box).s)*x2[i]/(*box).dy0) / log((*box).s) );
			//if(k>0)
			//	printf("x2[%d] = %e\n", i, x2[i]);
		}
		else if(x1[i]<(*box).Lx[0]+(*box).Lx[1]){
			j = (*box).Nx[0]+floor((x1[i]-(*box).Lx[0])/dx[1]);
			k = floor( log(1-(1-(*box).s)*(x2[i]-(x1[i]-(*box).Lx[0])*tan((*box).alpha[0]))/(*box).dy0) / log((*box).s) );
		}
		else{
			j = (*box).Nx[0]+(*box).Nx[1]+floor((x1[i]-(*box).Lx[0]-(*box).Lx[1])/dx[2]);
			k = floor( log(1-(1-(*box).s)*(x2[i]-(*box).Lx[1]*tan((*box).alpha[0])-(x1[i]-(*box).Lx[0]-(*box).Lx[1])*tan((*box).alpha[1]))/(*box).dy0) / log((*box).s) );
		}
		num =  k*Nx + j;
    		if(j<0)
      			printf("OOPS, particle is out from left, x1=%e x1<0\n",x1[i]);
    		if(j>(*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2]-1)
      			printf("OOPS, particle is out from right, x1=%e x1>L=%e\n",x1[i],(*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2]);
    		if(k<0)
      			printf("OOPS, particle is out from down, x2=%e x2<0\n",x2[i]);
    		if(k>(*box).Ny-1)
      			printf("OOPS, particle is out from up, x2=%e x2>L\n",x2[i]);
    		index[i] = num;
		cells[num].num_inside = cells[num].num_inside+1;
	}
  }
  else if((*box).problem == flatnose){
	for (i=0; i < (*gas).N; i++){
		num =  flatnose_index(x1[i], x2[i], cells, box, &j, &k);
    		if(j<0)
      			printf("OOPS, particle is out from left, x1=%e x1<0\n",x1[i]);
    		if(j>Nx-1)
      			printf("OOPS, particle is out from right, x1=%e x1>L\n",x1[i]);
    		if(k<0)
      			printf("OOPS, particle is out from down, x2=%e x2<0\n",x2[i]);
    		if(k>Ny-1)
      			printf("OOPS, particle is out from up, x2=%e x2>L\n",x2[i]);
		if(k<(*box).N2[0] && j>=(*box).N1[0])
			printf("OOPS, particle is out from right down, x1=%e x2=%e x2>L\n",x1[i],x2[i]);
    		index[i] = num;
		cells[num].num_inside = cells[num].num_inside+1;
	}
  }
  else{
  	for (i=0; i < (*gas).N; i++){
	    	if( (*box).direction[0] == 1 )
			j = floor(x1[i]/(*box).delta_dim[0]);
		else
			j = 0;
    		k = floor(x2[i]/(*box).delta_dim[1]);
    		if(k>0 && k<(*box).N[1]-1){
    			if(fabs(x2[i]-(*box).delta_dim[1]*k)<1e-15 || fabs(x2[i]-(*box).delta_dim[1]*(k+1))<1e-15){
				r = dis_r(gen);
				if(r>0.5)
					sgn = 1;
				else
					sgn = -1;
				k = floor( (x2[i]+1e-14*sgn) /(*box).delta_dim[1]);
     			}
    		}
    		if(j<0){
      			printf("OOPS, particle is out from left, j=%d x1<0\n",j);
			exit(1);
		}
    		if(j>(*box).N[0]-1){
      			printf("OOPS, particle is out from right, j=%d x1>L\n",j);
			exit(1);
		}
    		if(k<0){
      			printf("OOPS, particle is out from down, x2=%e x2<0\n",x2[i]);
			exit(1);
		}
    		if(k>(*box).N[1]-1){
      			printf("OOPS, particle is out from up, x2=%e x2>L\n",x2[i]);
			exit(1);
		}
    		num =  k*(*box).N[0] + j;
    		index[i] = num;
		cells[num].num_inside = cells[num].num_inside+1;

		if(x1[i]<cells[num].cell_center[0]){
			if(x2[i]<cells[num].cell_center[1]){
				cells[num].n_n[0] += (*gas).Fn/cells[num].volume;
				kk = k; jj = j-1;
				if(jj>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[1] += (*gas).Fn/cells[num2].volume;
				}
				kk = k-1; jj = j;
				if(kk>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[2] += (*gas).Fn/cells[num2].volume;
				}
				kk = k-1; jj = j-1;
				if(kk>=0 && jj>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[3] += (*gas).Fn/cells[num2].volume;
				}
			}
			else{
				cells[num].n_n[2] += (*gas).Fn/cells[num].volume;
				kk = k; jj = j-1;
				if(jj>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[3] += (*gas).Fn/cells[num2].volume;
				}
				kk = k+1; jj = j;
				if(kk<=(*box).N[1]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[0] += (*gas).Fn/cells[num2].volume;
				}
				kk = k+1; jj = j-1;
				if(kk<=(*box).N[1]-1 && jj>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[1] += (*gas).Fn/cells[num2].volume;
				}
			}
		}
		else{
			if(x2[i]<cells[num].cell_center[1]){
				cells[num].n_n[1] += (*gas).Fn/cells[num].volume;
				kk = k; jj = j+1;
				if(jj<=(*box).N[0]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[0] += (*gas).Fn/cells[num2].volume;
				}
				kk = k-1; jj = j;
				if(kk>=0){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[3] += (*gas).Fn/cells[num2].volume;
				}
				kk = k-1; jj = j+1;
				if(kk>=0 && jj<=(*box).N[0]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[2] += (*gas).Fn/cells[num2].volume;
				}
			}
			else{
				cells[num].n_n[3] += (*gas).Fn/cells[num].volume;
				kk = k; jj = j+1;
				if(jj<=(*box).N[0]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[2] += (*gas).Fn/cells[num2].volume;
				}
				kk = k+1; jj = j;
				if(kk<=(*box).N[1]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[1] += (*gas).Fn/cells[num2].volume;
				}
				kk = k+1; jj = j+1;
				if(kk<=(*box).N[1]-1 && jj<=(*box).N[0]-1){
					num2 =  kk*(*box).N[0] + jj;
					cells[num2].n_n[0] += (*gas).Fn/cells[num2].volume;
				}
			}
			}
		}
		for(jj=0; jj<(*box).N[0]; jj++){
			num2 =  0*(*box).N[0] + jj;
			cells[num2].n_n[0] = cells[num2].n_n[0]*2.0;
			cells[num2].n_n[1] = cells[num2].n_n[1]*2.0;
		}
		for(jj=0; jj<(*box).N[0]; jj++){
			num2 =  ((*box).N[1]-1)*(*box).N[0] + jj;
			cells[num2].n_n[2] = cells[num2].n_n[2]*2.0;
			cells[num2].n_n[3] = cells[num2].n_n[3]*2.0;
		}
		for(kk=0; kk<(*box).N[1]; kk++){
			num2 = kk*(*box).N[0] + 0;
			cells[num2].n_n[0] = cells[num2].n_n[0]*2.0;
			cells[num2].n_n[2] = cells[num2].n_n[2]*2.0;
		}
		for(kk=0; kk<(*box).N[1]; kk++){
			num2 = kk*(*box).N[0] + (*box).N[0]-1;
			cells[num2].n_n[1] = cells[num2].n_n[1]*2.0;
			cells[num2].n_n[3] = cells[num2].n_n[3]*2.0;
		}
  	}

   for (num=0; num<Ncells; num++){
	  free(cells[num].indices_inside);
	  cells[num].indices_inside = (int *) malloc( cells[num].num_inside * sizeof(int) );
    }

  for (i=0; i < (*gas).N; i++){
    num = index[i];
    cells[num].indices_inside[ numbering[num] ] = i;
    numbering[num]++;
  }
//if( (*gas).long_range == 2){
  double h,r_cut, mb;
	r_cut = 2.0*(*box).delta_dim[1];
	h = r_cut/4.0;
  int search, kkk;
  double r1, r2, rr, sqrt_pi;
  double r1r2;
  double W;
  sqrt_pi = sqrt(2.0*acos(-1.0));
  search = floor(r_cut/(*box).delta_dim[1]);
  mb = (*gas).Fn;// (*gas).n*(*box).Len[1]/4.0/(*gas).N;
  if((*box).problem != dcones && (*box).problem != flatnose){
    for(num=0; num< (*box).N[2]*(*box).N[1]*(*box).N[0]; num++){
	cells[num].n_smooth = 0.0;
	cells[num].dn2 = 0.0;
        cells[num].dn3 = 0.0;
    }
    double s2h2 = pow(1.0/h/h,2.0);
    for (j=0; j < (*box).N[0]; j++){
      for (k=0; k < (*box).N[1]; k++){
	num =  k*(*box).N[0] +j;
	r1 = cells[num].cell_center[1];
	for(kk = k-search; kk<= k+search; kk++){
		if(kk<0)
			kkk = 0;
		else if(kk>(*box).N[1]-1)
			kkk = (*box).N[1]-1;
		else
			kkk = kk;
		num2 = kkk*(*box).N[0] +j;
		for(i=0; i<cells[num2].num_inside; i++){
			id = cells[num2].indices_inside[i];
			r2 = x2[id]+(kk-kkk)*(*box).delta_dim[1];
			r1r2 = r2-r1;
			rr = fabs(r2-r1);
			if(rr<r_cut){
				cells[num].n_smooth += mb/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
				if(r2>r1)
					sgn = -1.0;
				else
					sgn = 1.0;
				//cells[num].dn2 += mb/sqrt_pi/h*exp(-rr*rr/h/h)*(-2.0*rr/h/h)*sgn;
				W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
				cells[num].dn2 += mb*W*r1r2*(-1.0/h/h);
				cells[num].dn3 += mb*W*s2h2*(3.0*r1r2-1.0/h/h*r1r2*r1r2*r1r2);
			}
		}
	}
      }
    }
  }
//}

	/*
	// Then find other properties
	for(i=0; i<(*gas).N; i++){
		r1 = x2[i];
		num = index[i];
		for(kk = num-search; kk<= num+search; kk++){
			if(kk<0)
				kkk = 0;
			else if(kk>(*box).N[1]-1)
				kkk = (*box).N[1]-1;
			else
				kkk = kk;
			for(id=0; id<cells[kkk].num_inside; id++){
				j = cells[kkk].indices_inside[id];
				r2 = x2[j]+(kk-kkk)*(*box).delta_dim[1];
				r1r2 = r2-r1;
				rr = fabs(r2-r1);
				if(rr<r_cut){
					W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
					Tnew[i] += mb*W*(*T)[j]/rho[j];
				}
			}
		}
	}
	double *temp;
	temp = *T;
	*T = Tnew;
	free(temp);
	*/
//}

  for(num=0; num< Ncells; num++){
    	cells[num].U_space[0] = 0.0;
	cells[num].U_space[1] = 0.0;
	cells[num].U_space[2] = 0.0;
	cells[num].weight = 0.0;
	for(j=0; j<6; j++)
	    cells[num].PIJ[j] = 0.0;
	for(j=0; j<10; j++)
	      cells[num].M3[j] = 0.0;
	for(j=0; j<6; j++)
	    cells[num].M4[j] = 0.0;
	for(j=0; j<3; j++)
	  cells[num].Q[j] = 0.0;
	for(j=0; j<3; j++)
	  cells[num].M5[j] = 0.0;
	cells[num].DM4 = 0.0;
	cells[num].DM2 = 0.0;
  for(j=0; j<3; j++)
    cells[num].Mp[j] = 0.0;
  }

  double w;
   //#pragma omp parallel for  private(i, id, num)
  for(num=0; num< Ncells; num++){
    for(i=0; i<cells[num].num_inside; i++){
      id = cells[num].indices_inside[i];
      w = 1.0;
      cells[num].U_space[0] = cells[num].U_space[0]+U1[id]*w;
      cells[num].U_space[1] = cells[num].U_space[1]+U2[id]*w;
      cells[num].U_space[2] = cells[num].U_space[2]+U3[id]*w;
      cells[num].weight += w;
    }
  }
  for(num=0; num<Ncells; num++){
    if(cells[num].num_inside >0){
	    	cells[num].U_space[0] = cells[num].U_space[0]/cells[num].weight;
    		cells[num].U_space[1] = cells[num].U_space[1]/cells[num].weight;
    		cells[num].U_space[2] = cells[num].U_space[2]/cells[num].weight;
    }
  }
  //#pragma omp parallel for  private(i, id, energy)
  for(num=0; num< Ncells; num++){
    for(i=0; i<cells[num].num_inside; i++){
      id = cells[num].indices_inside[i];
      w = 1.0;
      Mp1 = U1[id]-cells[num].U_space[0];
      Mp2 = U2[id]-cells[num].U_space[1];
      Mp3 = U3[id]-cells[num].U_space[2];

      cells[num].Mp[0] += Mp1*w;
      cells[num].Mp[1] += Mp2*w;
      cells[num].Mp[2] += Mp3*w;

      cells[num].PIJ[0] = cells[num].PIJ[0] + Mp1 * Mp1*w;//11
      cells[num].PIJ[1] = cells[num].PIJ[1] + Mp1 * Mp2*w;//12
      cells[num].PIJ[2] = cells[num].PIJ[2] + Mp1 * Mp3*w;//13

      cells[num].PIJ[3] = cells[num].PIJ[3] + Mp2 * Mp2*w;//22
      cells[num].PIJ[4] = cells[num].PIJ[4] + Mp2 * Mp3*w;//23

      cells[num].PIJ[5] = cells[num].PIJ[5] + Mp3 * Mp3*w;//33

    // < M^p_k * M^p_k * M^p_i > called Q
      energy = Mp1*Mp1 + Mp2*Mp2 + Mp3*Mp3;

    // < M^p_i * M^p_j * M^p_k > called M3
      cells[num].M3[0] = cells[num].M3[0] + Mp1 * Mp1 * Mp1*w;//111
      cells[num].M3[1] = cells[num].M3[1] + Mp1 * Mp1 * Mp2*w;//112
      cells[num].M3[2] = cells[num].M3[2] + Mp1 * Mp1 * Mp3*w;//113

      cells[num].M3[3] = cells[num].M3[3] + Mp1 * Mp2 * Mp2*w;//122
      cells[num].M3[4] = cells[num].M3[4] + Mp1 * Mp2 * Mp3*w;//123

      cells[num].M3[5] = cells[num].M3[5] + Mp1 * Mp3 * Mp3*w;//133

      cells[num].M3[6] = cells[num].M3[6] + Mp2 * Mp2 * Mp2*w;//222
      cells[num].M3[7] = cells[num].M3[7] + Mp2 * Mp2 * Mp3*w;//223
      cells[num].M3[8] = cells[num].M3[8] + Mp2 * Mp3 * Mp3*w;//233

      cells[num].M3[9] = cells[num].M3[9] + Mp3 * Mp3 * Mp3*w;//333

    // < M^p_k * M^p_k * M^p_i * M^p_j> called M4
      cells[num].M4[0] = cells[num].M4[0] + Mp1 * Mp1 * energy*w;//11
      cells[num].M4[1] = cells[num].M4[1] + Mp1 * Mp2 * energy*w;//12
      cells[num].M4[2] = cells[num].M4[2] + Mp1 * Mp3 * energy*w;//13

      cells[num].M4[3] = cells[num].M4[3] + Mp2 * Mp2 * energy*w;//22
      cells[num].M4[4] = cells[num].M4[4] + Mp2 * Mp3 * energy*w;//23

      cells[num].M4[5] = cells[num].M4[5] + Mp3 * Mp3 * energy*w;//33

//    < M^p_i * M^p_j * M^p_j *  M^p_k * M^p_k>    called M5
      cells[num].M5[0] = cells[num].M5[0] + Mp1 * energy * energy*w;//1
      cells[num].M5[1] = cells[num].M5[1] + Mp2 * energy * energy*w;//2
      cells[num].M5[2] = cells[num].M5[2] + Mp3 * energy * energy*w;//3

    // < (M^p_k * M^p_k)^2 > called DM2
      cells[num].DM2 = cells[num].DM2 + energy*w;

    }
  }
  double Y;
  for(num=0; num<Ncells; num++){
     if(cells[num].num_inside > 0){
    	for(j=0; j<6; j++)
	    cells[num].PIJ[j] = cells[num].PIJ[j]/cells[num].weight;
	for(j=0; j<10; j++)
	      cells[num].M3[j] = cells[num].M3[j]/cells[num].weight;
	for(j=0; j<6; j++)
	    cells[num].M4[j] = cells[num].M4[j]/cells[num].weight;
	for(j=0; j<3; j++)
	    cells[num].M5[j] = cells[num].M5[j]/cells[num].weight;
	for(j=0; j<3; j++)
	  cells[num].Q[j] = cells[num].Q[j]/cells[num].weight;
	for(j=0; j<3; j++)
          cells[num].Mp[j] = cells[num].Mp[j]/cells[num].weight;
        cells[num].DM2 = cells[num].DM2/cells[num].weight;
    	cells[num].DM4 = cells[num].M4[0]+ cells[num].M4[3]+cells[num].M4[5];
    	cells[num].Q[0] = cells[num].M3[0]+cells[num].M3[3]+cells[num].M3[5];
    	cells[num].Q[1] = cells[num].M3[1]+cells[num].M3[6]+cells[num].M3[8];
    	cells[num].Q[2] = cells[num].M3[2]+cells[num].M3[7]+cells[num].M3[9];
    }



    cells[num].T = (*gas).m*cells[num].DM2/(3.0*(*gas).kb);

    cells[num].n = cells[num].n_ratio*(*gas).Fn*cells[num].weight/cells[num].volume;
    Y = increase_collision_rate(cells[num].n, gas);
    if(cells[num].weight > 1e-14){
	    cells[num].gamma_st =  -0.001*( cells[num].n*(*gas).b*Y/( (*gas).kb*cells[num].T/(*gas).m ) );
    }
  }
//	 (*gas).n = 100;
	//for(num=0; num<Ncells; num++){
	//	cells[num].n = (*gas).n;
	//	for(i=0; i<4; i++)
	//		cells[num].n_n[i] = (*gas).n;
	//}

	if( (*gas).model == "SPH" || (*gas).model =="SPH_DSMC" ){
		//(*gas).N_SPH = 0;
		for(i=0; i<(*gas).N; i++){
			num = index[i];
			if( cells[num].model == 1 ){
				color[i] = 1;
				//(*gas).N_SPH ++;
			}
			else{
				color[i] = 4;
			}
		}

		// mean temperature for SPH cells only
	  for(num=0; num<Ncells; num++){
		 if( cells[num].model == 1 )// if it was SPH
				cells[num].T = 0.0;
		}
		for(i=0; i<(*gas).N; i++){
			num = index[i];
			if( cells[num].model == 1 )
				cells[num].T += (*T)[i];
		}
		for(num=0; num<Ncells; num++){
			cells[num].T = cells[num].T/cells[num].weight;
		}
}


if( (*gas).model == "SPH" || (*gas).model =="SPH_DSMC" ){
	int search, kk, kkk;
	double W, sqrt_pi, mb, rr, r1r2, r2, r1, h, r_cut, s2h2;

	double cv = 3.0*(*gas).kb/(*gas).m/2.0;
  double gama = 1.4;

	for(i=0; i<(*gas).N; i++){
		rho[i] = 0.0;
	}
	r_cut = (*box).delta_dim[1]/2.0;
	h = r_cut/10.0;
	(*gas).r_cut = r_cut; (*gas).h = h;
	search = floor(r_cut/(*box).delta_dim[1])+1;
	mb = (*gas).m*(*gas).Fn;//*(*box).Len[1]/(*gas).N;
	(*gas).mb = mb;
	sqrt_pi = sqrt(2.0*acos(-1.0));

	for(num=0; num<Ncells; num++){
	 if( cells[num].model == 1 ){// if it was SPH
			cells[num].PIJ[3] = 0.0;
			cells[num].Q[1] = 0.0;
		}
	}
	// first find rhoi
	double mu,k,DW,dudx,dTdx;
	for(i=0; i<(*gas).N; i++){
		r1 = x2[i];
		num = index[i];
		mu = (*gas).mu*pow((*T)[i]/(*gas).T0,0.5);
    k = (*gas).k/(*gas).mu*mu;
		for(kk = num-search; kk<= num+search; kk++){
			if(kk<0)
				kkk = 0;
			else if(kk>(*box).N[1]-1)
				kkk = (*box).N[1]-1;
			else
				kkk = kk;
			for(id=0; id<cells[kkk].num_inside; id++){
				j = cells[kkk].indices_inside[id];
				r2 = x2[j]+(kk-kkk)*(*box).delta_dim[1];
				r1r2 = r2-r1;
				rr = fabs(r2-r1);
				if(rr<r_cut){
					W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
					rho[i] += mb*W;
				}
			}
		}
	}
	for(i=0; i<(*gas).N; i++){
		r1 = x2[i];
		num = index[i];
		if(cells[num].model == 1){
			mu = (*gas).mu*pow((*T)[i]/(*gas).T0,0.5);
    	k = (*gas).k/(*gas).mu*mu;
			for(kk = num-search; kk<= num+search; kk++){
				if(kk<0)
					kkk = 0;
				else if(kk>(*box).N[1]-1)
					kkk = (*box).N[1]-1;
				else
					kkk = kk;
				for(id=0; id<cells[kkk].num_inside; id++){
					j = cells[kkk].indices_inside[id];
					r2 = x2[j]+(kk-kkk)*(*box).delta_dim[1];
					r1r2 = r2-r1;
					rr = fabs(r2-r1);
					if(rr<r_cut){

						W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
						if(i!=j){
            	DW = W*r1r2*(-1.0/h/h);
          	}
          	else{
            	DW = 0.0;
          	}
						dudx   += mb*U2[j]/rho[j]*DW;
          	dTdx   += mb*(*T)[j]/rho[j]*DW;
						cells[num].PIJ[3] += 4.0/3.0*mu*dudx;
          	cells[num].Q[1] += -2.0*k*dTdx;
					}
				}
			}
		}
	}
}

  double tau, mu;
for(num=0; num<Ncells; num++){
    if(cells[num].weight > 1e-14){
      scale_mu(cells[num].T, (*gas).T0, &mu, (*gas).mu, gas);
      tau = 2.0*mu/( cells[num].n*(*gas).kb*cells[num].T );
      cells[num].Lambda = (cells[num].PIJ[0]-cells[num].DM2/3.0)*( (cells[num].PIJ[3]-cells[num].DM2/3.0)*(cells[num].PIJ[5]-cells[num].DM2/3.0)-cells[num].PIJ[4]*cells[num].PIJ[4] );
      cells[num].Lambda = cells[num].Lambda - cells[num].PIJ[1]*( cells[num].PIJ[1]*(cells[num].PIJ[5]-cells[num].DM2/3.0)-cells[num].PIJ[4]*cells[num].PIJ[2] );
      cells[num].Lambda = cells[num].Lambda + cells[num].PIJ[2]*( cells[num].PIJ[1]*cells[num].PIJ[4]-cells[num].PIJ[2]*(cells[num].PIJ[3]-cells[num].DM2/3.0) );
      cells[num].Lambda = - cells[num].Lambda/( tau*pow(cells[num].DM2,4.0) );
   }
}
if((*gas).model ==  "FP_dense" || (*gas).model == "FP_ideal_gas" || (*gas).model == "FP_linear"){
     if((*box).step>0 && (*gas).model ==  "FP_dense")
        alpha_calculation( U1, U2, U3, gas, cells, box);
     c_gamma_calculation( U1, U2, U3, gas, cells, box);
  }
  for(num=0; num<Ncells; num++){
      cells[num].phi[0] =  cells[num].alpha[0]*cells[num].Mp[0]
                          +cells[num].alpha[1]*cells[num].Mp[1]
                          +cells[num].alpha[2]*cells[num].Mp[2];

      cells[num].phi[1] =  cells[num].alpha[3]*cells[num].Mp[0]
                          +cells[num].alpha[4]*cells[num].Mp[1]
                          +cells[num].alpha[5]*cells[num].Mp[2];

      cells[num].phi[2] =  cells[num].alpha[6]*cells[num].Mp[0]
                          +cells[num].alpha[7]*cells[num].Mp[1]
                          +cells[num].alpha[8]*cells[num].Mp[2];
  }
 free(dx);
 free(numbering);
}

double sign(double x){
      if(x>1e-16)
            return 1.0;
      else
            return -1.0;
}
void intersect_3D_x(double x, double *y, double *z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new){
      // input is x
      // output is corresponding y and z
      double n[3];
      int i;
      double size = 0.0;
      n[0] = x_new-x_old;
      n[1] = y_new-y_old;
      n[2] = z_new-z_old;
      for(i=0; i<3; i++)
            size += n[i]*n[i];
      size = sqrt(size);
      for(i=0; i<3; i++)
            n[i] = n[i]/size;
      // until here the code is the same for all directions
      if(fabs(n[0])>1e-16){
            *y = (n[1]/n[0])*(x-x_old)+y_old;
            *z = (n[2]/n[0])*(x-x_old)+z_old;
      }
      else{
            *y = y_new;
            *z = z_new;
      }
}
void intersect_3D_y(double *x, double y, double *z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new){
      // input is x
      // output is corresponding y and z
      double n[3];
      int i;
      double size = 0.0;
      n[0] = x_new-x_old;
      n[1] = y_new-y_old;
      n[2] = z_new-z_old;
      for(i=0; i<3; i++)
            size += n[i]*n[i];
      size = sqrt(size);
      for(i=0; i<3; i++)
            n[i] = n[i]/size;
      // until here the code is the same for all directions
      if(fabs(n[1])>1e-16){
            *x = (n[0]/n[1])*(y-y_old)+x_old;
            *z = (n[2]/n[1])*(y-y_old)+z_old;
      }
      else{
            *x = x_new;
            *z = z_new;
      }
}
void intersect_3D_z(double *x, double *y, double z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new){
      // input is x
      // output is corresponding y and z
      double n[3];
      int i;
      double size = 0.0;
      n[0] = x_new-x_old;
      n[1] = y_new-y_old;
      n[2] = z_new-z_old;
      for(i=0; i<3; i++)
            size += n[i]*n[i];
      size = sqrt(size);
      for(i=0; i<3; i++)
            n[i] = n[i]/size;
      // until here the code is the same for all directions
      if(fabs(n[2])>1e-16){
            *x = (n[0]/n[2])*(z-z_old)+x_old;
            *y = (n[1]/n[2])*(z-z_old)+y_old;
      }
      else{
            *x = x_new;
            *y = y_new;
      }
}
double line_interesect_y(double x, double x_old, double y_old, double x_new, double y_new){
      return (x-x_old)*(y_new-y_old)/(x_new-x_old)+y_old;
}
double line_interesect_x(double y, double x_old, double y_old, double x_new, double y_new){
      return (y-y_old)*(x_new-x_old)/(y_new-y_old)+x_old;
}
void update_face_info(double *U1,double *U2,double *U3, double *x1,double *x2, double *x1_old,double *x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index) {
      int id, i_old, i_new, j_old, j_new, i, j;
      double x,y, sgn;
      int num;
      double dt, energy;
      double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
      dt = 0.0;
      for(num=0; num<(*box).N[0]*(*box).N[1]*(*box).N[2]; num++){
          for(j=0; j<6; j++){
            cells[num].M_n[j] = 0.0;
            cells[num].M1_f[j] = 0.0;
            cells[num].M2_f[j] = 0.0;
            cells[num].M3_f[j] = 0.0;
            cells[num].T_f[j] = 0.0;
          }
					//for(j=0; j<4; j++){
					//	cells[num].sum_N[j] = 0.0;
					//}
          if( (*box).direction[0] == 0 ){
                cells[num].M_n[0] = 1.0;
                cells[num].M_n[1] = 1.0;
                cells[num].sum_M_n[0] = 1.0;
                cells[num].sum_M_n[1] = 1.0;
								cells[num].sum_M_np[0] = 0.0;
                cells[num].sum_M_np[1] = 0.0;
          }
          if( (*box).direction[1] == 0 ){
               cells[num].M_n[2] = 1.0;
               cells[num].M_n[3] = 1.0;
               cells[num].sum_M_n[2] = 1.0;
               cells[num].sum_M_n[3] = 1.0;
							 cells[num].sum_M_np[2] = 0.0;
               cells[num].sum_M_np[3] = 0.0;
          }
          if( (*box).direction[2] == 0 ){
               cells[num].M_n[4] = 1.0;
               cells[num].M_n[5] = 1.0;
               cells[num].sum_M_n[4] = 1.0;
               cells[num].sum_M_n[5] = 1.0;
							 cells[num].sum_M_np[4] = 0.0;
               cells[num].sum_M_np[5] = 0.0;
          }
      }

      for(id=0; id<(*gas).N; id++){
            // intersect with vertical lines
						if(	(*box).direction[0]==1 ){
            	i_old = floor(x1_old[id]/(*box).delta_dim[0]);
            	i_new = floor(x1[id]    /(*box).delta_dim[0]);
						}
						else{
							i_old = 0;	i_new =0;
						}
            if(i_old != i_new){
                  sgn = sign(i_new-i_old);
                  x = (i_old+sgn/2.0+1.0/2.0)*(*box).delta_dim[0];
                  y = line_interesect_y(x, x1_old[id], x2_old[id], x1[id], x2[id]);
                  while(sgn*(x-x1_old[id]) > 0.0 && sgn*(x-x1[id]) < 0.0 && x>-epsilon && x<(*box).Len[0]+epsilon && y>-epsilon && y<(*box).Len[1]+epsilon){

                        j = floor(y/(*box).delta_dim[1]);

                        i = floor((x+epsilon)/(*box).delta_dim[0]);
                        if(i<(*box).N[0] && j>-1 && j<(*box).N[1]){
                              num = j*(*box).N[0] + i;
                              energy =    (U1[id]-  cells[num].U_space[0])*(U1[id]-  cells[num].U_space[0])
                        		+ (U2[id]-  cells[num].U_space[1])*(U2[id]-  cells[num].U_space[1])
                        		+ (U3[id]-  cells[num].U_space[2])*(U3[id]-  cells[num].U_space[2]);
			      dt = 0.0;
                              dt = fabs(U1[id]);
			   if(fabs(dt)>1e-14){
                              cells[num].M_n[0] += 1.0/dt;

                              cells[num]    .M1_f[ 0 ] +=  U1[id]/dt;
                              cells[num]    .M2_f[ 0 ] +=  U2[id]/dt;
                              cells[num]    .M3_f[ 0 ] +=  U3[id]/dt;
                              cells[num]    .T_f [ 0 ] +=  energy/dt;

                              if( (*box).step > (*gas).times*(*box).after){

                                    cells[num]    .sum_M_n[0]    += 1.0/dt;

																		cells[num]		.sum_T_f[0]		 += energy/dt;

                                    cells[num]    .sum_M1_f[0]   += U1[id]/dt;
                                    cells[num]    .sum_M2_f[0]   += U2[id]/dt;
                                    cells[num]    .sum_M3_f[0]   += U3[id]/dt;

                                    cells[num]    .sum_MM_f_0[0] += U1[id]*U1[id]/dt;
                                    cells[num]    .sum_MM_f_0[1] += U1[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_0[2] += U1[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_0[3] += U2[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_0[4] += U2[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_0[5] += U3[id]*U3[id]/dt;
                              }
			   }
                        }

                        i = floor((x-epsilon)/(*box).delta_dim[0]);
                        if(i>-1 && j>-1 && j<(*box).N[1]){
                              num = j*(*box).N[0] + i;
                              energy =    (U1[id]-  cells[num].U_space[0])*(U1[id]-  cells[num].U_space[0])
                                    + (U2[id]-  cells[num].U_space[1])*(U2[id]-  cells[num].U_space[1])
                                    + (U3[id]-  cells[num].U_space[2])*(U3[id]-  cells[num].U_space[2]);
			      dt = 0.0;
                              dt = fabs(U1[id]);
			   if(fabs(dt)>1e-14){
                              cells[num].M_n[1] += 1.0/dt;

                              cells[num]    .M1_f[ 1 ] +=  U1[id]/dt;
                              cells[num]    .M2_f[ 1 ] +=  U2[id]/dt;
                              cells[num]    .M3_f[ 1 ] +=  U3[id]/dt;
                              cells[num]    .T_f [ 1 ] +=  energy/dt;

                              if( (*box).step > (*gas).times*(*box).after){

                                    cells[num]    .sum_M_n[1]    += 1.0/dt;

																		cells[num]		.sum_T_f[1]		 += energy/dt;

                                    cells[num]    .sum_M1_f[1]   += U1[id]/dt;
                                    cells[num]    .sum_M2_f[1]   += U2[id]/dt;
                                    cells[num]    .sum_M3_f[1]   += U3[id]/dt;

                                    cells[num]    .sum_MM_f_1[0] += U1[id]*U1[id]/dt;
                                    cells[num]    .sum_MM_f_1[1] += U1[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_1[2] += U1[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_1[3] += U2[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_1[4] += U2[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_1[5] += U3[id]*U3[id]/dt;
                              }
			   }
                        }

                        x = x +sgn*(*box).delta_dim[0];
                        y = line_interesect_y(x, x1_old[id], x2_old[id], x1[id], x2[id]);
                  }
            }
            // intersect with horizontal lines
            // intersect with vertical lines
            j_old = floor(x2_old[id]/(*box).delta_dim[1]);
            j_new = floor(x2[id]    /(*box).delta_dim[1]);
            if(j_old != j_new){
                  sgn = sign(j_new-j_old);
                  y = (j_old+sgn/2.0+1.0/2.0)*(*box).delta_dim[1];
                  x = line_interesect_x(y, x1_old[id], x2_old[id], x1[id], x2[id]);
                  while(sgn*(y-x2_old[id]) > 0.0 && sgn*(y-x2[id]) < 0.0 && x>-epsilon && x<(*box).Len[0]+epsilon && y>-epsilon && y<(*box).Len[1]+epsilon){

                        i = floor(x/(*box).delta_dim[0]);

                        j = floor((y+epsilon)/(*box).delta_dim[1]);
                        if(j<(*box).N[1] && i>-1 && i<(*box).N[0]){
                              num = j*(*box).N[0] + i;
                              energy =    (U1[id]-  cells[num].U_space[0])*(U1[id]-  cells[num].U_space[0])
                                    + (U2[id]-  cells[num].U_space[1])*(U2[id]-  cells[num].U_space[1])
                                    + (U3[id]-  cells[num].U_space[2])*(U3[id]-  cells[num].U_space[2]);
			      dt = 0.0;
                              //dt = fabs(U2[id]);
															dt = fabs( (x2[id]-x2_old[id])/(*gas).delta_t );
			   if(fabs(dt)>1e-14){
                              cells[num].M_n[2] += 1.0/dt;
                              cells[num]    .M1_f[ 2 ] +=  U1[id]/dt;
                              cells[num]    .M2_f[ 2]  +=  U2[id]/dt;
                              cells[num]    .M3_f[ 2 ] +=  U3[id]/dt;
                              cells[num]    .T_f [ 2 ] +=  energy/dt;

                              if( (*box).step > (*gas).times*(*box).after){

																		cells[num]    .sum_N[2]    += (x2[id]-x2_old[id])/(*gas).delta_t/dt;

                                    cells[num]    .sum_M_n[2]    += 1.0/dt;
																		if( (x2[id]-x2_old[id])/(*gas).delta_t > 0.0){
																			cells[num].sum_M_np[2] += 1.0/dt;
																		}
																		cells[num]		.sum_T_f[2]		 += energy/dt;

                                    cells[num]    .sum_M1_f[2]   += U1[id]/dt;
                                    cells[num]    .sum_M2_f[2]   += U2[id]/dt;
                                    cells[num]    .sum_M3_f[2]   += U3[id]/dt;

                                    cells[num]    .sum_MM_f_2[0] += U1[id]*U1[id]/dt;
                                    cells[num]    .sum_MM_f_2[1] += U1[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_2[2] += U1[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_2[3] += U2[id]*U2[id]/dt;
                                    cells[num]    .sum_MM_f_2[4] += U2[id]*U3[id]/dt;
                                    cells[num]    .sum_MM_f_2[5] += U3[id]*U3[id]/dt;
                              }
			   }
                        }

                        j = floor((y-epsilon)/(*box).delta_dim[1]);
                        if(j>-1 && i>-1 && i<(*box).N[0]){
                              num = j*(*box).N[0] + i;
                              energy =    (U1[id]-  cells[num].U_space[0])*(U1[id]-  cells[num].U_space[0])
                                    + (U2[id]-  cells[num].U_space[1])*(U2[id]-  cells[num].U_space[1])
                                    + (U3[id]-  cells[num].U_space[2])*(U3[id]-  cells[num].U_space[2]);
			      dt = 0.0;
                              //dt = fabs(U2[id]);
															dt = fabs( (x2[id]-x2_old[id])/(*gas).delta_t );
			      if(fabs(dt)>1e-14){
	                              cells[num].M_n[3] += 1.0/dt;
        	                      cells[num]    .M1_f[ 3 ] +=  U1[id]/dt;
        	                      cells[num]    .M2_f[ 3 ] +=  U2[id]/dt;
        	                      cells[num]    .M3_f[ 3 ] +=  U3[id]/dt;
        	                      cells[num]    .T_f [ 3 ] +=  energy/dt;

        	                      if( (*box).step > (*gas).times*(*box).after){
																			cells[num]    .sum_N[3]    += (x2[id]-x2_old[id])/(*gas).delta_t/dt;

																			if( (x2[id]-x2_old[id])/(*gas).delta_t > 0.0){
																				cells[num].sum_M_np[3] += 1.0/dt;
																			}

        	                            cells[num]    .sum_M_n[3]    += 1.0/dt;

																			cells[num]		.sum_T_f[3]		 += energy/dt;

        	                            cells[num]    .sum_M1_f[3]   += U1[id]/dt;
        	                            cells[num]    .sum_M2_f[3]   += U2[id]/dt;
        	                            cells[num]    .sum_M3_f[3]   += U3[id]/dt;

        	                            cells[num]    .sum_MM_f_3[0] += U1[id]*U1[id]/dt;
        	                            cells[num]    .sum_MM_f_3[1] += U1[id]*U2[id]/dt;
        	                            cells[num]    .sum_MM_f_3[2] += U1[id]*U3[id]/dt;
        	                            cells[num]    .sum_MM_f_3[3] += U2[id]*U2[id]/dt;
        	                            cells[num]    .sum_MM_f_3[4] += U2[id]*U3[id]/dt;
        	                            cells[num]    .sum_MM_f_3[5] += U3[id]*U3[id]/dt;
        	                      }
			      }
                        }

                        y = y +sgn*(*box).delta_dim[1];
                        x = line_interesect_x(y, x1_old[id], x2_old[id], x1[id], x2[id]);
                  }
            }
      }
}

void update_face_info_BC(double U1,double U2,double U3, double x1,double x2, double x1_old,double x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells) {
      int i_old, i_new, j_old, j_new, i, j;
      double x,y, sgn;
      int num;
      double energy, dt;
      double epsilon = max((*box).delta_dim[0],(*box).delta_dim[1])/100;

      i_old = floor(x1_old/(*box).delta_dim[0]);
      i_new = floor(x1    /(*box).delta_dim[0]);
      if(i_old != i_new){
            sgn = sign(i_new-i_old);
            x = (i_old+sgn/2.0+1.0/2.0)*(*box).delta_dim[0];
            y = line_interesect_y(x, x1_old, x2_old, x1, x2);
            while(sgn*(x-x1_old) > 0.0 && sgn*(x-x1) < 0.0 && x>-epsilon && x<(*box).Len[0]+epsilon && y>-epsilon && y<(*box).Len[1]+epsilon){

                  j = floor(y/(*box).delta_dim[1]);

                  i = floor((x+epsilon)/(*box).delta_dim[0]);
                  if(i<(*box).N[0] && j>-1 && j<(*box).N[1]){
                        num = j*(*box).N[0] + i;
                        energy =    (U1-  cells[num].U_space[0])*(U1-  cells[num].U_space[0])
                              + (U2-  cells[num].U_space[1])*(U2-  cells[num].U_space[1])
                              + (U3-  cells[num].U_space[2])*(U3-  cells[num].U_space[2]);

                        dt = fabs(U1);
                        cells[num].M_n[0] += 1.0/dt;

                        cells[num]    .M1_f[ 0 ] +=  U1/dt;
                        cells[num]    .M2_f[ 0 ] +=  U2/dt;
                        cells[num]    .M3_f[ 0 ] +=  U3/dt;
                        cells[num]    .T_f [ 0 ] +=  energy/dt;

                        if( (*box).step > (*gas).times*(*box).after){
                              cells[num]    .sum_M_n[0]    += 1.0/dt;

                              cells[num]    .sum_M1_f[0]   += U1/dt;
                              cells[num]    .sum_M2_f[0]   += U2/dt;
                              cells[num]    .sum_M3_f[0]   += U3/dt;

                              cells[num]    .sum_MM_f_0[0] += U1*U1/dt;
                              cells[num]    .sum_MM_f_0[1] += U1*U2/dt;
                              cells[num]    .sum_MM_f_0[2] += U1*U3/dt;
                              cells[num]    .sum_MM_f_0[3] += U2*U2/dt;
                              cells[num]    .sum_MM_f_0[4] += U2*U3/dt;
                              cells[num]    .sum_MM_f_0[5] += U3*U3/dt;
                        }
                  }

                  i = floor((x-epsilon)/(*box).delta_dim[0]);
                  if(i>-1 && j>-1 && j<(*box).N[1]){
                        num = j*(*box).N[0] + i;
                        energy =    (U1-  cells[num].U_space[0])*(U1-  cells[num].U_space[0])
                              + (U2-  cells[num].U_space[1])*(U2-  cells[num].U_space[1])
                              + (U3-  cells[num].U_space[2])*(U3-  cells[num].U_space[2]);

                        dt = fabs(U1);
                        cells[num].M_n[1] += 1.0/dt;

                        cells[num]    .M1_f[ 1 ] +=  U1/dt;
                        cells[num]    .M2_f[ 1 ] +=  U2/dt;
                        cells[num]    .M3_f[ 1 ] +=  U3/dt;
                        cells[num]    .T_f [ 1 ] +=  energy/dt;

                        if( (*box).step > (*gas).times*(*box).after){
                              cells[num]    .sum_M_n[1]    += 1.0/dt;

                              cells[num]    .sum_M1_f[1]   += U1/dt;
                              cells[num]    .sum_M2_f[1]   += U2/dt;
                              cells[num]    .sum_M3_f[1]   += U3/dt;

                              cells[num]    .sum_MM_f_1[0] += U1*U1/dt;
                              cells[num]    .sum_MM_f_1[1] += U1*U2/dt;
                              cells[num]    .sum_MM_f_1[2] += U1*U3/dt;
                              cells[num]    .sum_MM_f_1[3] += U2*U2/dt;
                              cells[num]    .sum_MM_f_1[4] += U2*U3/dt;
                              cells[num]    .sum_MM_f_1[5] += U3*U3/dt;
                        }
                  }

                  x = x +sgn*(*box).delta_dim[0];
                  y = line_interesect_y(x, x1_old, x2_old, x1, x2);
            }
      }
      // intersect with horizontal lines
      // intersect with vertical lines
      j_old = floor(x2_old/(*box).delta_dim[1]);
      j_new = floor(x2    /(*box).delta_dim[1]);
      if(j_old != j_new){
            sgn = sign(j_new-j_old);
            y = (j_old+sgn/2.0+1.0/2.0)*(*box).delta_dim[1];
            x = line_interesect_x(y, x1_old, x2_old, x1, x2);
            while(sgn*(y-x2_old) > 0.0 && sgn*(y-x2) < 0.0 && x>-epsilon && x<(*box).Len[0]+epsilon && y>-epsilon && y<(*box).Len[1]+epsilon){

                  i = floor(x/(*box).delta_dim[0]);

                  j = floor((y+epsilon)/(*box).delta_dim[1]);
                  if(j<(*box).N[1] && i>-1 && i<(*box).N[0]){
                        num = j*(*box).N[0] + i;
                        energy =    (U1-  cells[num].U_space[0])*(U1-  cells[num].U_space[0])
                              + (U2-  cells[num].U_space[1])*(U2-  cells[num].U_space[1])
                              + (U3-  cells[num].U_space[2])*(U3-  cells[num].U_space[2]);

                        dt = fabs(U2);
                        cells[num].M_n[2] += 1.0/dt;
                        cells[num]    .M1_f[ 2 ] +=  U1/dt;
                        cells[num]    .M2_f[ 2 ] +=  U2/dt;
                        cells[num]    .M3_f[ 2 ] +=  U3/dt;
                        cells[num]    .T_f [ 2 ] +=  energy/dt;

                        if( (*box).step > (*gas).times*(*box).after){
                              cells[num]    .sum_M_n[2]    += 1.0/dt;

                              cells[num]    .sum_M1_f[2]   += U1/dt;
                              cells[num]    .sum_M2_f[2]   += U2/dt;
                              cells[num]    .sum_M3_f[2]   += U3/dt;

                              cells[num]    .sum_MM_f_2[0] += U1*U1/dt;
                              cells[num]    .sum_MM_f_2[1] += U1*U2/dt;
                              cells[num]    .sum_MM_f_2[2] += U1*U3/dt;
                              cells[num]    .sum_MM_f_2[3] += U2*U2/dt;
                              cells[num]    .sum_MM_f_2[4] += U2*U3/dt;
                              cells[num]    .sum_MM_f_2[5] += U3*U3/dt;
                        }
                  }

                  j = floor((y-epsilon)/(*box).delta_dim[1]);
                  if(j>-1 && i>-1 && i<(*box).N[0]){
                        num = j*(*box).N[0] + i;

                        energy =    (U1-  cells[num].U_space[0])*(U1-  cells[num].U_space[0])
                              + (U2-  cells[num].U_space[1])*(U2-  cells[num].U_space[1])
                              + (U3-  cells[num].U_space[2])*(U3-  cells[num].U_space[2]);

                        dt = fabs(U2);
                        cells[num].M_n[3] += 1.0/dt;
                        cells[num]    .M1_f[ 3 ] +=  U1/dt;
                        cells[num]    .M2_f[ 3 ] +=  U2/dt;
                        cells[num]    .M3_f[ 3 ] +=  U3/dt;
                        cells[num]    .T_f [ 3 ] +=  energy/dt;

                        if( (*box).step > (*gas).times*(*box).after){
                              cells[num]    .sum_M_n[3]    += 1.0/dt;

                              cells[num]    .sum_M1_f[3]   += U1/dt;
                              cells[num]    .sum_M2_f[3]   += U2/dt;
                              cells[num]    .sum_M3_f[3]   += U3/dt;

                              cells[num]    .sum_MM_f_3[0] += U1*U1/dt;
                              cells[num]    .sum_MM_f_3[1] += U1*U2/dt;
                              cells[num]    .sum_MM_f_3[2] += U1*U3/dt;
                              cells[num]    .sum_MM_f_3[3] += U2*U2/dt;
                              cells[num]    .sum_MM_f_3[4] += U2*U3/dt;
                              cells[num]    .sum_MM_f_3[5] += U3*U3/dt;
                        }
                  }

                  y = y +sgn*(*box).delta_dim[1];
                  x = line_interesect_x(y, x1_old, x2_old, x1, x2);
            }
      }
}

void shuffle(double *x1,double *x2,double *x3, struct GAS *gas){
  // inputs: 	positions x1, x2, x3
  //		gas data
  //		dphi_dx_p
  double d = (*gas).sigma;
  // Now calculate the dphi for each particle
  double r[3],R;
//   double pi = acos(-1.0);
  int done = 1;
  while(done){
  done = 0;
  int number = 0;
  for(int ii=0; ii<(*gas).N; ii++){
    for(int jj=0; jj<(*gas).N; jj++){
      if(ii != jj){
	r[0] = x1[jj] - x1[ii];
	r[1] = x2[jj] - x2[ii];
	r[2] = x3[jj] - x3[ii];
	R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	if ( R - d < -0.01*d ){
	  done = 1;
	  x1[ii] = x1[ii] - (d/2.0)*(r[0]/R);
	  x2[ii] = x2[ii] - (d/2.0)*(r[1]/R);
	  x3[ii] = x3[ii] - (d/2.0)*(r[2]/R);

	  x1[jj] = x1[jj] + (d/2.0)*(r[0]/R);
	  x2[jj] = x2[jj] + (d/2.0)*(r[1]/R);
	  x3[jj] = x3[jj] + (d/2.0)*(r[2]/R);
	  number ++;
	}
      }
    }
  }
  printf("number = %d\n",number);
}
}

void min_rij(int n, double *rx, double *ry, double *rz, int *nn, double *rr){
      int i;
      double r;
      *nn = 0;
      *rr = 1e15;
      for(i=0; i<n; i++){
            r = sqrt(rx[i]*rx[i]+ry[i]*ry[i]+rz[i]*rz[i]);
            if(r < *rr){
                  *rr = r;
                  *nn = i;
            }
      }
}

void poisson_solver(double *U2, double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
	//Thompson_3d(double *a,double *b,double *c,double *x, double *RHS,int n)
	int nf = (*box).N[1]+1;
	double *a = (double *) malloc( nf * sizeof(double) );
	double *b = (double *) malloc( (nf-1) * sizeof(double) );
	double *c = (double *) malloc( (nf-1) * sizeof(double) );
	double *x = (double *) malloc( nf * sizeof(double) );
	double *xx = (double *) malloc( nf * sizeof(double) );

	double *nj = (double *) malloc( nf * sizeof(double) );
	double *rhs = (double *) malloc( nf * sizeof(double) );

	//double *rhs0 = (double *) malloc( nf * sizeof(double) );

	int i,fc,ii;
	double y;
	double V_lw, V_uw;
	int num;
	double pi = acos(-1);



  int NNr;//(*box).N[1];//500;
  int NNphi;

   int nn = 50;
    NNr = nn;
    NNphi = nn;
    //double *xi = (double *) malloc( nn * sizeof(double) );
    //double *wi = (double *) malloc( nn * sizeof(double) );
    double bb, aa, cc, dd;

      double *rr = (double *) malloc( NNr * sizeof(double) );
      double *phi = (double *) malloc( NNphi * sizeof(double) );
      double *ddphi = (double *) malloc( NNphi * sizeof(double) );
      double *delt_r = (double *) malloc( NNr * sizeof(double) );

   bb = 3.0*(*gas).eqdist;
   aa = 0.0*(*gas).eqdist;

   cc = 0.0;
   dd = pi;

	 double xi[50] = {-0.9988664,  -0.99403197, -0.98535408, -0.97286439, -0.95661096, -0.93665662,
										-0.91307856, -0.88596798, -0.85542977, -0.82158207, -0.78455583, -0.7444943,
										-0.70155247, -0.65589647, -0.60770293, -0.5571583 , -0.50445814, -0.44980633,
										-0.39341431, -0.33550025, -0.27628819, -0.21600724, -0.15489059, -0.0931747,
										-0.03109834,  0.03109834,  0.0931747 ,  0.15489059,  0.21600724,  0.27628819,
										0.33550025 , 0.39341431 , 0.44980633 , 0.50445814 , 0.5571583  , 0.60770293,
										0.65589647 , 0.70155247 , 0.7444943  , 0.78455583 , 0.82158207 , 0.85542977,
										0.88596798 , 0.91307856 , 0.93665662 , 0.95661096 , 0.97286439 , 0.98535408,
										0.99403197 , 0.9988664};
	 double wi[50] = {0.00290862, 0.0067598,  0.01059055, 0.01438082, 0.01811556, 0.02178024,
	 0.02536067, 0.02884299, 0.03221373, 0.03545984, 0.03856876, 0.04152846,
	 0.0443275 , 0.04695505, 0.04940094, 0.0516557 , 0.05371062, 0.05555774,
	 0.05718993, 0.05860085, 0.05978506, 0.06073797, 0.0614559 , 0.06193607,
	 0.06217662, 0.06217662, 0.06193607, 0.0614559 , 0.06073797, 0.05978506,
	 0.05860085, 0.05718993, 0.05555774, 0.05371062, 0.0516557 , 0.04940094,
	 0.04695505, 0.0443275 , 0.04152846, 0.03856876, 0.03545984, 0.03221373,
	 0.02884299, 0.02536067, 0.02178024, 0.01811556, 0.01438082, 0.01059055,
	 0.0067598 , 0.00290862};

double n_out;
if( (*box).problem!=inverted )
	n_out = 0.5*(cells[0].n+cells[(*box).N[1]-1].n);

  // approximate number density on the faces
	for(i=0; i<nf; i++){
		if(i==0){
			if( (*box).problem == inverted ){
				nj[i] = cells[0].n;
			}
			else{
				nj[i] = cells[0].n;//0.5*( n_out + cells[0].n );
			}
		}
		else if(i==(*box).N[1]){
			if( (*box).problem == inverted ){
				nj[i] = cells[(*box).N[1]-1].n;
			}
			else{
				nj[i] = cells[(*box).N[1]-1].n;//0.5*( n_out+cells[(*box).N[1]-1].n );//0.01*(*gas).n;
			}
		}
		else{
			nj[i] = (cells[i].n+cells[i-1].n)/2.0;
		}
	}

	// Build the tridiagonal matrix
	for(i=0; i<nf; i++){
		a[i] = (2.0+pow( (*box).delta_dim[1]*(*gas).bst, 2.0) );
	}
	for(i=0; i<nf-1; i++){
		b[i] = -1.0;
		c[i] = -1.0;
	}
	for(i=0; i<nf; i++){
		rhs[i] = nj[i]*pow((*box).delta_dim[1],2.0);
	}


  for(i=0; i<NNr; i++){
		 rr[i] = (bb-aa)/2.0*xi[i]+(aa+bb)/2.0;
		 delt_r[i] = wi[i]*(bb-aa)/2.0;
  }
  for(i=0; i<NNr; i++){
		phi[i] = (dd-cc)/2.0*xi[i]+(cc+dd)/2.0;
		ddphi[i] = wi[i]*(dd-cc)/2.0;
  }

double radial_part;
double n = 0.0;

V_lw = 0.0; V_uw = 0.0;
	for(fc = 0; fc<nf; fc= fc+nf-1){
		for(i=0; i<NNr; i++){
			radial_part = exp(-(*gas).bst*rr[i])/(2.0)*rr[i]*delt_r[i];
			for(ii=0; ii<NNphi; ii++){
				y = fc*(*box).delta_dim[1] + rr[i]*cos(phi[ii]);
				num = floor(y/(*box).delta_dim[1]);
				if(num>-1 && num<(*box).N[1]){
					n = cells[num].n;
				}
				if(num<0){
					if( (*box).problem == inverted ){
						n = cells[0].n;
					}
					else{
							n =  cells[0].n;//n_out;
					}
				}
				else if(num > (*box).N[1]-1){
					if( (*box).problem == inverted ){
						n = cells[(*box).N[1]-1].n;
					}
					else{
						n = cells[(*box).N[1]-1].n;//n_out;
					}
				}
				if(fc==0)
					V_lw += radial_part*n*sin(phi[ii])*ddphi[ii];
				if(fc==nf-1)
					V_uw += radial_part*n*sin(phi[ii])*ddphi[ii];
			}
		}
	}


// my code, working with dirichlet BC

	c[0] = 0.0;
	a[0] = 1.0;
	rhs[0] = V_lw;//- nj[0]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);

	a[nf-1]    = 1.0;
	rhs[nf-1]  = V_uw;//- nj[nc-1]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);
	b[nf-2]    = 0.0;


// let's try zero dirichlet
/*
	c[0] = 0.0;
	a[0] = 1.0;
	rhs[0] = 0.0;//- nj[0]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);

	a[nf-1]    = 1.0;
	rhs[nf-1]  = 0.0;//- nj[nc-1]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);
	b[nf-2]    = 0.0;
*/

// Let's try Neumann
/*
c[0] = -1.0;
a[0] = 1.0;
rhs[0] = 0.0;//- nj[0]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);

a[nf-1]    = 1.0;
rhs[nf-1]  = 0.0;//- nj[nc-1]*( 1.0/(*gas).bst/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) + (*gas).eqdist/(*gas).bst*exp(-(*gas).bst*(*gas).eqdist) - 1.0/(*gas).bst/(*gas).bst);
b[nf-2]    = -1.0;
*/
	// solve the system
	//for(i=0; i<nf; i++){
	//	rhs0[i] = rhs[i];
	//}
	Thompson_3d(a, b, c, x, rhs, nf);

	printf("1D:\n");
  for(i=0; i<nf; i++){
 		printf("x[%d]=%e\n", i,x[i]);
   }

	/*
	FILE *fp;
	fp = fopen("sol_poisson.txt", "w");
	for(i=0; i<nf; i++){
		fprintf(fp,"%d   %e   %e\n", i, x[i], nj[i]);
	}
	fclose(fp);

	double *res = (double *) malloc( nf * sizeof(double) );
	for(i=0; i<nf; i++){
		if(i==0){
			res[i] = x[i] - 10000.0;
		}
		else if(i==nf-1){
			res[i] = x[i] - 0.0;
		}
		else{
			res[i] = rhs0[i] - (2.0+pow( (*box).delta_dim[1]*(*gas).bst, 2.0) )*x[i] + x[i-1] + x[i+1];
		}
	}

	double l2norm_res=0.0, l2norm_x = 0.0;
	for(i=0; i<nf; i++){
		l2norm_x += pow( x[i], 2.0 );
		l2norm_res += pow( res[i], 2.0);
	}
	l2norm_x = sqrt(l2norm_x);
	l2norm_res = sqrt(l2norm_res);
	printf("||err||_2/||x||_2 = %e \n", l2norm_res/l2norm_x);
	*/
// Direct way of computing int_0^sigma (attraction)
   bb = 1.0*(*gas).eqdist; aa=0.0;
   for(i=0; i<NNr; i++){
		 rr[i] = (bb-aa)/2.0*xi[i]+(aa+bb)/2.0;
		 delt_r[i] = wi[i]*(bb-aa)/2.0;
   }

  for(fc = 0; fc<nf; fc++)
		xx[fc] = 0.0;

/*
  for(fc = 0; fc<nf; fc++){
		for(i=0; i<NNr; i++){
			radial_part = exp(-(*gas).bst*rr[i])*rr[i]*delt_r[i]/2.0;
			for(ii=0; ii<NNphi; ii++){
				y = fc*(*box).delta_dim[1] + rr[i]*cos(phi[ii]);
				num = floor(y/(*box).delta_dim[1]);
				if(num>-1 && num<(*box).N[1]){
					n = cells[num].n;
				}
				else if(num<0)
					n = cells[0].n;
				else if(num > (*box).N[1]-1)
					n = cells[(*box).N[1]-1].n;
				xx[fc] += radial_part*n*sin(phi[ii])*ddphi[ii];
		 }
	 }
 }
*/

 // analytic way of computing int_0^sigma (attraction)

 double dummy = 1.0-((*gas).sigma*(*gas).bst+1.0)*exp( -(*gas).bst*(*gas).sigma );
 dummy = dummy/(*gas).bst/(*gas).bst;
 for(i=0; i<nf; i++){
		xx[i] =  nj[i]*dummy;
 }


// storing potential of the left face for each cell, postproc purposes
for(num = 0; num<(*box).N[1]; num++){
	 cells[num].UU = x[num];///(*gas).m*(*gas).ast;
	 cells[num].UU0 = nj[num];//xx[num];///(*gas).m*(*gas).ast;

	 //cells[num].UUR = x[num+1]/(*gas).m*(*gas).ast;
	 cells[num].UU0R = xx[num+1]/(*gas).m*(*gas).ast;
 }
/// the following lines are neceesary. uncomment them ASAP

 for(i=0; i<nf; i++){
		x[i] = x[i] - xx[i];
  }


// test
/*
for(i=0; i<nf; i++){
		x[i] = xx[i];
  }
*/

	// calculate the force at y, for cell num
	for(num = 0; num<(*box).N[1]; num++){
		cells[num].UUR = (nj[num+1]-nj[num])/(*box).delta_dim[1];
		cells[num].F2 = -(x[num+1]-x[num])/(*box).delta_dim[1];
		cells[num].F2 = cells[num].F2/(*gas).m*(*gas).ast;
	}
   for(i=0; i<(*gas).N; i++){
	   num = index[i];
	   U2[i] = U2[i] + cells[num].F2*(*gas).delta_t;
   }

   free(x); free(xx); free(rhs); free(a); free(b); free(c);
   free(rr); free(phi); free(delt_r); free(nj);
}
void vlasov_taking_the_integral(double *U2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *x1, double *x2, double *x3){
  std::random_device rd;
  std::mt19937 gen(rd());


  std::uniform_real_distribution<> uniform(0.0, 1.0);

      int nc = (*box).N[1];
      int nf = (*box).N[1] + 1;
      int fc, id;
      double *V = (double *) malloc( nf * sizeof(double) );
      double *N = (double *) malloc( nf * sizeof(double) );
      int i, num, ii;
      double r;
      double pi = acos(-1);

      double y;
      double V_r,V_phi;
      int j, j0, neigh;
      double r0 = (*gas).eqdist;
      double cut_off = 3.0*r0;
      double epsilon = 1e-15;
      neigh = floor((cut_off)/(*box).delta_dim[1])+1;

      int NNr;
      int NNphi;

	//	rr[i] = 1.122*(*gas).sigma + (cut_off-1.122*(*gas).sigma)*uniform(gen);
	//	phi[i] = 0.0 + pi*uniform(gen);

      double sum_c = 0.0;
double n_out;
if( (*box).problem == inverted ||  (*box).problem == evaporation)
 	n_out = 0.5*(cells[0].n+cells[(*box).N[1]-1].n);
else if( (*box).problem == vacuum  ){
	n_out = 0.5*(cells[0].n+cells[(*box).N[1]-1].n);
}
/*
      for(i=0; i<NNr; i++)
		sum_c += pow(c,i);
      delt_r[0] = (cut_off-r0)/sum_c;
      for(i=1; i<NNr; i++)
		delt_r[i] = delt_r[0]*pow(c,i);
      rr[0] = r0 + delt_r[0]/2.0;
      for(i=1; i<NNr; i++)
		rr[i] = rr[i-1]+(delt_r[i-1]+delt_r[i])/2.0;
*/

			 int nn = 50;
    NNr = nn;
    NNphi = nn;

    double b, a, c, d;

      double *rr = (double *) malloc( NNr * sizeof(double) );
      double *phi = (double *) malloc( NNphi * sizeof(double) );
      double *ddphi = (double *) malloc( NNphi * sizeof(double) );
      double *delt_r = (double *) malloc( NNr * sizeof(double) );

   b = 3.0*(*gas).eqdist;
   a = 1.0*(*gas).eqdist;

   c = 0.0;
   d = pi;


   double xi[50] = {-0.9988664,  -0.99403197, -0.98535408, -0.97286439, -0.95661096, -0.93665662,
                    -0.91307856, -0.88596798, -0.85542977, -0.82158207, -0.78455583, -0.7444943,
                    -0.70155247, -0.65589647, -0.60770293, -0.5571583 , -0.50445814, -0.44980633,
                    -0.39341431, -0.33550025, -0.27628819, -0.21600724, -0.15489059, -0.0931747,
                    -0.03109834,  0.03109834,  0.0931747 ,  0.15489059,  0.21600724,  0.27628819,
                    0.33550025 , 0.39341431 , 0.44980633 , 0.50445814 , 0.5571583  , 0.60770293,
                    0.65589647 , 0.70155247 , 0.7444943  , 0.78455583 , 0.82158207 , 0.85542977,
                    0.88596798 , 0.91307856 , 0.93665662 , 0.95661096 , 0.97286439 , 0.98535408,
                    0.99403197 , 0.9988664};
   double wi[50] = {0.00290862, 0.0067598,  0.01059055, 0.01438082, 0.01811556, 0.02178024,
   0.02536067, 0.02884299, 0.03221373, 0.03545984, 0.03856876, 0.04152846,
   0.0443275 , 0.04695505, 0.04940094, 0.0516557 , 0.05371062, 0.05555774,
   0.05718993, 0.05860085, 0.05978506, 0.06073797, 0.0614559 , 0.06193607,
   0.06217662, 0.06217662, 0.06193607, 0.0614559 , 0.06073797, 0.05978506,
   0.05860085, 0.05718993, 0.05555774, 0.05371062, 0.0516557 , 0.04940094,
   0.04695505, 0.0443275 , 0.04152846, 0.03856876, 0.03545984, 0.03221373,
   0.02884299, 0.02536067, 0.02178024, 0.01811556, 0.01438082, 0.01059055,
   0.0067598 , 0.00290862};

/*
    int nn = 10;
		double *xi = (double *) malloc( nn * sizeof(double) );
		double *wi = (double *) malloc( nn * sizeof(double) );
   xi[0] = -0.1488743389816312;
   xi[1] =  0.1488743389816312;
   xi[2] = -0.4333953941292472;
   xi[3] =  0.4333953941292472;
   xi[4] = -0.6794095682990244;
   xi[5] =  0.6794095682990244;
   xi[6] = -0.8650633666889845;
   xi[7] =  0.8650633666889845;
   xi[8] = -0.9739065285171717;
   xi[9] =  0.9739065285171717;

   wi[0] = 0.2955242247147529;
   wi[1] = 0.2955242247147529;
   wi[2] = 0.2692667193099963;
   wi[3] = 0.2692667193099963;
   wi[4] = 0.2190863625159820;
   wi[5] = 0.2190863625159820;
   wi[6] = 0.1494513491505806;
   wi[7] = 0.1494513491505806;
   wi[8] = 0.0666713443086881;
   wi[9] = 0.0666713443086881;
*/
/*
   xi[0] =  0.0;
   xi[1] =  (1.0/3.0)*sqrt( 5.0-2.0*sqrt(10.0/7.0) );
   xi[2] = -(1.0/3.0)*sqrt( 5.0-2.0*sqrt(10.0/7.0) );
   xi[3] =  (1.0/3.0)*sqrt( 5.0+2.0*sqrt(10.0/7.0) );
   xi[4] = -(1.0/3.0)*sqrt( 5.0+2.0*sqrt(10.0/7.0) );
   wi[0] = 128.0/225.0;
   wi[1] = ( 322.0+13.0*sqrt(70.0) )/900.0;
   wi[2] = ( 322.0+13.0*sqrt(70.0) )/900.0;
   wi[3] = ( 322.0-13.0*sqrt(70.0) )/900.0;
   wi[4] = ( 322.0-13.0*sqrt(70.0) )/900.0;
*/
/*
      for (i=0; i<NNr; ++i){
		rr[i] = r0 + ( (cut_off-r0)*(i+0.5) )/( 1.0* NNr);
		delt_r[i] = (cut_off-r0)/( 1.0* NNr);
      }
*/
   for(i=0; i<NNr; i++){
	rr[i] = (b-a)/2.0*xi[i]+(a+b)/2.0;
	delt_r[i] = wi[i]*(b-a)/2.0;
   }
   for(i=0; i<NNr; i++){
	phi[i] = (d-c)/2.0*xi[i]+(c+d)/2.0;
	ddphi[i] = wi[i]*(d-c)/2.0;
   }
//      for (i=0; i<NNphi; ++i)
//		phi[i] = 0.0 + (pi*(i))/(1.0*NNphi);

double radial_part = 0.0;
//double coeff = (cut_off-1.122*(*gas).sigma)/(1.0*NNr)*pi/(2.0*NNphi);
double coeff = pi/(2.0*NNphi);
double n = 0.0;
double int_theta = 2.0*pi;
double dphi = pi/(1.0*NNphi);
int num1;

/*
for(num=0; num<(*box).N[1]; num++){
	cells[num].F2 = 0.0;
	for(i=0; i<NNr; i++){
		radial_part = 2.0*pi*6.0*pow((*gas).eqdist/rr[i],5.0)*(*gas).eqdist*delt_r[i];
		for(ii=0; ii<NNphi; ii++){
			y = cells[num].cell_center[1] + rr[i]*cos(phi[ii]);
			num1 = floor(y/(*box).delta_dim[1]);
			if(num1>-1 && num1<(*box).N[1]){
				n = cells[num1].n;
			}
			else
				n = n_out;
			cells[num].F2 += radial_part*n*0.5*sin(2.0*phi[ii])*ddphi[ii];
		}
	}
}
*/

// working, for the paper was used
	for(fc = 0; fc<nf; fc++){
		V[fc] = 0.0;
		V_r=0.0; V_phi=0.0;
		for(i=0; i<NNr; i++){
			//radial_part = -(pow((*gas).sigma/rr[i],12.0) - pow((*gas).sigma/rr[i],6.0) )*rr[i]*rr[i]*int_theta*delt_r[i];
			radial_part = -( - pow((*gas).sigma/rr[i],6.0) )*rr[i]*rr[i]*int_theta*delt_r[i];
			for(ii=0; ii<NNphi; ii++){
				y = fc*(*box).delta_dim[1] + rr[i]*cos(phi[ii]);
				/*
				if( y > (*box).Len[1])
					y = y - (*box).Len[1];
				if( y < 0.0)
					y = y + (*box).Len[1];
				*/
				num = floor(y/(*box).delta_dim[1]);
				if(num>-1 && num<(*box).N[1]){
					n = cells[num].n;
				}
				else if(num<0){
					//if( (*box).problem == inverted ||  (*box).problem == evaporation)
					if(  (*box).problem == evaporation ){
						//n = (*gas).nv1;
						n = n_out;
					}
					else if( (*box).problem == inverted ){
						n  =  cells[0].n;
					}
					else{
						n  =  cells[0].n;
					}
				}
				else if(num > (*box).N[1]-1){
					if(  (*box).problem == evaporation ){
						//n = (*gas).nv2;
						n = n_out;
					}
					else if( (*box).problem == inverted ){
						n =  cells[(*box).N[1]-1].n;
					}
					else{
						n =  cells[(*box).N[1]-1].n;
					}
				}
				V[fc] += -radial_part*n*sin(phi[ii])*ddphi[ii];
			}
		}
	}

/* //working
	for(fc = 0; fc<nf; fc++){
		V[fc] = 0.0;
		V_r=0.0; V_phi=0.0;
		for(i=0; i<NNr; i++){
			radial_part = exp(-(*gas).bst*rr[i])*rr[i]*coeff*delt_r[i];
			for(ii=0; ii<NNphi; ii++){
				y = fc*(*box).delta_dim[1] + rr[i]*cos(phi[ii]);
				num = floor(y/(*box).delta_dim[1]);
				if(num>-1 && num<(*box).N[1]){
					n = cells[num].n_smooth;
				}
				else if(num<0)
					n = n_out;
				else if(num > (*box).N[1]-1)
					n = n_out;
				V[fc] += radial_part*n*sin(phi[ii]);
			}
		}
	}
*/
/*
	for(fc = 0; fc<nf; fc++){
		V[fc] =  (*gas).n*( exp(-(*gas).bst*1.122*(*gas).sigma)/(*gas).bst*(1.122*(*gas).sigma + 1.0/(*gas).bst)
				    -exp(-(*gas).bst*cut_off)/(*gas).bst*(cut_off + 1.0/(*gas).bst) );
	}
*/
/*

            	for(j=num-neigh; j<=num+neigh+1; j++){
			if(j>-1 && j<(*box).N[1])
				j0 = j;
			else if(j<0)
				j0 = 0;
			else if(j > (*box).N[1]-1)
				j0 = (*box).N[1]-1;
	                	for(i=0; i<cells[j0].num_inside; i++){
					id = cells[j0].indices_inside[i];
					//r = sqrt( pow(x1[id]-0.5*(*box).Len[0],2.0) + pow(x2[id]+(j-j0)*(*box).delta_dim[1]-fc*(*box).delta_dim[1],2.0) + pow(x3[id]-0.5*(*box).Len[2],2.0));
					r = fabs(x2[id]+(j-j0)*(*box).delta_dim[1]-fc*(*box).delta_dim[1]);
					if(r > 1.122*(*gas).sigma && r < cut_off){
						V[fc] += exp(-(*gas).bst*r)*r*cells[j0].n;
						N[fc] += 1.0;
					}
				}
		}
	}
	for(fc = 0; fc<nf; fc++){
		if(N[fc]>0.5)
			V[fc] = V[fc]*(cut_off-1.122*(*gas).sigma)/N[fc];
		else
			printf("Problem: no particle was found in  1.122*sigma<r<cut_off\n\n");
	}
*/
/*
	for(fc = 0; fc<nf; fc++){
		V[fc] = 0.0;
		N[fc] = 0.0;
		num = fc-1;
            	for(j=num-neigh; j<=num+neigh+1; j++){
			if(j>-1 && j<(*box).N[1])
				j0 = j;
			else if(j<0)
				j0 = 0;
			else if(j > (*box).N[1]-1)
				j0 = (*box).N[1]-1;
	                	for(i=0; i<cells[j0].num_inside; i++){
					id = cells[j0].indices_inside[i];
					//r = sqrt( pow(x1[id]-0.5*(*box).Len[0],2.0) + pow(x2[id]+(j-j0)*(*box).delta_dim[1]-fc*(*box).delta_dim[1],2.0) + pow(x3[id]-0.5*(*box).Len[2],2.0));
					r = fabs(x2[id]+(j-j0)*(*box).delta_dim[1]-fc*(*box).delta_dim[1]);
					if(r > 1.122*(*gas).sigma && r < cut_off){
						V[fc] += exp(-(*gas).bst*r)*r*cells[j0].n;
						N[fc] += 1.0;
					}
				}
		}
	}
	for(fc = 0; fc<nf; fc++){
		if(N[fc]>0.5)
			V[fc] = V[fc]*(cut_off-1.122*(*gas).sigma)/N[fc];
		else
			printf("Problem: no particle was found in  1.122*sigma<r<cut_off\n\n");
	}
*/

/*
	for(fc = 0; fc<nf; fc++){
		V[fc] = 0.0;

		//lower cells
		num = fc-1;

            	for(j=num; j>=num-neigh; j--){
                  y=(j-num-0.5)*(*box).delta_dim[1];
                  if(j<0)
                        nj = 0.0;//cells[0].n;
                  else if(j>(*box).N[1]-1)
                        nj = 0.0;//cells[(*box).N[1]-1].n;
                  else
                        nj = cells[j].n;
                  if(fabs(y)<1.123*(*gas).sigma)
                        V[fc] += 0.0;//4.0*(*gas).epsilon*pi*pow( (*gas).sigma,3.0 )*nj/(*gas).m;
                  else
                        V[fc] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,3.0 )*nj*pow((*gas).sigma/y,4.0)/(*gas).m;
            	}

		//upper cells
		num = fc;

            	for(j=num; j<=num+neigh; j++){
                  y=(j-num+0.5)*(*box).delta_dim[1];
                  if(j<0)
                        nj = 0.0;//cells[0].n;
                  else if(j>(*box).N[1]-1)
                        nj = 0.0;//cells[(*box).N[1]-1].n;
                  else
                        nj = cells[j].n;
                  if(fabs(y)<1.123*(*gas).sigma)
                        V[fc] += 0.0;//4.0*(*gas).epsilon*pi*pow( (*gas).sigma,3.0 )*nj/(*gas).m;
                  else
                        V[fc] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,3.0 )*nj*pow((*gas).sigma/y,4.0)/(*gas).m;
            	}
		y = fc*(*box).delta_dim[1] + 1.123*(*gas).sigma;
		V[fc] += 4.0*(*gas).epsilon*( pow((*gas).sigma/y,12.0)- pow((*gas).sigma/y,6.0) )/(*gas).m;
	}
*/
	for(num=0; num<(*box).N[1]; num++){
		fc = num+1;
		cells[num].F2 = -(V[fc]-V[fc-1])/(*box).delta_dim[1];
		cells[num].F2 = cells[num].F2/(*gas).m*(*gas).phi;
	//	cells[num].F2 = 4.0*(*gas).epsilon*cells[num].F2/(*gas).m;
	}
double dummyF;
      for(id=0; id<(*gas).N; id++){
	    num = index[id];
	    U2[id] = U2[id] + cells[num].F2*(*gas).delta_t;
	  //  if(x2[id]>2.0*(*gas).sigma){
	//	    dummyF = 24.0*(*gas).epsilon/x2[id]*( -2.0*pow((*gas).sigma/x2[id],12.0) + pow((*gas).sigma/x2[id],6.0) )/(*gas).m;
	  //          U2[id] = U2[id] - dummyF*(*gas).delta_t;
	  //  }
      }

	free(V); free(N); free(rr); free(phi);  free(delt_r); free(ddphi);
	// free(xi); free(wi);
}
void vlasov_term(double *U2, double *V, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *x1, double *x2, double *x3){
      int i, num, id;
      double F;
      double r;
      int id2;
      double epsilon = 4.47478200e-21;
    if( (*gas).long_range!= 0 ){
      //if( (*gas).long_range == 1 )
      //    calculate_potential(V, gas, box, cells);
      for(id=0; id<(*gas).N; id++){
           num = index[id];
           F = 0.0;

        //   }
  //         F = F*epsilon/((*gas).Fn*(*gas).m);
  //         F = 0.0;
  /*
           G = 0.0;
           if(step > 1){
             if(num == 0){
             G = (*gas).Fn*cells[num].M_n[3]/(*gas).delta_t/((*box).delta_dim[0]*(*box).delta_dim[2]);
             G = G - (*gas).n;
             }
             else if(num == (*box).N[1]-1){
             G = (*gas).n;
             G = G - (*gas).Fn*cells[num].M_n[2]/(*gas).delta_t/((*box).delta_dim[0]*(*box).delta_dim[2]);
             }
             else{
             G = (*gas).Fn*cells[num].M_n[3]/(*gas).delta_t/((*box).delta_dim[0]*(*box).delta_dim[2]);
             G = G - (*gas).Fn*cells[num].M_n[2]/(*gas).delta_t/((*box).delta_dim[0]*(*box).delta_dim[2]);}
             G = G/((*box).delta_dim[1]/2.0);
             G = G*pow((*gas).sigma,3.0)*4.47478200e-21/(*gas).m;
          }
          */
  /*
          G = 0.0;
           if(num == 0){
             F = (cells[num+1].n-(*gas).n+cells[num+2].n/128.0-(*gas).n/128.0)/(*box).delta_dim[1];
           }
           else if(num == 1){
                 F = (cells[num+1].n-cells[num-1].n+cells[num+2].n/128.0-(*gas).n/128.0)/(*box).delta_dim[1];
           }
           else if(num == (*box).N[1]-2){
             F = (cells[num+1].n-cells[num-1].n+(*gas).n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
           else if(num == (*box).N[1]-1){
             F = ((*gas).n-cells[num-1].n+(*gas).n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
           else{
             F = (cells[num+1].n-cells[num-1].n+cells[num+2].n/128.0-cells[num-2].n/128.0)/(*box).delta_dim[1];}
  //4.47478200e-21*
            F = 4.0*3.14*F*pow((*gas).sigma,3.0)*4.47478200e-21/(*gas).m;
  */
            //F = F+G;
  //          F = 0.0;
  /*
        F = 0.0;
           if(num == 0)
             F = (cells[num+1].n-cells[0].n)/(2.0*(*box).delta_dim[1]);
           else if(num == (*box).N[1]-1)
             F = (cells[(*box).N[1]-1].n-cells[num].n)/(2.0*(*box).delta_dim[1]);
           else
             F = (cells[num+1].n-cells[num-1].n)/(2.0*(*box).delta_dim[1]);

            F = 2.0*F*(4.0*3.14)*pow((*gas).sigma,3.0)*4.47478200e-21/(*gas).m;
  */
        //   F = 0.0;


        if((*gas).long_range==1){
              i = num+1;
              F =  2.0*(V[i+1]-V[i-1])/(*box).delta_dim[1]/2.0;
        }

	F = 24.0*(*gas).epsilon/x2[id]*( -2.0*pow((*gas).sigma/x2[id],12.0) + pow((*gas).sigma/x2[id],6.0) )/(*gas).m;
            //Mp2[id] = Mp2[id] - F*(*gas).delta_t;
            U2[id] = U2[id] - F*(*gas).delta_t;
      }
     }
}

void calculate_potential(double *V, struct GAS *gas, struct BOX *box, struct CELLS *cells){
      int i, num;
      double a,b,c,x,y;
      double pi = acos(-1);
      double p_m1 = -0.5*(*box).delta_dim[1];
      double p_m2 = -1.5*(*box).delta_dim[1];
      double p_m3 = -2.5*(*box).delta_dim[1];
      double n_m1 = cells[0].n;
      double n_m2 = cells[0].n;
      double n_m3 = cells[0].n;
      double r;
      double p_N   = (*box).Len[1] + 0.5*(*box).delta_dim[1];
      double p_Np1 = (*box).Len[1] + 1.5*(*box).delta_dim[1];
      double p_Np2 = (*box).Len[1] + 2.5*(*box).delta_dim[1];
      double n_N   = cells[(*box).N[1]-1].n;
      double n_Np1 = cells[(*box).N[1]-1].n;
      double n_Np2 = cells[(*box).N[1]-1].n;

      double nj,yj;
      int j, neigh;
      double cut_off = 2.0*(*gas).sigma;
      double epsilon = 1e-15;
      neigh = floor((cut_off+epsilon)/(*box).delta_dim[1]);
      for(num=-1; num<=(*box).N[1]; num++){
            // current cell
            j=num;
            if(j<0)
                  nj = cells[0].n;
            else if(j>(*box).N[1]-1)
                  nj = cells[(*box).N[1]-1].n;
            else
                  nj = cells[j].n;
            V[num+1] = 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,2.0 )*(*box).delta_dim[1]*nj/(*gas).m;

            // upper cells
            for(j=num+1; j<=num+neigh; j++){
                  y=(j-num)*(*box).delta_dim[1];
                  if(j<0)
                        nj = cells[0].n;
                  else if(j>(*box).N[1]-1)
                        nj = cells[(*box).N[1]-1].n;
                  else
                        nj = cells[j].n;
                  if(fabs(y)<(*gas).sigma)
                        V[num+1] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,2.0 )*nj/(*gas).m*(*box).delta_dim[1];
                  else
                        V[num+1] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,2.0 )*nj*pow((*gas).sigma/y,4.0)/(*gas).m*(*box).delta_dim[1];
            }

            for(j=num-1; j>=num-neigh; j--){
                  y=(j-num)*(*box).delta_dim[1];
                  if(j<0)
                        nj = cells[0].n;
                  else if(j>(*box).N[1]-1)
                        nj = cells[(*box).N[1]-1].n;
                  else
                        nj = cells[j].n;
                  if(fabs(y)<(*gas).sigma)
                        V[num+1] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,2.0 )*nj/(*gas).m*(*box).delta_dim[1];
                  else
                        V[num+1] += 4.0*(*gas).epsilon*pi*pow( (*gas).sigma,2.0 )*nj*pow((*gas).sigma/y,4.0)/(*gas).m*(*box).delta_dim[1];
            }

      }
/*
      for(num=-1; num<=(*box).N[1]; num++){
            i = num+1;
            if(num == -1){
                  x = p_m1;

                  a = -(cells[num+1].n -cells[num+2].n)/(*box).delta_dim[1];
                  b = cells[num+1].n- a*(*gas).sigma;//(cells[num+1].cell_center[1]-x);
                  r = 2.0*(*gas).sigma;
                  V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r =  (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*r+4.0*b)/(20*pow(r,5.0));

                  a = 0.0;//(n_m2 - n_m3)/(*box).delta_dim[1];
                  b = n_m2;// - a*(*gas).sigma;//(x-p_m2);
                  y = p_m2;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(x-y)+4.0*b)/(20*pow(x-y,5.0));
                  y = p_m3;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(x-y)+4.0*b)/(20*pow(x-y,5.0));
            }
            else if(num == 0){
                  x = cells[num].cell_center[1];

                  a = -(cells[num+1].n -cells[num+2].n)/(*box).delta_dim[1];
                  b = cells[num+1].n- a*(*gas).sigma;//(cells[num+1].cell_center[1]-x);
                  r = 2.0*(*gas).sigma;
                  V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r =  (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*r+4.0*b)/(20*pow(r,5.0));

                  a = 0.0;//(n_m1 - n_m2)/(*box).delta_dim[1];
                  b = n_m1;// - a*(*gas).sigma;//(x-p_m1);
                  y = p_m1;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(x-y)+4.0*b)/(20*pow(x-y,5.0));
                  y = p_m2;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(x-y)+4.0*b)/(20*pow(x-y,5.0));
            }
            else if(num == 1){
                 x = cells[num].cell_center[1];

                 a = -(cells[num+1].n -cells[num+2].n)/(*box).delta_dim[1];
                 b = cells[num+1].n- a*(*gas).sigma;//(cells[num+1].cell_center[1]-x);
                 r = 2.0*(*gas).sigma;
                 V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                 r =  (*gas).sigma;
                 V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*r+4.0*b)/(20*pow(r,5.0));

                 a = -(cells[num-1].n - n_m1)/(*box).delta_dim[1];
                 b = cells[num-1].n - a*(*gas).sigma;//(x-cells[num-1].cell_center[1]);
                 r = 2.0*(*gas).sigma;
                 V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                 r = (*gas).sigma;
                 V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));

            }
            else if(num == (*box).N[1]-2){
                  x = cells[num].cell_center[1];

                  a = -(cells[num+1].n -n_N)/(*box).delta_dim[1];
                  b = cells[num+1].n- a*(*gas).sigma;//(cells[num+1].cell_center[1]-x);
                  r = 2.0*(*gas).sigma;
                  V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r =  (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*r+4.0*b)/(20*pow(r,5.0));

                  a = -(cells[num-1].n - cells[num-2].n)/(*box).delta_dim[1];
                  b = cells[num-1].n - a*(*gas).sigma;//(x-cells[num-1].cell_center[1]);
                  r = 2.0*(*gas).sigma;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r = (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
            }
            else if(num == (*box).N[1]-1){
                  x = cells[num].cell_center[1];

                  a = 0.0;//(n_Np1 - n_N)/(*box).delta_dim[1];
                  b = n_N;// - a*(*gas).sigma;//(p_N-x);
                  y = p_Np1;
                  V[i] = -pow((*gas).sigma,8.0)*(5.0*a*(y-x)+4.0*b)/(20*pow(y-x,5.0));
                  y = p_N;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(y-x)+4.0*b)/(20*pow(y-x,5.0));

                  a = -(cells[num-1].n - cells[num-2].n)/(*box).delta_dim[1];
                  b = cells[num-1].n - a*(*gas).sigma;//(x-cells[num-1].cell_center[1]);
                  r = 2.0*(*gas).sigma;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r = (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
            }
            else if(num == (*box).N[1]){
                  x = p_N;

                  a = 0.0;//(n_Np2 - n_Np1)/(*box).delta_dim[1];
                  b = n_Np1;// - a*(*gas).sigma;//(p_Np1-x);
                  y = p_Np2;
                  V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(y-x)+4.0*b)/(20*pow(y-x,5.0));
                  y = p_Np1;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(y-x)+4.0*b)/(20*pow(y-x,5.0));

                  a = -(cells[num-1].n - cells[num-2].n)/(*box).delta_dim[1];
                  b = cells[num-1].n - a*(*gas).sigma;//(x-cells[num-1].cell_center[1]);
                  r = 2.0*(*gas).sigma;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r = (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
            }
            else{
                  x = cells[num].cell_center[1];

                  a = -(cells[num+1].n -cells[num+2].n)/(*box).delta_dim[1];
                  b = cells[num+1].n- a*(*gas).sigma;//(cells[num+1].cell_center[1]-x);
                  r = 2.0*(*gas).sigma;
                  V[i]  = -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r =  (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*r+4.0*b)/(20*pow(r,5.0));

                  a = -(cells[num-1].n - cells[num-2].n)/(*box).delta_dim[1];
                  b = cells[num-1].n - a*(*gas).sigma;//(x-cells[num-1].cell_center[1]);
                  r = 2.0*(*gas).sigma;
                  V[i] += -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
                  r = (*gas).sigma;
                  V[i] -= -pow((*gas).sigma,8.0)*(5.0*a*(r)+4.0*b)/(20*pow(r,5.0));
            }
            V[i] = V[i]*4.0*(*gas).epsilon/(*gas).m;
      }
      V[0] = V[2];
      V[(*box).N[1]+1] = V[(*box).N[1]-1];
*/
      if ( (*box).step > (*gas).times*(*box).after){
            for(num=0; num<(*box).N[1]; num++){
                  i = num+1;
                  cells[num].G_sum += V[i];
            }
      }
}

void initialization(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells){

  std::random_device rd;
  std::mt19937 gen(rd());

//   std::uniform_int_distribution<> dis(1, 6);
  std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::uniform_real_distribution<> dis_xx(-1.0, 1.0);
int i,j,k, id;
int count = 0;
int num;

(*gas).count_SP_c = 0.0;  (*gas).count_reem_c=0.0;  (*gas).count_evap_c=0.0;
(*gas).count_SP_h = 0.0;  (*gas).count_reem_h=0.0;  (*gas).count_evap_h=0.0;
if((*box).problem == dcones){
	double Lx = (*box).Lx[0] + (*box).Lx[1] + (*box).Lx[2];
	int Nx = (*box).Nx[0] + (*box).Nx[1] + (*box).Nx[2];
	double yu,yd;
	double dx[3];
	int num;
	int Nc;
	dx[0] = (*box).Lx[0]/(1.0*(*box).Nx[0]);
	dx[1] = (*box).Lx[1]/(1.0*(*box).Nx[1]);
	dx[2] = (*box).Lx[2]/(1.0*(*box).Nx[2]);
	id = 0;
	(*gas).N = 0;
	for(i=0; i<Nx; i++)
	{
		for(k=0; k<(*box).Ny; k++){
			num = k*Nx + i;
			Nc = floor( cells[num].volume*(*gas).n/(*gas).Fn );
			(*gas).N += Nc;
			for(j=0; j<Nc; j++){
				yd = (*box).dy0*geometric_sum(k-1,  (*box).s);
				yu = (*box).dy0*geometric_sum(k,(*box).s);
				x2[id] = yd + dis_x(gen)*(yu-yd);
				if(i<(*box).Nx[0]){
					x1[id] = i*dx[0] + dis_x(gen)*dx[0];
				}
				else if(i>=(*box).Nx[0] && i<(*box).Nx[0]+(*box).Nx[1]){
					x1[id] = (*box).Lx[0] + (i-(*box).Nx[0])*dx[1]+ dis_x(gen)*dx[1];
				}
				else{
					x1[id] = (*box).Lx[0] + (*box).Lx[1] + (i-(*box).Nx[0]-(*box).Nx[1])*dx[2] + dis_x(gen)*dx[2];
				}
				id++;
			}
		}
	}
    	for ( i = 0; i < (*gas).N; ++i){
		if(x1[i]<(*box).Lx[0]){
		}
		else if(x1[i]<(*box).Lx[0]+(*box).Lx[1]){
			x2[i] = x2[i] + tan((*box).alpha[0])*(x1[i]-(*box).Lx[0]);
		}
		else{
			x2[i] = x2[i] + (*box).Lx[1]*tan((*box).alpha[0]) + tan((*box).alpha[1])*(x1[i]-(*box).Lx[0]-(*box).Lx[1]);
		}
		num = dcone_index(x1[i], x2[i], cells,  box, &j, &k);
    	}

}
else if((*box).problem == flatnose){
	int Nx = (*box).N1[0] + (*box).N1[1];
	int Ny = (*box).N2[0] + (*box).N2[1];
	double yu,yd;
	int Nc;
        double dummy;
	for(num=0; num<Nx*Ny; num++){
		if(cells[num].in == 1)
			dummy += cells[num].volume;
	}
	id = 0;
	(*gas).N = 0;
	for(i=0; i<Nx; i++)
	{
		for(k=0; k<Ny; k++){
			num = k*Nx + i;
			if(cells[num].in == 1){
				Nc = floor( cells[num].volume*(*gas).n/(*gas).Fn );
				(*gas).N += Nc;
				for(j=0; j<Nc; j++){
					yd = cells[num].dim[0];
					yu = cells[num].dim[2];
					x2[id] = yd + dis_x(gen)*(yu-yd);
					x1[id] = cells[num].cell_center[0]-cells[num].dx/2.0 + dis_x(gen)*cells[num].dx;
					if(x1[id] > (*box).L1[0] && x2[id] < (*box).L2[0])
						printf("It's out from right left!!!\n");
					id++;
				}
			}
		}
	}
}
else{
	if((*gas).model == "MD" && (*box).problem == vacuum){
            for(i=0; i<(*box).N_grid[0]; i++){
                  for(j=0; j<(*box).N_grid[1]; j++){
                        for(k=0; k<(*box).N_grid[2]; k++){
                              if(count > (*gas).N){
                                    i = (*box).N_grid[0];
                                    j = (*box).N_grid[1];
                                    k = (*box).N_grid[2];
                                    break;
                              }
			      else{
	                              x1[count] = (i+0.5)*(*box).Len[0]/(1.0*(*box).N_grid[0]);
	                              x2[count] = (*box).Len[1]*3.0/8.0+(j+0.5)*(*box).Len[1]*0.25/(1.0*(*box).N_grid[1]);
				      if( (*box).direction[2] == 1)
		                              x3[count] = (k+0.5)*(*box).Len[2]/(1.0*(*box).N_grid[2]);
	                              count++;
				}
                        }
                  }
            }
	}
	else if((*gas).model == "MD" && (*box).problem == inverted){
	    ///  The Nc part
            for(i=0; i<(*box).N_grid[0]; i++){
                  for(j=0; j<(*box).N_grid[1]; j++){
                        for(k=0; k<(*box).N_grid[2]; k++){
                              if(count > (*gas).Nc){
                                    i = (*box).N_grid[0];
                                    j = (*box).N_grid[1];
                                    k = (*box).N_grid[2];
                                    break;
                              }
			      else{
	                              x1[count] = (i+0.5)*(*box).Len[0]/(1.0*(*box).N_grid[0]);
	                              x2[count] = (*gas).Lv1+(j+0.5)*(*gas).Lc/(1.0*(*box).N_grid[1]);
	                              x3[count] = (k+0.5)*(*box).Len[2]/(1.0*(*box).N_grid[2]);
	                              count++;
				}
                        }
                  }
            }
	    ///  The Nh part
            for(i=0; i<(*box).N_grid[0]; i++){
                  for(j=0; j<(*box).N_grid[1]; j++){
                        for(k=0; k<(*box).N_grid[2]; k++){
                              if(count > (*gas).Nc+(*gas).Nh){
                                    i = (*box).N_grid[0];
                                    j = (*box).N_grid[1];
                                    k = (*box).N_grid[2];
                                    break;
                              }
			      else{
	                              x1[count] = (i+0.5)*(*box).Len[0]/(1.0*(*box).N_grid[0]);
	                              x2[count] = (*gas).Lv1+(*gas).Lc+(*gas).Lv+(j+0.5)*(*gas).Lh/(1.0*(*box).N_grid[1]);
	                              x3[count] = (k+0.5)*(*box).Len[2]/(1.0*(*box).N_grid[2]);
	                              count++;
				}
                        }
                  }
            }
	    /// The vapour parts
	    // Dynamically allocate a 3D array arr[height][width][depth],
	 /*  Note the parenthesis at end of new. These cause the allocated memory's
    		value to be set to zero a la calloc (value-initialize). */
	    int done;
	    int height = (*gas).N13;
	    int width  = (*gas).LvN;
	    int depth  = (*gas).N13;
    	    int ***arr = new int **[height]();
    	    for (i = 0; i < height; i++)
    	    {
        	arr[i] = new int *[width]();
        	for (j = 0; j < width; j++)
            		arr[i][j] = new int [depth]();
    	    }
 	    /// The Nv1
	    for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			for(k=0; k<depth; k++)
				arr[i][j][k] = 0;
	    for(id=0; id<(*gas).Nv1; id++){
		done = 0;
		while(done == 0){
			x1[count] =  dis_x(gen)*(*box).Len[0];
			x2[count] =  0.5*dis_x(gen)*(*gas).Lc;
			x3[count] =  dis_x(gen)*(*box).Len[2];

			i = floor(x1[count]/(*gas).eqdist);
			j = floor(x2[count]/(*gas).eqdist);
			k = floor(x3[count]/(*gas).eqdist);
			if( arr[i][j][k] == 0 ){
				done = 1;
				arr[i][j][k] = 1;
				count ++;
			}
		}
	    }
	    // deallocate arr
 	    for (i = 0; i < height; i++)
	    {
	    	for (j = 0; j < width; j++)
        		delete[] arr[i][j];
    	    	delete[] arr[i];
	    }
	    delete[] arr;
 	    /// The Nv2
	    width  = (*gas).LvN;
    	    arr = new int **[height]();
    	    for (i = 0; i < height; i++)
    	    {
        	arr[i] = new int *[width]();
        	for (j = 0; j < width; j++)
            		arr[i][j] = new int [depth]();
    	    }
	    for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			for(k=0; k<depth; k++)
				arr[i][j][k] = 0;
	    double d = pow(1.0/(*gas).nv2, 1.0/3.0);
	    for(id=0; id<(*gas).Nv2; id++){
		done = 0;
		while(done == 0){
			x1[count] =  dis_x(gen)*(*box).Len[0];
			x2[count] =  (*gas).Lv1+(*gas).Lc+dis_x(gen)*(*gas).Lv;
			x3[count] =  dis_x(gen)*(*box).Len[2];

			i = floor(x1[count]/(*gas).eqdist);
			j = floor((x2[count]-((*gas).Lv1+(*gas).Lc))/d);
			k = floor(x3[count]/(*gas).eqdist);
			if( arr[i][j][k] == 0 ){
				done = 1;
				arr[i][j][k] = 1;
				count ++;
			}
		}
	    }
 	    /// The Nv3
	    width  = floor(0.5*(*gas).LhN);
    	    arr = new int **[height]();
    	    for (i = 0; i < height; i++)
    	    {
        	arr[i] = new int *[width]();
        	for (j = 0; j < width; j++)
            		arr[i][j] = new int [depth]();
    	    }
	    for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			for(k=0; k<depth; k++)
				arr[i][j][k] = 0;
	    for(id=0; id<(*gas).Nv3; id++){
		done = 0;
		while(done == 0){
			x1[count] =  dis_x(gen)*(*box).Len[0];
			x2[count] =  (*gas).Lv1+(*gas).Lc+(*gas).Lv+(*gas).Lh+0.5*dis_x(gen)*(*gas).Lv3;
			x3[count] =  dis_x(gen)*(*box).Len[2];

			i = floor(x1[count]/(*gas).eqdist);
			j = floor((x2[count]-((*gas).Lv1+(*gas).Lc+(*gas).Lv+(*gas).Lh) )/(*gas).eqdist);
			k = floor(x3[count]/(*gas).eqdist);
			if( arr[i][j][k] == 0 ){
				done = 1;
				arr[i][j][k] = 1;
				count ++;
			}
		}
	    }
	    // deallocate arr
 	    for (i = 0; i < height; i++)
	    {
	    	for (j = 0; j < width; j++)
        		delete[] arr[i][j];
    	    	delete[] arr[i];
	    }
	    delete[] arr;
	}
	else{
  		if((*box).problem == evaporation || (*box).problem == wall){
			for ( i = 0; i < (*gas).N; ++i){


				int sg1;
				if((*gas).LvN1>0)
					sg1 = 1;
				else
					sg1 = 0;


	   			if( (*box).direction[0] == 1)
	           			x1[i] =  dis_x(gen)*(*box).Len[0];
				else
					x1[i] =  0.5*(*box).Len[0];

				if(i<(*gas).Nv1)
					x2[i] = dis_x(gen)*(*gas).Lv1;

				else if(i<(*gas).Nv1 + (*gas).Nc1)
		   			x2[i] = (*gas).Lv1 + dis_x(gen)*( ((*gas).Lc-(*gas).Lc_ratio)/2.0 );
				else if(i<(*gas).Nv1 + (*gas).Nc1+(*gas).Nc2)
		   			x2[i] = (*gas).Lv1 + sg1*((*gas).Lc-(*gas).Lc_ratio)/2.0 +dis_x(gen)*( (*gas).Lc_ratio + ((*gas).Lc-(*gas).Lc_ratio)/2.0*(1-sg1) );
				else if(i<(*gas).Nv1 + (*gas).Nc1+(*gas).Nc2+(*gas).Nc3)
		   			x2[i] = (*gas).Lv1 + sg1*((*gas).Lc-(*gas).Lc_ratio)/2.0+ (*gas).Lc_ratio + ((*gas).Lc-(*gas).Lc_ratio)/2.0*(1-sg1) +dis_x(gen)*( ((*gas).Lc-(*gas).Lc_ratio)/2.0 );

				else if(i<(*gas).Nv1 + (*gas).Nc + (*gas).Nv2)
					x2[i] = (*gas).Lv1 + (*gas).Lc + dis_x(gen)*(*gas).Lv;


	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
				else
					x3[i] =  0.5*(*box).Len[2];
        		}
/*
			int N1 = ceil( (*gas).N/(1+(*gas).n_ratio) );
			int N2 = (*gas).N - N1;
        		for ( i = 0; i < N1; ++i){
           			x1[i] =  dis_x(gen)*(*box).Len[0];
	   			x2[i] =  dis_x(gen)*(*box).Len[1]*0.5;
	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
        		}
			for ( i = N1; i < (*gas).N; ++i){
           			x1[i] =  dis_x(gen)*(*box).Len[0];
	   			x2[i] =  (*box).Len[1]*0.5+dis_x(gen)*(*box).Len[1]*0.5;
	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
        		}
*/
  		}
  		else if((*box).problem == vacuum){
        		for ( i = 0; i < (*gas).N; ++i){
	   			if( (*box).direction[0] == 1)
	           			x1[i] =  dis_x(gen)*(*box).Len[0];
				else
					x1[i] =  0.5*(*box).Len[0];
				if( i<floor((*gas).Nv1/2.0) )
					x2[i] =  dis_x(gen)*(*box).Len[1]*3.0/8.0;
				else if(i<floor((*gas).Nv1) )
					x2[i] =  (*box).Len[1]*5.0/8.0 + dis_x(gen)*(*box).Len[1]*3.0/8.0;
				else
		   			x2[i] =  (*box).Len[1]*3.0/8.0 + dis_x(gen)*(*box).Len[1]*0.25;
	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
				else
					x3[i] =  0.5*(*box).Len[2];
        		}
  		}
  		else if((*box).problem == inverted){
        		for ( i = 0; i < (*gas).N; ++i){


				int sg1, sg3;
				if((*gas).LvN1>0)
					sg1 = 1;
				else
					sg1 = 0;

				if((*gas).LvN3>0)
					sg3 = 1;
				else
					sg3 = 0;


	   			if( (*box).direction[0] == 1)
	           			x1[i] =  dis_x(gen)*(*box).Len[0];
				else
					x1[i] =  0.5*(*box).Len[0];

				if(i<(*gas).Nv1)
					x2[i] = dis_x(gen)*(*gas).Lv1;

				else if(i<(*gas).Nv1 + (*gas).Nc1)
		   			x2[i] = (*gas).Lv1 + dis_x(gen)*( ((*gas).Lc-(*gas).Lc_ratio)/2.0 );
				else if(i<(*gas).Nv1 + (*gas).Nc1+(*gas).Nc2)
		   			x2[i] = (*gas).Lv1 + sg1*((*gas).Lc-(*gas).Lc_ratio)/2.0 +dis_x(gen)*( (*gas).Lc_ratio + ((*gas).Lc-(*gas).Lh_ratio)/2.0*(1-sg1) );
				else if(i<(*gas).Nv1 + (*gas).Nc1+(*gas).Nc2+(*gas).Nc3)
		   			x2[i] = (*gas).Lv1 + sg1*((*gas).Lc-(*gas).Lc_ratio)/2.0+ (*gas).Lc_ratio + ((*gas).Lc-(*gas).Lh_ratio)/2.0*(1-sg1) +dis_x(gen)*( ((*gas).Lc-(*gas).Lc_ratio)/2.0 );

				else if(i<(*gas).Nv1 + (*gas).Nc + (*gas).Nv2)
					x2[i] = (*gas).Lv1 + (*gas).Lc + dis_x(gen)*(*gas).Lv;

				else if(i<(*gas).Nv1 + (*gas).Nc + (*gas).Nv2 + (*gas).Nh1)
					x2[i] = (*gas).Lv1 + (*gas).Lc + (*gas).Lv +  dis_x(gen)*( ((*gas).Lh-(*gas).Lh_ratio)/2.0 );
				else if(i<(*gas).Nv1 + (*gas).Nc + (*gas).Nv2 + (*gas).Nh1 + (*gas).Nh2)
					x2[i] = (*gas).Lv1 + (*gas).Lc + (*gas).Lv + ((*gas).Lh-(*gas).Lh_ratio)/2.0 +dis_x(gen)*( (*gas).Lh_ratio + ((*gas).Lh-(*gas).Lh_ratio)/2.0*(1-sg3) );
				else if(i<(*gas).Nv1 + (*gas).Nc + (*gas).Nv2 + (*gas).Nh1 + (*gas).Nh2 + (*gas).Nh3)
					x2[i] = (*gas).Lv1 + (*gas).Lc + (*gas).Lv + ((*gas).Lh-(*gas).Lh_ratio)/2.0+ (*gas).Lh_ratio +dis_x(gen)*( ((*gas).Lh-(*gas).Lh_ratio)/2.0 );

				else
					x2[i] = (*gas).Lv1 + (*gas).Lc + (*gas).Lv + (*gas).Lh + dis_x(gen)*(*gas).Lv3;

	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
				else
					x3[i] =  0.5*(*box).Len[2];
        		}
  		}
		else if((*box).problem == shock){
			for ( i = 0; i < (*gas).N; ++i){

					if( (*box).direction[0] == 1)
	           x1[i] =  dis_x(gen)*(*box).Len[0];
					else
						x1[i] =  0.5*(*box).Len[0];

					if((*gas).model == "SPH" || (*gas).model == "Hybrid"){
						if(i<(*gas).Nv1){
							x2[i] = (*gas).Lv1/(*gas).Nv1*i;
						}
						else{
							x2[i] = (*gas).Lv1 + (*gas).Lv3/(*gas).Nv3*(i-(*gas).Nv1);
						}
					}

					if((*gas).model == "DSMC_VHS"){
						if(i<(*gas).Nv1){
							x2[i] = dis_x(gen)*(*gas).Lv1;
						}
						else{
							x2[i] = (*gas).Lv1 + dis_x(gen)*(*gas).Lv3;
						}
					}
			}
		}
		else if((*box).problem == injection){
			double r, theta;
				for ( i = 0; i < (*gas).N; ++i){
					if(i<(*gas).N-(*gas).Nh){
									x1[i] =  dis_x(gen)*(*box).Len[0];
									x2[i] =  dis_x(gen)*(*box).Len[1];
					}
					else{
									//printf("here");
									r = ((*box).drop_dia/2.0) * sqrt(dis_x(gen));
									theta = dis_x(gen) * 2.0 * acos(-1.0);
									x1[i] = 0.5*(*box).Len[0] + r * cos(theta);
									x2[i] = 0.5*(*box).Len[1] + r * sin(theta);
					}
					x3[i] =  0.5*(*box).Len[2];
				}
		}
 	 	else{
    		for ( i = 0; i < (*gas).N; ++i){
	   			if( (*box).direction[0] == 1)
	           			x1[i] =  dis_x(gen)*(*box).Len[0];
					else
									x1[i] =  0.5*(*box).Len[0];
	   	 		if( (*box).direction[1] == 1)
	           			x2[i] =  dis_x(gen)*(*box).Len[1];
					else
						x2[i] =  0.5*(*box).Len[1];
	   			if( (*box).direction[2] == 1)
	           			x3[i] =  dis_x(gen)*(*box).Len[2];
				else
					x3[i] =  0.5*(*box).Len[2];
    			}
	 	}
 	}
}

if((*box).problem == dcones || (*box).problem == flatnose){
  double std = sqrt((*gas).kb*(*gas).T/(*gas).m);
  std::normal_distribution<> dis_u((*gas).U0[0],std);
  std::normal_distribution<> dis_v((*gas).U0[1],std);
  std::normal_distribution<> dis_w((*gas).U0[2],std);
    for ( i = 0; i < (*gas).N; ++i){
        U1[i] =  dis_u(gen);
	U2[i] =  dis_v(gen);
	U3[i] =  dis_w(gen);}
}
else if((*box).problem == relaxation){
  double T1 = 300.0, T2 = 300;
  double U1_1 = 100.0, U2_1 = -100.0;
  double U1_2 = 50.0, U2_2 = -50.0;
  double U1_3 = 200.0, U2_3 = -200.0;
 // double std1 = sqrt((*gas).kb*T1/(*gas).m), std2 = sqrt((*gas).kb*T2/(*gas).m);
  double std1 = 100.0, std2 = 300.0;
  std::normal_distribution<> dis_u1(U1_1,std1);
  std::normal_distribution<> dis_v1(U1_2,std1);
  std::normal_distribution<> dis_w1(U1_3,std1);

  std::normal_distribution<> dis_u2(U2_1,std2);
  std::normal_distribution<> dis_v2(U2_2,std2);
  std::normal_distribution<> dis_w2(U2_3,std2);

   for ( i = 0; i < (*gas).N; ++i){
	if( i< (*gas).N/2){
		U1[i] =  dis_u1(gen);
		U2[i] =  dis_v1(gen);
		U3[i] =  dis_w1(gen);
	}
	else{
		U1[i] =  dis_u2(gen);
		U2[i] =  dis_v2(gen);
		U3[i] =  dis_w2(gen);
	}
  }
}
else if((*box).problem == shock){
	double stdT;
	std::normal_distribution<> dis_norm(0.0,1.0);
	for ( i = 0; i < (*gas).N ; ++i){
		if( (*gas).model == "SPH" || (*gas).model == "Hybrid"){
			U1[i] = 0.0;
			U2[i] = 0.0;
			U3[i] = 0.0;
		}
		else if( (*gas).model == "DSMC_VHS" ){
			if(i<(*gas).Nv1){
				stdT = sqrt( (*gas).kb*(*gas).Tv1/(*gas).m );
				U1[i] = dis_norm(gen)*stdT;
				U2[i] = dis_norm(gen)*stdT;
				U3[i] = dis_norm(gen)*stdT;
			}
			else{
				stdT = sqrt( (*gas).kb*(*gas).Tv3/(*gas).m );
				U1[i] = dis_norm(gen)*stdT;
				U2[i] = dis_norm(gen)*stdT;
				U3[i] = dis_norm(gen)*stdT;
			}
		}
	}
}
else{
  if((*box).problem == inverted){
	double std = sqrt((*gas).kb*(*gas).Tv1/(*gas).m);
  	std::normal_distribution<> dis_u1((*gas).U0[0],std);
  	std::normal_distribution<> dis_v1((*gas).U0[1],std);
  	std::normal_distribution<> dis_w1((*gas).U0[2],std);
  	for ( i = 0; i < (*gas).Nv1 ; ++i){
        	U1[i] =  dis_u1(gen);
		U2[i] =  dis_v1(gen);
		U3[i] =  dis_w1(gen);}

	std = sqrt((*gas).kb*(*gas).Tc/(*gas).m);
  	std::normal_distribution<> dis_uc((*gas).U0[0],std);
  	std::normal_distribution<> dis_vc((*gas).U0[1],std);
  	std::normal_distribution<> dis_wc((*gas).U0[2],std);
  	for ( i = (*gas).Nv1; i < (*gas).Nv1 + (*gas).Nc; ++i){
        	U1[i] =  dis_uc(gen);
		U2[i] =  dis_vc(gen);
		U3[i] =  dis_wc(gen);}


	std = sqrt((*gas).kb*((*gas).Tv2)/(*gas).m);
  	std::normal_distribution<> dis_u2((*gas).U0[0],std);
  	std::normal_distribution<> dis_v2((*gas).U0[1],std);
  	std::normal_distribution<> dis_w2((*gas).U0[2],std);
  	for ( i = (*gas).Nv1 + (*gas).Nc; i < (*gas).Nv1 + (*gas).Nc + (*gas).Nv2; ++i){
        	U1[i] =  dis_u2(gen);
		U2[i] =  dis_v2(gen);
		U3[i] =  dis_w2(gen);}


	std = sqrt((*gas).kb*(*gas).Th/(*gas).m);
  	std::normal_distribution<> dis_uh((*gas).U0[0],std);
  	std::normal_distribution<> dis_vh((*gas).U0[1],std);
  	std::normal_distribution<> dis_wh((*gas).U0[2],std);
  	for ( i = (*gas).Nv1 + (*gas).Nc+ (*gas).Nv2; i < (*gas).Nv1 + (*gas).Nc+ (*gas).Nv2 + (*gas).Nh ; ++i){
        	U1[i] =  dis_uh(gen);
		U2[i] =  dis_vh(gen);
		U3[i] =  dis_wh(gen);}

	std = sqrt((*gas).kb*(*gas).Tv3/(*gas).m);
  	std::normal_distribution<> dis_u3((*gas).U0[0],std);
  	std::normal_distribution<> dis_v3((*gas).U0[1],std);
  	std::normal_distribution<> dis_w3((*gas).U0[2],std);
  	for ( i = (*gas).Nv1 + (*gas).Nc+ (*gas).Nv2 + (*gas).Nh; i < (*gas).N ; ++i){
        	U1[i] =  dis_u3(gen);
		U2[i] =  dis_v3(gen);
		U3[i] =  dis_w3(gen);}


  }
  else{
	double T =  max( (*box).T_wall_1+(*box).T_wall_2 , (*box).T_wall_3+(*box).T_wall_4 ) /2.0;
  	T = max((*gas).T,T);
  	double std = sqrt((*gas).kb*T/(*gas).m);
  	std::normal_distribution<> dis_u((*gas).U0[0],std);
  	std::normal_distribution<> dis_v((*gas).U0[1],std);
  	std::normal_distribution<> dis_w((*gas).U0[2],std);
  	for ( i = 0; i < (*gas).N; ++i){
        	U1[i] =  dis_u(gen);
		U2[i] =  dis_v(gen);
		U3[i] =  dis_w(gen);}


   	for ( i = 0; i < (*gas).N; ++i){
		if( (*box).problem == evaporation || (*box).problem == wall){
			j = floor(x1[i]/(*box).delta_dim[0]);
    			k = floor(x2[i]/(*box).delta_dim[1]);
    			num =  k*(*box).N[0] + j;
		}
		else{
		}
  	}
  }
}
}

void BC_wall(double *U1,double *U2,double *U3, double *x2, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x2_old){

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
double xw;
double yw;
double zw;
int num_hitting_wall = 0;
// double Tol = 0.0;
// double Tol2 = 0.0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j;
double U_eff[3];
int done;
int dummy;
// #pragma omp parallel for  private(dt_remain, num_hitting_wall, num, temp_U20, energy,  U_dummy)
for(i=0; i<(*gas).N; i++){
  dt = (*gas).delta_t;
  done = 0;
  dummy = 0;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 5){
      printf("*** TOO MANY BC OPERATIONS ***\n");
      printf(" point %d \n", i);
      printf(" with coord  %e,  with dt = %e\n",  x2[i] , dt);
      if(dummy > 20){
	exit(1);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    if(dt < 1e-15){
          if(x2[i] < epsilon)
              x2[i] = epsilon;
          else if(x2[i] > (*box).Len[1] - epsilon)
              x2[i] = (*box).Len[1] - epsilon;
   }
    //////////////////
    //////////////////              UPPER WALL
    //////////////////
    else if(x2[i] > (*box).Len[1]) {
      // reflect
			U2[i] = - U2[i];
			x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] );
    }
    //////////////////
    //////////////////              LOWER WALL
    //////////////////
    else if(x2[i] < 0.0){
      yw = epsilon;
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[1];
      else
            temp_U0 = U2[i];


      U1[i] = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal(gen) + (1.0*(*box).U_wall_1);
      U2[i] =   sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) * Normal(gen);

     if(dt >1e-15){
           dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
           dt_remain = min(dt, dt_remain);
           x2[i] =  yw + dt_remain*U2[i];
           dt = dt - dt_remain;
           x2_old[i] = yw;
      }
      else{
            x2[i] = yw;
      }
    }
    if(x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0)
      done = 1;
    dummy ++;
  }
}
}


void BC_thermal_wall_3D(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, double *x3_old){

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
double xw;
double yw;
double zw;
int num_hitting_wall = 0;
// double Tol = 0.0;
// double Tol2 = 0.0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j;
double U_eff[3];
int done;
int dummy;
// #pragma omp parallel for  private(dt_remain, num_hitting_wall, num, temp_U20, energy,  U_dummy)
for(i=0; i<(*gas).N; i++){
  dt = (*gas).delta_t;
  done = 0;
  dummy = 0;
  U_eff[0] = fabs(x1[i]-x1_old[i])/(*gas).delta_t;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  U_eff[2] = fabs(x3[i]-x3_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 5){
      printf("*** TOO MANY BC OPERATIONS ***\n");
      printf(" point %d \n", i);
      printf(" with coord  %e,  %e, %e with dt = %e\n", x1[i], x2[i], x3[i], dt);
      if(dummy > 20){
	exit(1);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    if(dt < 1e-15){
          if(x1[i] < epsilon)
            x1[i] = epsilon;
          else if(x1[i] > (*box).Len[0] - epsilon)
            x1[i] = (*box).Len[0] - epsilon;
          if(x2[i] < epsilon)
              x2[i] = epsilon;
          else if(x2[i] > (*box).Len[1] - epsilon)
              x2[i] = (*box).Len[1] - epsilon;
          if(x3[i] < epsilon)
                x3[i] = epsilon;
          else if(x3[i] > (*box).Len[2] - epsilon)
                x3[i] = (*box).Len[2] - epsilon;
   }
    /*
    if(x1[i] - (*box).Len[0] > 0.0 &&  x2[i] - (*box).Len[1] > 0.0 && dt < 1e-15){
          x1[i] = (*box).Len[0] - epsilon;
          x2[i] = (*box).Len[1] - epsilon;
    }
    else if(x1[i] - (*box).Len[0] > 0.0 &&  x2[i]  < 0.0 && dt < 1e-15){
          x1[i] = (*box).Len[0] - epsilon;
          x2[i] =  epsilon;
    }
    else if(x1[i] < 0.0 &&  x2[i] - (*box).Len[1] > 0.0 && dt < 1e-15){
          x1[i] =  epsilon;
          x2[i] = (*box).Len[1] - epsilon;
    }
    else if(x1[i] < 0.0 &&  x2[i] < 0.0 && dt < 1e-15){
          x1[i] = epsilon;
          x2[i] = epsilon;
    }
    */
    //////////////////
    //////////////////              RIGHT WALL
    //////////////////
     if(x1[i] - (*box).Len[0] > 0.0){
       xw = (*box).Len[0]-epsilon;
       if(fabs(x1_old[i]-x1[i])>1e-15)
            intersect_3D_x(xw, &yw, &zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
      else{
            yw = 0.5*(x2_old[i]+x2[i]);
            zw = 0.5*(x3_old[i]+x3[i]);
      }
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];

      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[0];
      else
            temp_U0 = U1[i];

      U1[i] =  - sqrt(2.0* (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen )+ (1.0*(*box).U_wall_4);
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );

      if(dt >1e-15 ){
	           dt_remain = fabs ( (x1[i]-xw)/temp_U0 );
                 dt_remain = min(dt, dt_remain);
	           x1[i] =  xw + dt_remain*U1[i];
                 x2[i] =  yw + dt_remain*U2[i];
                 x3[i] =  zw + dt_remain*U3[i];
                 dt = dt - dt_remain;
                 x1_old[i] = xw;
                 x2_old[i] = yw;
                 x3_old[i] = zw;
      }
      else{
      	x1[i] = xw;
      }
      M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      flag[i] += 1;
    }
    //////////////////
    //////////////////              LEFT WALL
    //////////////////
    else if(x1[i] < 0.0){
         xw = epsilon;
         if(fabs(x1_old[i]-x1[i])>1e-15)
             intersect_3D_x(xw, &yw, &zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
         else{
                  yw = 0.5*(x2_old[i]+x2[i]);
                  zw = 0.5*(x3_old[i]+x3[i]);
          }
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[0];
      else
            temp_U0 = U1[i];
      U1[i] = sqrt(2.0* (*gas).kb*(*box).T_wall_3/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_3/ (*gas).m ) * Normal(gen) + (1.0*(*box).U_wall_3);
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_3/ (*gas).m ) * Normal(gen);

     if(dt >1e-15 ){
           dt_remain = fabs ( (x1[i]-xw)/temp_U0 );
           dt_remain = min(dt, dt_remain);
           x1[i] =  xw + dt_remain*U1[i];
           x2[i] =  yw + dt_remain*U2[i];
           x3[i] =  zw + dt_remain*U3[i];
           dt = dt - dt_remain;
           x1_old[i] = xw;
           x2_old[i] = yw;
           x3_old[i] = zw;
      }
      else{
            x1[i] = xw;
      }
      M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      flag[i] += 1000;
    }
    //////////////////
    //////////////////              UPPER WALL
    //////////////////
    else if(x2[i] > (*box).Len[1] -(*gas).sigma-epsilon) {
      yw = (*box).Len[1]-(*gas).sigma-epsilon;
      intersect_3D_y(&xw, yw, &zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[1];
      else
            temp_U0 = U2[i];

      U1[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen ) + (1.0*(*box).U_wall_2);
      U2[i] =  - sqrt(2.0* (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen );

      if(dt >1e-15){
            dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
            dt_remain = min(dt, dt_remain);
            x1[i] =  xw + dt_remain*U1[i];
            x2[i] =  yw + dt_remain*U2[i];
            x3[i] =  zw + dt_remain*U3[i];
            dt = dt - dt_remain;
            x1_old[i] = xw;
            x2_old[i] = yw;
            x3_old[i] = zw;
      }
      else{
            x2[i] = yw;
      }
      M[2] = M[2] + (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1;
    }
    //////////////////
    //////////////////              LOWER WALL
    //////////////////
    else if(x2[i] < (*gas).sigma+epsilon){
      yw = (*gas).sigma + epsilon;
      intersect_3D_y(&xw, yw, &zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[1];
      else
            temp_U0 = U2[i];


      U1[i] = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal(gen) + (1.0*(*box).U_wall_1);
      U2[i] =   sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) * Normal(gen);

     if(dt >1e-15){
           dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
           dt_remain = min(dt, dt_remain);
           x1[i] =  xw + dt_remain*U1[i];
           x2[i] =  yw + dt_remain*U2[i];
           x3[i] =  zw + dt_remain*U3[i];
           dt = dt - dt_remain;
           x1_old[i] = xw;
           x2_old[i] = yw;
           x3_old[i] = zw;
      }
      else{
            x2[i] = yw;
      }
      M[3] = M[3] +  (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1000;
    }
    //////////////////
    //////////////////              BEHIND WALL
    //////////////////
    else if(x3[i] - (*box).Len[2] > 0.0) {
      zw = (*box).Len[2]-epsilon;
      intersect_3D_z(&xw, &yw, zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[2];
      else
            temp_U0 = U3[i];

      U1[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen );
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen );
      U3[i] =  - sqrt(2.0* (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * sqrt(-log(uniform(gen2)));


      if(dt >1e-15){
            dt_remain = fabs ( (x3[i]-zw)/temp_U0 );
            dt_remain = min(dt, dt_remain);
            x1[i] =  xw + dt_remain*U1[i];
            x2[i] =  yw + dt_remain*U2[i];
            x3[i] =  zw + dt_remain*U3[i];
            dt = dt - dt_remain;
            x1_old[i] = xw;
            x2_old[i] = yw;
            x3_old[i] = zw;
      }
      else{
            x3[i] = zw;
      }
      M[4] = M[4] + (*gas).Fn*(*gas).m*( fabs(U3[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1;
    }
    //////////////////
    //////////////////              FRONT WALL
    //////////////////
    else if(x3[i] < 0.0){
      zw = epsilon;
      intersect_3D_z(&xw, &yw, zw, x1_old[i], x2_old[i], x3_old[i], x1[i], x2[i], x3[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[2];
      else
            temp_U0 = U3[i];


      U1[i] = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal(gen);
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) * Normal(gen);
      U3[i] =   sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gen2)));

     if(dt >1e-15){
           dt_remain = fabs ( (x3[i]-zw)/temp_U0 );
           dt_remain = min(dt, dt_remain);
           x1[i] =  xw + dt_remain*U1[i];
           x2[i] =  yw + dt_remain*U2[i];
           x3[i] =  zw + dt_remain*U3[i];
           dt = dt - dt_remain;
           x1_old[i] = xw;
           x2_old[i] = yw;
           x3_old[i] = zw;
      }
      else{
            x3[i] = zw;
      }
      M[5] = M[5] +  (*gas).Fn*(*gas).m*( fabs(U3[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1000;
    }
    if(x1[i] - (*box).Len[0] < 0.0 && x1[i] > 0.0 && x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0 & x3[i] - (*box).Len[2] < 0.0 && x3[i] > 0.0)
      done = 1;
    dummy ++;
  }
}

  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
}


void BC_inverted(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, double *x1_old, double *x2_old, double *x3_old){
//   printf("BC Starts\n");

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::uniform_real_distribution<> uniform2(-1.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
  double std = sqrt((*gas).kb*(*box).T_wall_1/(*gas).m);

  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
double xw;
double yw;
int num_hitting_wall = 0;
// double Tol = 0.0;
// double Tol2 = 0.0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j,k;
int done;
int dummy;
double U_eff[2];

double Tol1=0.0,Tol2=0.0;
  // checking 1st dimension
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  dt = (*gas).delta_t;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 20){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %f \n", x2[i]/(*box).Len[1] );
      if(dummy > 30){
	exit(1);
      }
//      scanf("%d", &dummy);
    }
    if(	(*box).direction[0]==1 ){
	    while(x1[i] - (*box).Len[0] > Tol2){
      		x1[i] = x1[i] - (*box).Len[0];
      		M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    	    while(x1[i] < Tol1){
      		x1[i] = x1[i] + (*box).Len[0];
      		M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    }
    if(	(*box).direction[1]==1 ){
    	//while(x2[i] > (*box).Len[1]  - (*box).ghost  && x2_old[i] <  (*box).Len[1]  - (*box).ghost ){
	//while(x2[i] > (*box).Len[1] -epsilon && x2_old[i] <  (*box).Len[1] - epsilon ){
	while(x2[i] > (*box).Len[1] -epsilon ){

      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}
/*
		//xw = (*box).Len[1] - (*box).ghost - epsilon ;
		xw = (*box).Len[1]   ;
		x2[i] = xw - abs(x2[i] - xw) -1.0*epsilon;
		x2_old[i] = xw;
      		k = floor(x2[i]/(*box).delta_dim[1]);
      		num =  k;
      		index[i] = num;
		U2[i] = -U2[i];
*/

     		//yw = (*box).Len[1]-epsilon- (*box).ghost;
		yw = (*box).Len[1]-epsilon;

      		temp_U1 = U1[i];
      		temp_U2 = U2[i];
      		temp_U3 = U3[i];
      		if(fabs(dt-(*gas).delta_t)<1e-15)
        	    	temp_U0 = U_eff[1];
      		else
            		temp_U0 = U2[i];

      		U1[i] = sqrt( (*gas).kb*(*gas).Th/( (*gas).m) ) * Normal ( gen );// + 1000.0;
      		U2[i] =  - sqrt(2.0* (*gas).kb*(*gas).Th/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      		U3[i] = sqrt( (*gas).kb*(*gas).Th/( (*gas).m) ) * Normal ( gen );

      		if(dt >1e-15){
            		dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
            		x2[i] =  yw + dt_remain*U2[i] ;
		        dt = dt - dt_remain;
            		x2_old[i] = yw;
      		}
      		else{
	        	    x2[i] = yw-epsilon;
      		}


	        M[3] = M[3] + (*gas).Fn*fabs((*gas).m*U2[i]);
      		num_hitting_wall++;
    	}
	//while(x2[i] <    (*box).ghost  && x2_old[i] >   (*box).ghost){
//	while(x2[i] <  epsilon   && x2_old[i] > + epsilon ){
		while(x2[i] <  epsilon   ){
      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}
/*
		//xw = epsilon+ (*box).ghost ;
		xw = 0.0;
		x2[i] = xw + abs(x2[i] - xw)+1.0*epsilon;
		x2_old[i] = xw;
      		k = floor(x2[i]/(*box).delta_dim[1]);
      		num =  k;
      		index[i] = num;
		U2[i] = -U2[i];
*/

		//yw = (*box).ghost + epsilon;
		yw = epsilon;

      		temp_U1 = U1[i];
      		temp_U2 = U2[i];
      		temp_U3 = U3[i];
      		if(fabs(dt-(*gas).delta_t)<1e-15)
        	    	temp_U0 = U_eff[1];
      		else
            		temp_U0 = U2[i];

      		U1[i] = sqrt( (*gas).kb*(*gas).Tc/( (*gas).m) ) * Normal ( gen );// - 1000.0;
      		U2[i] =  + sqrt(2.0* (*gas).kb*(*gas).Tc/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      		U3[i] = sqrt( (*gas).kb*(*gas).Tc/( (*gas).m) ) * Normal ( gen );

      		if(dt >1e-15){
            		dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
            		x2[i] =  yw + dt_remain*U2[i] ;
		        dt = dt - dt_remain;
            		x2_old[i] = yw;
      		}
      		else{
	        	    x2[i] = yw+epsilon ;
      		}

      		M[3] = M[3] +  (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      		num_hitting_wall++;
    	}
    }
    if(	(*box).direction[2]==1 ){
    	while(x3[i] - (*box).Len[2] > Tol2) {
      		x3[i] = x3[i] - (*box).Len[2];
      		M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    	while(x3[i] < Tol1){
      		x3[i] = x3[i] + (*box).Len[2];
      		M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    }
    done =1;
    if(	(*box).direction[0]==1 ){
	if(x1[i] - (*box).Len[0] > 0.0 || x1[i] < 0.0)
		done = 0;
    }
    if(	(*box).direction[1]==1 ){
	//if ( ( x2[i] <    (*box).ghost && x2_old[i] >   (*box).ghost ) || (x2[i] > (*box).Len[1]  - (*box).ghost && x2_old[i] <  (*box).Len[1]  - (*box).ghost))
	if ( ( x2[i] <    0.0  ) || (x2[i] > (*box).Len[1]  ) )
			done = 0;
    }
    if(	(*box).direction[2]==1 ){
	if(x3[i] - (*box).Len[2] > 0.0 || x3[i] < 0.0)
		done = 0;
    }
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
}



void BC_evaporation(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, double *x3_old){
//   printf("BC Starts\n");

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::uniform_real_distribution<> uniform2(-1.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
  double std = sqrt((*gas).kb*(*box).T_wall_1/(*gas).m);


//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
double xw;
double yw;
int num_hitting_wall = 0;
// double Tol = 0.0;
// double Tol2 = 0.0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j,k;
double U_eff[2];
int done;
int dummy;
double Tol1=0.0,Tol2=0.0;
  // checking 1st dimension
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  dt = (*gas).delta_t;
  U_eff[0] = fabs(x1[i]-x1_old[i])/(*gas).delta_t;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 20){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %f,  %f, %f \n", x1[i]/(*box).Len[0], x2[i]/(*box).Len[1], x3[i]/(*box).Len[2]);
      if(dummy > 30){
	exit(1);
      }
//      scanf("%d", &dummy);
    }
    if(	(*box).direction[0]==1 ){
	    while(x1[i] - (*box).Len[0] > Tol2){
      		x1[i] = x1[i] - (*box).Len[0];
      		M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    	    while(x1[i] < Tol1){
      		x1[i] = x1[i] + (*box).Len[0];
      		M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    }
    if(	(*box).direction[1]==1 ){
    	while(x2[i] > (*box).Len[1] -epsilon){
      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}
/*
      		yw = (*box).Len[1]-epsilon;
      		xw = line_interesect_x(yw, x1_old[i], x2_old[i], x1[i], x2[i]);
      		num = index[i];

      		temp_U1 = U1[i];
      		temp_U2 = U2[i];
      		temp_U3 = U3[i];
      		if(fabs(dt-(*gas).delta_t)<1e-15)
            		temp_U0 = U_eff[1];
      		else
            		temp_U0 = U2[i];

      		M[2] = M[2] + (*gas).Fn*fabs((*gas).m*U2[i]);
*/
		//x2[i] = (*box).Len[1]*3.0/8.0+uniform(gen)*(*box).Len[1]/4.0;
		x2[i] = (*box).Len[1]/2.0+uniform2(gen)*(*box).Len[1]/4.0;
		//x2[i] = x2[i] - (*box).Len[1];
      		k = floor(x2[i]/(*box).delta_dim[1]);
      		num =  k;
      		index[i] = num;

      		U1[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[0];
      		U2[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[1];
      		U3[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[2];

//      		U1[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[0];
//      		U2[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[1];
//      		U3[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[2];

	        M[3] = M[3] + (*gas).Fn*fabs((*gas).m*U2[i]);
      		num_hitting_wall++;
      		flag[i] += 1;
    	}
	while(x2[i] < +epsilon){
      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}
/*
      		yw = epsilon;
      		xw = line_interesect_x(yw, x1_old[i], x2_old[i], x1[i], x2[i]);
      		num = index[i];

      		temp_U1 = U1[i];
      		temp_U2 = U2[i];
      		temp_U3 = U3[i];
      		if(fabs(dt-(*gas).delta_t)<1e-15)
            		temp_U0 = U_eff[1];
      		else
            		temp_U0 = U2[i];
*/
		//x2[i] = (*box).Len[1]*3.0/8.0+uniform(gen)*(*box).Len[1]/4.0;
		x2[i] = (*box).Len[1]/2.0+uniform2(gen)*(*box).Len[1]/4.0;
	//	x2[i] = x2[i] + (*box).Len[1];
	        k = floor(x2[i]/(*box).delta_dim[1]);
      		num =  k;
      //		index[i] = num;

      		U1[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[0];
      		U2[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[1];
      		U3[i] = sqrt( (*gas).kb*(*gas).T/(*gas).m ) * Normal(gen);// + cells[num].U_space[2];

 //     		U1[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[0];
 //     		U2[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[1];
 //     		U3[i] = sqrt( (*gas).kb*cells[num].T/(*gas).m ) * Normal(gen);// + cells[num].U_space[2];

      		M[3] = M[3] +  (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      		num_hitting_wall++;
      		flag[i] += 1000;
    	}
    }
    if(	(*box).direction[2]==1 ){
    	while(x3[i] - (*box).Len[2] > Tol2) {
      		x3[i] = x3[i] - (*box).Len[2];
      		M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    	while(x3[i] < Tol1){
      		x3[i] = x3[i] + (*box).Len[2];
      		M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    }
    done =1;
    if(	(*box).direction[0]==1 ){
	if(x1[i] - (*box).Len[0] > 0.0 || x1[i] < 0.0)
		done = 0;
    }
    if(	(*box).direction[1]==1 ){
	if(x2[i] - (*box).Len[1] > 0.0 || x2[i] < 0.0)
		done = 0;
    }
    if(	(*box).direction[2]==1 ){
	if(x3[i] - (*box).Len[2] > 0.0 || x3[i] < 0.0)
		done = 0;
    }
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
}


void BC_evaporation2(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, double *x1_old, double *x2_old, double *x3_old){
//   printf("BC Starts\n");

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::uniform_real_distribution<> uniform2(-1.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
  double std = sqrt((*gas).kb*(*box).T_wall_1/(*gas).m);

  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;
double xw;
double yw;
int num_hitting_wall = 0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j,k;
int done;
int dummy;
double U_eff[2];

double Tol1=0.0,Tol2=0.0;
  // checking 1st dimension
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  dt = (*gas).delta_t;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 20){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %f \n", x2[i]/(*box).Len[1] );
      if(dummy > 30){
	exit(1);
      }
//      scanf("%d", &dummy);
    }
    if(	(*box).direction[0]==1 ){
	    while(x1[i] - (*box).Len[0] > Tol2){
      		x1[i] = x1[i] - (*box).Len[0];
      		M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    	    while(x1[i] < Tol1){
      		x1[i] = x1[i] + (*box).Len[0];
      		M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      		num_hitting_wall++;
    	    }
    }
    if(	(*box).direction[1]==1 ){
    	while(x2[i] > (*box).Len[1] -epsilon){
      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}

		x2[i] = x2[i] - (*box).Len[1] + epsilon;

	        M[3] = M[3] + (*gas).Fn*fabs((*gas).m*U2[i]);
      		num_hitting_wall++;
    	}
	while(x2[i] < epsilon){
      		if((*box).step > (*box).after){
			(*box).out += 1.0;
      		}

		x2[i] = x2[i] + (*box).Len[1] - epsilon;

      		M[3] = M[3] +  (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      		num_hitting_wall++;
    	}
    }
    if(	(*box).direction[2]==1 ){
    	while(x3[i] - (*box).Len[2] > Tol2) {
      		x3[i] = x3[i] - (*box).Len[2];
      		M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    	while(x3[i] < Tol1){
      		x3[i] = x3[i] + (*box).Len[2];
      		M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      		num_hitting_wall++;
    	}
    }
    done =1;
    if(	(*box).direction[0]==1 ){
	if(x1[i] - (*box).Len[0] > 0.0 || x1[i] < 0.0)
		done = 0;
    }
    if(	(*box).direction[1]==1 ){
	if(x2[i] - (*box).Len[1] > 0.0 || x2[i] < 0.0)
		done = 0;
    }
    if(	(*box).direction[2]==1 ){
	if(x3[i] - (*box).Len[2] > 0.0 || x3[i] < 0.0)
		done = 0;
    }
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
}


void store_old_positions(double *x1,double *x2, double *x3, double *x1_old, double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box){
	int i;
	if( (*box).direction[0] == 1 ){
        	for(i=0; i<(*gas).N; i++){
			x1_old[i] = x1[i];}
	}
	if( (*box).direction[1] == 1 ){
		for(i=0; i<(*gas).N; i++){
			x2_old[i] = x2[i];}
	}
	if( (*box).direction[2] == 1 ){
		for(i=0; i<(*gas).N; i++){
			x3_old[i] = x3[i];}
        }

	// Update time for postproc
	if ( (*box).step > (*gas).times*(*box).after ){

		(*gas).avgeraging_time_till_now +=  (*gas).delta_t;
	}

	if ((*box).step % (*box).every == 0)
		printf("stored old positions.\n");
}

void BC_thermal_wall(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x1_old, double *x2_old){

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);
//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;

double epsilon = 1e-15;//max((*box).delta_dim[0],(*box).delta_dim[1])/100;
double xw;
double yw;
int num_hitting_wall = 0;
// double Tol = 0.0;
// double Tol2 = 0.0;
double dt_remain;
int num;
double U_dummy, temp_U0, temp_U1, temp_U2, temp_U3;
double energy;
int id_0[3], dist, id_temp1, num_temp1, id_temp2, num_temp2, id_temp, num_temp;
double dt;
  // checking 1st dimension
int i,j;
double U_eff[2];
int done;
int dummy;
// #pragma omp parallel for  private(dt_remain, num_hitting_wall, num, temp_U20, energy,  U_dummy)
for(i=0; i<(*gas).N; i++){
  dt = (*gas).delta_t;
  done = 0;
  dummy = 0;
  U_eff[0] = fabs(x1[i]-x1_old[i])/(*gas).delta_t;
  U_eff[1] = fabs(x2[i]-x2_old[i])/(*gas).delta_t;
  while (done == 0){
    dt_remain = 0.0;
    if (dummy > 10){
      printf("*** TOO MANY BC OPERATIONS ***\n");
      printf(" point %d \n", i);
      printf(" with coord  %e,  %e\n", x1[i], x2[i]);
      if(dummy > 20){
	exit(1);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    if(x1[i] - (*box).Len[0] > 0.0 &&  x2[i] - (*box).Len[1] > 0.0 && dt < 1e-15){
          x1[i] = (*box).Len[0] - epsilon;
          x2[i] = (*box).Len[1] - epsilon;
    }
    else if(x1[i] - (*box).Len[0] > 0.0 &&  x2[i]  < 0.0 && dt < 1e-15){
          x1[i] = (*box).Len[0] - epsilon;
          x2[i] =  epsilon;
    }
    else if(x1[i] < 0.0 &&  x2[i] - (*box).Len[1] > 0.0 && dt < 1e-15){
          x1[i] =  epsilon;
          x2[i] = (*box).Len[1] - epsilon;
    }
    else if(x1[i] < 0.0 &&  x2[i] < 0.0 && dt < 1e-15){
          x1[i] = epsilon;
          x2[i] = epsilon;
    }
     if(x1[i] - (*box).Len[0] > 0.0){
       xw = (*box).Len[0]-epsilon;
       if(fabs(x1_old[i]-x1[i])>1e-15)
            yw = line_interesect_y(xw, x1_old[i], x2_old[i], x1[i], x2[i]);
      else
            yw = 0.5*(x2_old[i]+x2[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];

      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[0];
      else
            temp_U0 = U1[i];

      U1[i] =  - sqrt(2.0* (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen )+ (1.0*(*box).U_wall_4);
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );

      if(dt >1e-15 ){
	           dt_remain = fabs ( (x1[i]-xw)/temp_U0 );
	           x1[i] =  xw + dt_remain*U1[i];
                 x2[i] =  yw + dt_remain*U2[i];
                 dt = dt - dt_remain;
                 x1_old[i] = xw;
                 x2_old[i] = yw;
      }
      else{
            //x1[i] = x1[i] - 2.0*fabs( x1[i] - (*box).Len[0] );
      	x1[i] = xw;
            //x2[i] = yw;//line_interesect_y(xw, x1_old[i], x2_old[i], x1[i], x2[i]);
      }
      if(x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0){
            yw = line_interesect_y(xw+2.0*epsilon, x1_old[i], x2_old[i], x1[i], x2[i]);
            update_face_info_BC(U1[i],U2[i],U3[i], x1[i], x2[i], xw+2.0*epsilon, yw, gas, box, cells);
      }
      M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      flag[i] += 1;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    else if(x1[i] < 0.0){
         xw = epsilon;
         if(fabs(x1_old[i]-x1[i])>1e-15)
             yw = line_interesect_y(xw, x1_old[i], x2_old[i], x1[i], x2[i]);
         else
             yw = 0.5*(x2_old[i]+x2[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[0];
      else
            temp_U0 = U1[i];
      U1[i] = sqrt(2.0* (*gas).kb*(*box).T_wall_3/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U2[i] = sqrt( (*gas).kb*(*box).T_wall_3/ (*gas).m ) * Normal(gen) + (1.0*(*box).U_wall_3);
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_3/ (*gas).m ) * Normal(gen);

     if(dt >1e-15 ){
           dt_remain = fabs ( (x1[i]-xw)/temp_U0 );
           x1[i] =  xw + dt_remain*U1[i];
           x2[i] =  yw + dt_remain*U2[i];
           dt = dt - dt_remain;
           x1_old[i] = xw;
           x2_old[i] = yw;
      }
      else{
            //x1[i] = x1[i] + 2.0*fabs( x1[i]  );
            x1[i] = xw;
            //x2[i] = yw;
      }
      if( x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0){
            yw = line_interesect_y(xw-2.0*epsilon, x1_old[i], x2_old[i], x1[i], x2[i]);
            update_face_info_BC(U1[i],U2[i],U3[i], x1[i], x2[i], xw-2.0*epsilon, yw, gas, box, cells);
      }
      M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      flag[i] += 1000;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    else if(x2[i] - (*box).Len[1] > 0.0) {
      yw = (*box).Len[1]-epsilon;
      xw = line_interesect_x(yw, x1_old[i], x2_old[i], x1[i], x2[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[1];
      else
            temp_U0 = U2[i];

      U1[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen ) + (1.0*(*box).U_wall_2);
      U2[i] =  - sqrt(2.0* (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gen );

      if(dt >1e-15){
            dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
            x1[i] =  xw + dt_remain*U1[i];
            x2[i] =  yw + dt_remain*U2[i];
            //x2[i] =  yw + dt_remain*U2[i];
            dt = dt - dt_remain;
            x1_old[i] = xw;
            x2_old[i] = yw;
      }
      else{
            //x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] );
            x2[i] = yw;
            //x1[i] = xw;
      }
      if(x1[i] - (*box).Len[0] < 0.0 && x1[i] > 0.0){
            xw = line_interesect_x(yw+2.0*epsilon, x1_old[i], x2_old[i], x1[i], x2[i]);
            update_face_info_BC(U1[i],U2[i],U3[i], x1[i], x2[i], xw, yw+2.0*epsilon, gas, box, cells);
      }
      M[2] = M[2] + (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1;
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    else if(x2[i] < 0.0){
      yw = epsilon;
      xw = line_interesect_x(yw, x1_old[i], x2_old[i], x1[i], x2[i]);
      num = index[i];

      temp_U1 = U1[i];
      temp_U2 = U2[i];
      temp_U3 = U3[i];
      if(fabs(dt-(*gas).delta_t)<1e-15)
            temp_U0 = U_eff[1];
      else
            temp_U0 = U2[i];


      U1[i] = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal(gen) + (1.0*(*box).U_wall_1);
      U2[i] =   sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
      U3[i] = sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) * Normal(gen);

     if(dt >1e-15){
           dt_remain = fabs ( (x2[i]-yw)/temp_U0 );
           x1[i] =  xw + dt_remain*U1[i];
           x2[i] =  yw + dt_remain*U2[i];
           dt = dt - dt_remain;
           x1_old[i] = xw;
           x2_old[i] = yw;
      }
      else{
//            x2[i] = x2[i] + 2.0*fabs( x2[i] );
            x2[i] = yw;
            //x1[i] = xw;
      }
      if(x1[i] - (*box).Len[0] < 0.0 && x1[i] > 0.0){
            xw = line_interesect_x(yw-2.0*epsilon, x1_old[i], x2_old[i], x1[i], x2[i]);
            update_face_info_BC(U1[i],U2[i],U3[i], x1[i], x2[i], xw, yw-2.0*epsilon, gas, box, cells);
      }
      M[3] = M[3] +  (*gas).Fn*(*gas).m*( fabs(U2[i])+fabs(temp_U0) );
      num_hitting_wall++;
      flag[i] += 1000;
    }

    if(x1[i] - (*box).Len[0] < 0.0 && x1[i] > 0.0 && x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0)
      done = 1;
//     else
//       printf("one particle hits both walls in time step\n");
    dummy ++;
//     if(dt_remain > (*gas).delta_t && ((*gas).model == "FP" || (*gas).model == "FP_ideal_gas") )
//       printf("ooops, dt_remain>delta_t\n");
  }
}

  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
//   printf("***   number of wall hits = %d\n",num_hitting_wall);
}

void BC_periodic(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells){

//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;
int num_hitting_wall = 0;
double Tol1 = 0.0;//1e-16;
double Tol2 = 0.0;//-1e-16;
double epsilon = 1e-14;

  std::random_device rd;
    std::random_device rd2;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);

  int N1, N2, N3;
  int ix,iy,iz;
  int num_L, num_U, id, num, num2;
  int i,j,k,ii;
  //double eqdist = 2.2*pow(1.0/(*gas).n,1.0/3.0);//
  double h = 1.01*2.0*(*gas).eqdist;
  double fmid = 1.0/4.0;
  N1 = floor((*box).Len[0]/h);//ceil((*box).Len[0]/(eqdist));
  N2 = floor((*box).Len[1]/h*fmid);//ceil(0.5*(*box).Len[1]/(eqdist));
  N3 = floor((*box).Len[2]/h);//ceil((*box).Len[2]/(eqdist));
  //eqdist = (*box).Len[0]/(1.0*N1);
  double hh[3];
  hh[0] = (*box).Len[0]/(1.0*N1);
  hh[1] = fmid*(*box).Len[1]/(1.0*N2);
  hh[2] = (*box).Len[2]/(1.0*N3);
  int *N = (int *) malloc( N1*N2*N3 * sizeof(int) );
  for(num=0; num<N1*N2*N3; num++){
	N[num] = 0;
  }
  num_L = floor( (3.0/8.0)*(*box).Len[1]/(*box).delta_dim[1] );
  num_U = ceil(  (5.0/8.0)*(*box).Len[1]/(*box).delta_dim[1] );
  int not_in=0;
  for(num = num_L; num <= num_U; num++){
	for(ii=0; ii<cells[num].num_inside; ii++){
		id = cells[num].indices_inside[ii];
		if(x1[id] > 0.0 && x1[id]<(*box).Len[0] && x2[id] > 0.5*(*box).Len[1] - N2*hh[1]/2.0 + 1e-15&& x2[id] < 0.5*(*box).Len[1] + N2*hh[1]/2.0 - 1e-15&& x3[id]>0.0 && x3[id]<(*box).Len[2]){
			i = floor( x1[id]/(hh[0]));
			j = floor((x2[id]-(3.0/8.0)*(*box).Len[1])/(hh[1]));
			k = floor( x3[id]/(hh[2]) );
			num2 = i*N3*N2+j*N3+k;
			if(num2 <0 || num2 > N1*N2*N3-1)
				printf("num2 is out, num2 = %d, i=%d, j=%d, k=%d; N1=%d, N2=%d, N3=%d\n",num2, i, j, k, N1, N2, N3);
			else
				N[num2] = N[num2]+1;
		}
		else
			not_in ++;
	}
  }

  // checking 1st dimension
int num_zeros = 0;
for(i=0; i<N1*N2*N3; i++){
	if(N[i] == 0){
		num_zeros += 1.0;
	}
}

int count = 0;
int *index_zeros;

index_zeros = (int *) malloc( num_zeros * sizeof(int) );
if(num_zeros > 0){
	for(i=0; i < N1*N2*N3; i++){
		if(N[i] == 0){
			index_zeros[count] = i;
			count ++;
		}
	}
}

int num_zero_now = num_zeros;
int done;
int many;
ix = 0; iy = 0; iz = 0;
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  while (done == 0){
    if (dummy > 10){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %e,  %e, %e \n", x1[i], x2[i], x3[i]);
      printf(" L = %e,  %e, %e \n", (*box).Len[0], (*box).Len[1], (*box).Len[2]);
      if(dummy > 30){
	exit(1);
      }
    }
    many = 0;
    while(x1[i] > (*box).Len[0] - 1.0e-15 ){
//      U1[i] = - U1[i];
//      x1[i] = x1[i] - (*box).Len[0];
      x1[i] = x1[i] - (*box).Len[0]*floor(x1[i]/(*box).Len[0]) + epsilon;//*floor(x1[i]/(*box).Len[0]);
//      x1[i] = x1[i] - floor( fabs(x1[i]/(*box).Len[0]) )*(*box).Len[0];
      //x1[i] = x1[i] - (*box).Len[0];
      M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      many++;
      if(many>20)
	x1[i] = epsilon;
    }
    many = 0;
    while(x1[i] < 1.0e-15 ){
//      U1[i] = - U1[i];
      //x1[i] = x1[i] + (*box).Len[0];
      x1[i] = x1[i] + (*box).Len[0]*(floor(fabs(x1[i]/(*box).Len[0]))+1) - epsilon;
      //x1[i] = x1[i] + (floor( fabs(x1[i]/(*box).Len[0]) )+1)*(*box).Len[0];
      //x1[i] = x1[i] + (*box).Len[0];
      M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      many ++;
      if(many>20)
	x1[i] = (*box).Len[0] - epsilon;
    }
    while(x2[i] > (*box).Len[1] - 1e-15) {
      if((*box).step > (*box).after){
		(*box).out += 1.0;
      }
//      U2[i] = sqrt(2.0* (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
//      U1[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );
 //     U3[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );
 //     x2[i] = (*box).Len[1]+U2[i]*(*gas).delta_t;
      //U2[i] = - U2[i];
//      x2[i] = x2[i] - (*box).Len[1];
	if( num_zero_now>0 ){
		done = 0;
		while(done == 0){
			std::uniform_int_distribution<> rnd_cell(0, num_zero_now-1);
			ii = rnd_cell(gen);
		  	num2 = index_zeros[ ii ];
			if(N[num2]==0){
				done = 1;
				N[num2] = 1;
				for(j=ii; j<num_zero_now-1; j++){
					index_zeros[j] = index_zeros[j+1];
				}
				num_zero_now --;
			}
			else
				printf("something wrong in removing element\n");
		}
		iz = num2%N3;
		iy = ( num2%(N2*N3)-iz )/N3;
		ix = floor( ( num2-iz-iy*N3)/(N2*N3) );
		//printf("ix = %d\n",ix);
		//if(ix == 30){
		//	printf("num2=%d, ix=%d, iy=%d, iz=%d\n",num2, ix,iy,iz);
		//	printf("num2 30 = %d\n", num2%30);
		//	printf("N1, N2, N3 = %d %d %d\n", N1, N2, N3);
		//	exit(1);
		//}
		x1[i] = (ix+0.5)*hh[0];
		x2[i] = (3.0/8.0)*(*box).Len[1] + (iy+0.5)*hh[1];
		x3[i] = (iz+0.5)*hh[2];

		num = floor(x2[i]/(*box).delta_dim[1]);
		if( x1[i]>(*box).Len[0] - epsilon ){
			printf("HERE, particle %d is out from right, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
			printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    		}
		if( x1[i] <  epsilon ){
			printf("HERE, particle %d is out from left, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
			printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    		}

		U1[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		U2[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		U3[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		//printf("reinitialization successful\n");
	}
	else{
	      x2[i] = x2[i] - (*box).Len[1]*floor(x2[i]/(*box).Len[1]) +epsilon;
		//printf("reinitialization unsuccessful\n");
	}
	//x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] );
      //x2[i] = x2[i] - floor( fabs(x2[i]/(*box).Len[1]) )*(*box).Len[1];
      //x2[i] = x2[i] - (*box).Len[1];
      M[2] = M[2] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
    while(x2[i] <  1e-15){
      if((*box).step > (*box).after){
		(*box).out += 1.0;
      }
//U2[i] = -sqrt(2.0* (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * sqrt(-log(uniform(gen2)));
//      U1[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );
//      U3[i] = sqrt( (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * Normal ( gen );
     // x2[i] = 0.0+U2[i]*(*gas).delta_t;
//      U2[i] = - U2[i];
//      x2[i] = x2[i] + (*box).Len[1];
//x2[i] = x2[i] + 2.0*fabs( x2[i]  );
	if( num_zero_now>0 ){
		done = 0;
		while(done == 0){
			std::uniform_int_distribution<> rnd_cell(0, num_zero_now-1);
			ii = rnd_cell(gen);
		  	num2 = index_zeros[ ii ];
			if(N[num2]==0){
				done = 1;
				N[num2] = 1;
				for(j=ii; j<num_zero_now-1; j++){
					index_zeros[j] = index_zeros[j+1];
				}
				num_zero_now --;
			}
			else
				printf("something wrong in removing element\n");
		}
		iz = num2%N3;
		iy = ( num2%(N2*N3)-iz )/N3;
		ix = floor( ( num2-iz-iy*N3)/(N2*N3) );
		//printf("ix = %d\n",ix);
		//if(ix == 30){
		//	printf("num2=%d, ix=%d, iy=%d, iz=%d\n",num2, ix,iy,iz);
		//	printf("num2 30 = %d\n", num2%30);
		//	printf("N1, N2, N3 = %d %d %d\n", N1, N2, N3);
		//	exit(1);
		//}
		x1[i] = (ix+0.5)*hh[0];
		x2[i] = (3.0/8.0)*(*box).Len[1] + (iy+0.5)*hh[1];
		x3[i] = (iz+0.5)*hh[2];
		num = floor(x2[i]/(*box).delta_dim[1]);

		if( x1[i]>(*box).Len[0]- epsilon ){
			printf("HERE, particle %d is out from right, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
			printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    		}
		if( x1[i] <  epsilon ){
			printf("HERE, particle %d is out from left, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
			printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    		}
		U1[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		U2[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		U3[i] = sqrt( (*gas).kb*cells[num].T/( (*gas).m) ) * Normal ( gen );
		//printf("reinitialization successful\n");
	}
	else{
      		x2[i] = x2[i] + (*box).Len[1]*(floor(fabs(x2[i]/(*box).Len[1]))+1) -epsilon;
		//printf("reinitialization unsuccessful\n");
	}
      //x2[i] = x2[i] + ( floor( fabs(x2[i]/(*box).Len[1]) )+1)*(*box).Len[1];
      //x2[i] = x2[i] + (*box).Len[1];
      M[3] = M[3] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
    while(x3[i] > (*box).Len[2] -1e-15) {
      //U3[i] = - U3[i];
//      x3[i] = x3[i] - (*box).Len[2];
      x3[i] = x3[i] - (*box).Len[2]*floor(x3[i]/(*box).Len[2])+epsilon;
      //x2[i] = x2[i] - floor( fabs(x2[i]/(*box).Len[1]) )*(*box).Len[1];
      //x2[i] = x2[i] - (*box).Len[1];
      M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
    while(x3[i] < 1e-15){
      //U3[i] = - U3[i];
      //x3[i] = x3[i] + (*box).Len[2];
      x3[i] = x3[i] + (*box).Len[2]*(floor(fabs(x3[i]/(*box).Len[2]))+1)-epsilon;
      //x2[i] = x2[i] + ( floor( fabs(x2[i]/(*box).Len[1]) )+1)*(*box).Len[1];
      //x2[i] = x2[i] + (*box).Len[1];
      M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
    if(x1[i] - (*box).Len[0] < -1e-15 && x1[i] > 1e-15 && x2[i] - (*box).Len[1] < -1.0e-15 && x2[i] > 1.0e-15 && x3[i] - (*box).Len[2] < -1.0e-15 && x3[i] > 1.0e-15)
      done = 1;
    if( x1[i]>(*box).Len[0] ){
		printf("in Periodic, OOPS, particle %d is out from right, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
		printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    }
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
//   printf("***   number of wall hits = %d\n",num_hitting_wall);

free(N); free(index_zeros);
}


void BC_specular_periodic_reflection(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells){

//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;
int num_hitting_wall = 0;
double Tol1 = 0.0;//1e-16;
double Tol2 = 0.0;//-1e-16;
double epsilon = 1e-14;


  int N1, N2, N3;
  int ix,iy,iz;
  int num_L, num_U, id, num, num2;
  int i,j,k,ii;

int count = 0;



int done;
int many;
ix = 0; iy = 0; iz = 0;
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  while (done == 0){
    if (dummy > 10){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %e,  %e, %e \n", x1[i], x2[i], x3[i]);
      printf(" L = %e,  %e, %e \n", (*box).Len[0], (*box).Len[1], (*box).Len[2]);
      if(dummy > 30){
	exit(1);
      }
    }
    many = 0;
    while(x1[i] > (*box).Len[0] - 1.0e-15 ){
      x1[i] = x1[i] - (*box).Len[0]*floor(x1[i]/(*box).Len[0]) + epsilon;
      M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      many++;
      if(many>20)
	x1[i] = epsilon;
    }
    many = 0;
    while(x1[i] < 1.0e-15 ){
      x1[i] = x1[i] + (*box).Len[0]*(floor(fabs(x1[i]/(*box).Len[0]))+1) - epsilon;
      M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
      many ++;
      if(many>20)
	x1[i] = (*box).Len[0] - epsilon;
    }
    while(x2[i] > (*box).Len[1] - 1e-15) {
      if((*box).step > (*box).after){
		(*box).out += 1.0;
      }
      x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] );
      //x2[i] = x2[i] - (*box).Len[1]*floor(x2[i]/(*box).Len[1]) +epsilon;
      M[2] = M[2] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
    while(x2[i] <  1e-15){
      if((*box).step > (*box).after){
		(*box).out += 1.0;
      }
      x2[i] = x2[i] + 2.0*fabs( x2[i]  );
      //x2[i] = x2[i] + (*box).Len[1]*(floor(fabs(x2[i]/(*box).Len[1]))+1) -epsilon;
      M[3] = M[3] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
    while(x3[i] > (*box).Len[2] -1e-15) {
      //U3[i] = - U3[i];
//      x3[i] = x3[i] - (*box).Len[2];
      x3[i] = x3[i] - (*box).Len[2]*floor(x3[i]/(*box).Len[2])+epsilon;
      //x2[i] = x2[i] - floor( fabs(x2[i]/(*box).Len[1]) )*(*box).Len[1];
      //x2[i] = x2[i] - (*box).Len[1];
      M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
    while(x3[i] < 1e-15){
      //U3[i] = - U3[i];
      //x3[i] = x3[i] + (*box).Len[2];
      x3[i] = x3[i] + (*box).Len[2]*(floor(fabs(x3[i]/(*box).Len[2]))+1)-epsilon;
      //x2[i] = x2[i] + ( floor( fabs(x2[i]/(*box).Len[1]) )+1)*(*box).Len[1];
      //x2[i] = x2[i] + (*box).Len[1];
      M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
    if(x1[i] - (*box).Len[0] < -1e-15 && x1[i] > 1e-15 && x2[i] - (*box).Len[1] < -1.0e-15 && x2[i] > 1.0e-15 && x3[i] - (*box).Len[2] < -1.0e-15 && x3[i] > 1.0e-15)
      done = 1;
    if( x1[i]>(*box).Len[0] ){
		printf("in Periodic, OOPS, particle %d is out from right, x1=%e x1>L=%e\n",i, x1[i],(*box).Len[0]);
		printf("num = %d, num2=%d, ix=%d, iy=%d, iz=%d\n",num, num2, ix,iy,iz);
    }
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
//   printf("***   number of wall hits = %d\n",num_hitting_wall);

}


void BC_specular_reflection(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box){

//   printf("BC Starts\n");
  double M[NSD*2];
for(int i=0; i<NSD*2; i++)
  M[i]= 0.0;
for(int i=0; i<NSD*2; i++)
  (*gas).M[i]= 0.0;
int num_hitting_wall = 0;
double Tol1 = 1e-16;
double Tol2 = -1e-16;
double epsilon = 1e-15;
  // checking 1st dimension
for(int i=0; i<(*gas).N; i++){
  int done = 0;
  int dummy = 0;
  while (done == 0){
    if (dummy > 5){
      printf("*** TOO MANY BC OPERATIONS ***       %d\n", dummy);
      printf(" point %d \n", i);
      printf(" with coord  %e,  %e, %e \n", x1[i], x2[i], x3[i]);
//      scanf("%d", &dummy);
    }
    /*
    if(x1[i] - (*box).Len[0] > -epsilon && dummy > 5){
          x1[i] = (*box).Len[0] - epsilon;
    }
    else if(x1[i] < epsilon && dummy > 5){
          x1[i] = epsilon;
    }
    if(x2[i] - (*box).Len[1] > -epsilon && dummy > 5){
          x2[i] = (*box).Len[1] - epsilon;
    }
    else if(x2[i] < epsilon && dummy > 5){
          x2[i] = epsilon;
    }
    if(x3[i] - (*box).Len[2] > -epsilon && dummy > 5){
          x3[i] = (*box).Len[2] - epsilon;
    }
    else if(x3[i] < epsilon && dummy > 5){
          x3[i] = epsilon;
    }
*/
  if( (*box).direction[0] == 1 ){
    if(x1[i] - (*box).Len[0] > Tol2){
      U1[i] = - U1[i];
//      x1[i] = x1[i] - (*box).Len[0];
      x1[i] = x1[i] - 2.0*fabs( x1[i] - (*box).Len[0] );
//      x1[i] = x1[i] - floor( fabs(x1[i]/(*box).Len[0]) )*(*box).Len[0];
      //x1[i] = x1[i] - (*box).Len[0];
      M[0] = M[0] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
    }
    else if(x1[i] < Tol1){
      U1[i] = - U1[i];
      //x1[i] = x1[i] + (*box).Len[0];
      x1[i] = x1[i] + 2.0*fabs( x1[i]  );
      //x1[i] = x1[i] + (floor( fabs(x1[i]/(*box).Len[0]) )+1)*(*box).Len[0];
      //x1[i] = x1[i] + (*box).Len[0];
      M[1] = M[1] + (*gas).Fn*fabs(2*(*gas).m*U1[i]);
      num_hitting_wall++;
    }
  }
  if( (*box).direction[1] == 1 ){
    if(x2[i] - (*box).Len[1] > Tol2) {
      U2[i] = - U2[i];
			U1[i] =  U1[i];// + 100.0;
//      x2[i] = x2[i] - (*box).Len[1];
      x2[i] = x2[i] - 2.0*fabs( x2[i] - (*box).Len[1] );
      //x2[i] = x2[i] - floor( fabs(x2[i]/(*box).Len[1]) )*(*box).Len[1];
      //x2[i] = x2[i] - (*box).Len[1];
      M[2] = M[2] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
    else if(x2[i] < Tol1){
      U2[i] = - U2[i];
			U1[i] =  U1[i];// - 100.0;
//      x2[i] = x2[i] + (*box).Len[1];
      x2[i] = x2[i] + 2.0*fabs( x2[i]  );
      //x2[i] = x2[i] + ( floor( fabs(x2[i]/(*box).Len[1]) )+1)*(*box).Len[1];
      //x2[i] = x2[i] + (*box).Len[1];
      M[3] = M[3] + (*gas).Fn*fabs(2*(*gas).m*U2[i]);
      num_hitting_wall++;
    }
 }
 if( (*box).direction[2] == 1 ){
    if(x3[i] - (*box).Len[2] > Tol2) {
      U3[i] = - U3[i];
//      x3[i] = x3[i] - (*box).Len[2];
      x3[i] = x3[i] - 2.0*fabs( x3[i] - (*box).Len[2] );
      //x2[i] = x2[i] - floor( fabs(x2[i]/(*box).Len[1]) )*(*box).Len[1];
      //x2[i] = x2[i] - (*box).Len[1];
      M[4] = M[4] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
    else if(x3[i] < Tol1){
      U3[i] = - U3[i];
      //x3[i] = x3[i] + (*box).Len[2];
      x3[i] = x3[i] + 2.0*fabs( x3[i]  );
      //x2[i] = x2[i] + ( floor( fabs(x2[i]/(*box).Len[1]) )+1)*(*box).Len[1];
      //x2[i] = x2[i] + (*box).Len[1];
      M[5] = M[5] + (*gas).Fn*fabs(2*(*gas).m*U3[i]);
      num_hitting_wall++;
    }
  }
  done = 1;
  if( (*box).direction[0] == 1 ){
	if(x1[i] - (*box).Len[0] > 0.0 || x1[i] < 0.0)
		done = 0;
  }
  if( (*box).direction[1] == 1 ){
	if(x2[i] - (*box).Len[1] > 0.0 || x2[i] < 0.0)
		done = 0;
  }
  if( (*box).direction[2] == 1 ){
	if(x3[i] - (*box).Len[2] > 0.0 || x3[i] < 0.0)
		done = 0;
  }
 //   if(x1[i] - (*box).Len[0] < 0.0 && x1[i] > 0.0 && x2[i] - (*box).Len[1] < 0.0 && x2[i] > 0.0 && x3[i] - (*box).Len[2] < 0.0 && x3[i] > 0.0)
 //     done = 1;
    dummy ++;
  }
}
  for(int i=0; i<NSD*2; i++)
    (*gas).M[i] +=  M[i];
//   printf("***   number of wall hits = %d\n",num_hitting_wall);
}

void find_inside_cells(double *x1,double *x2,double *x3,double gas_N, struct BOX *box, struct CELLS *cells, int *index){

//   int number_insides[(*box).N[0]*(*box).N[1]*(*box).N[2]];
  int i, j,k,l, num;
  for (i=0; i < (*box).N[0]*(*box).N[1]*(*box).N[2]; i++)
	  cells[i].num_inside = 0;
  for (num=0; num < (*box).N[0]*(*box).N[1]*(*box).N[2]; num++)
    for (i=0; i < gas_N; i++)
      cells[num].indices_inside[ i ] = -1;
  for (i=0; i < gas_N; i++){
    j = floor(x1[i]/(*box).delta_dim[0]);
    k = floor(x2[i]/(*box).delta_dim[1]);
    l = floor(x3[i]/(*box).delta_dim[2]);
    if(j<0 || k<0 || l<0 || j>(*box).N[0]-1 || k>(*box).N[1]-1 || l>(*box).N[2]-1)
      printf("OOPS, particle is still outside the box, probelm with BC");
    num = l*(*box).N[0]*(*box).N[1] + k*(*box).N[0] + j;
    index[i] = num;
    cells[num].indices_inside[ cells[num].num_inside ] = i;
    cells[num].num_inside = cells[num].num_inside+1;
  }
//     for (int j=0; j < (*box).N[0]; j++)
//       for (int k=0; k < (*box).N[1]; k++)
//         for (int l=0; l < (*box).N[2]; l++){
// 	  int num = l*(*box).N[0]*(*box).N[1] + k*(*box).N[0] + j;
// 	  if(x1[i] > cells[num].dim[0] && x1[i] < cells[num].dim[1])
// 	    if(x2[i] > cells[num].dim[2] && x2[i] < cells[num].dim[3])
// 	      if(x3[i] > cells[num].dim[4] && x3[i] < cells[num].dim[5]){
// 		cells[num].indices_inside[ number_insides[num] ] = i;
// 		number_insides[num]++;}
// 	}
//   for ( j=0; j < (*box).N[0]; j++)
//     for ( k=0; k < (*box).N[1]; k++)
//        for ( l=0; l < (*box).N[2]; l++){
// 	  num = l*(*box).N[0]*(*box).N[1] + k*(*box).N[0] + j;
// 	  cells[num].num_inside = number_insides[num];}
//
// }
}
