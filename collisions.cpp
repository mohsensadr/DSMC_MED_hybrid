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
#include <omp.h>
#include "omp.h"
*/

void SPH_simple(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  int N = (*gas).N, i, j, id;
  int search, kk, kkk, num;
	double W, sqrt_pi, mb, rr, r1r2, r2, r1, h, r_cut;
  double R = (*gas).kb/(*gas).m;
  double dpddx, D2W, DW, k, mu, p, dusdx, gama, dTdt, d2Tdx2, dudx, drhodx, dmudx, d2dx2, ddx;
  double dpdx, dTdx, dkdx, mu0;
  double cv;
  cv = 3.0*(*gas).kb/(*gas).m/2.0;
  mu0 = 5.0/( 16.0*(*gas).sigma*(*gas).sigma )*sqrt( (*gas).m*(*gas).kb*(*gas).T0/acos(-1.0) );
  gama = 1.4;
  r_cut = (*gas).r_cut;
  h = (*gas).h;
  mb = (*gas).mb;
  search = floor(r_cut/(*box).delta_dim[1])+1;
  sqrt_pi = sqrt(2.0*acos(-1.0));

  double *F = (double *) malloc( N * sizeof(double) );
  double *G = (double *) malloc( N * sizeof(double) );

  for(i=0; i<N; i++){
    F[i] = 0.0;
    G[i] = 0.0;
    r1 = x2[i];
    num = index[i];
    p = rho[i]/R*T[i];
    mu = mu0*pow(T[i]/(*gas).T0,0.5);
    k = 15.0*(*gas).kb*mu/(4.0*(*gas).m);
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
            D2W = W*( r1r2*r1r2/pow(h,4.0) - r1r2/rr/h/h );
          }
          else{
            DW = 0.0; D2W = 0.0;
          }

          F[i] +=  mb*R*(T[i]/rho[i]+T[j]/rho[j])*DW;
          G[i] +=  0.5*mb*R*( T[i]/rho[i] + T[j]/rho[j] )*( U2[j]-U2[i] )*DW;
    		}
    	}
    }
  }



  for(i=0; i<N; i++){
      U2[i] = U2[i] + F[i]*(*gas).delta_t;
      T[i]  = T[i]  + G[i]*(*gas).delta_t/cv;
  }
  free(F); free(G);
}


void SPH_velocity_temperature_update(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  // first solve v*
  //double **A0;
  //double **A;
  int N = (*gas).N, i, j, id;
  int search, kk, kkk, num;
	double W, sqrt_pi, mb, rr, r1r2, r2, r1, h, r_cut;
  double R = (*gas).kb/(*gas).m;
  double dpddx, D2W, DW, k, mu, p,  gama,  dTdt, dkdx;
  double  mu0;
  double tol = 1e-15;
  mu0 = 5.0/( 16.0*(*gas).sigma*(*gas).sigma )*sqrt( (*gas).m*(*gas).kb*(*gas).T0/acos(-1.0) );
  gama = 1.4;
  r_cut = (*gas).r_cut;
  h = (*gas).h;
  mb = (*gas).mb;
  search = floor(r_cut/(*box).delta_dim[1])+1;
  sqrt_pi = sqrt(2.0*acos(-1.0));
  //A=new double *[N];
  //A0=new double *[N];
  //for(i=0;i<N;i++){
  //    A0[i]=new double [N];/
//
//      A[i]=new double [N];
//  }

  double *rA0 = (double *) malloc( N * sizeof(double) );
  double *rA = (double *) malloc( N * sizeof(double) );
  double *vs   = (double *) malloc( N * sizeof(double) );
  double *pnp1   = (double *) malloc( N * sizeof(double) );
  double *dmudx   = (double *) malloc( N * sizeof(double) );


  double res;
  double *drhodx = (double *) malloc( N * sizeof(double) );
  double *dusdx = (double *) malloc( N * sizeof(double) );
  double *dudx = (double *) malloc( N * sizeof(double) );
  double *d2Tdx2 = (double *) malloc( N * sizeof(double) );
  double *dTdx = (double *) malloc( N * sizeof(double) );
  double *dpdx = (double *) malloc( N * sizeof(double) );

  double *A = (double *) malloc( N *N* sizeof(double) );
  double *A0 = (double *) malloc( N *N* sizeof(double) );
  double *ddx = (double *) malloc( N *N* sizeof(double) );
  double *d2dx2 = (double *) malloc( N*N * sizeof(double) );

  for(i=0; i<N; i++){
      for(j=0; j<N; j++){
          //A[i][j] = 0.0;
          A[i+j*N] = 0.0;
          ddx[i+j*N] = 0.0;
          d2dx2[i+j*N] = 0.0;
      }
  }
  for(i=0; i<N; i++){
      rA[i] = 0.0;
  }
  for(i=0; i<N; i++)
      vs[i] = 0.0;
  for(i=0; i<N; i++)
      pnp1[i] = 0.0;

  int nn = 0;


  for(i=0; i<N; i++){
    ddx[i] = 0.0;
    d2dx2[i] = 0.0;
    dmudx[i] = 0.0;
    drhodx[i] = 0.0;
    dudx[i] = 0.0;
    d2Tdx2[i] = 0.0;
    dTdx[i] = 0.0;
    r1 = x2[i];
    num = index[i];
    p = rho[i]/R*T[i];
    mu = (*gas).mu*pow(T[i]/(*gas).T0,0.5);
    k = (*gas).k/(*gas).mu*mu;
    //k = 15.0*(*gas).kb*mu/(4.0*(*gas).m);
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

          D2W = W*( -1.0+rr*rr/h/h )/h/h;
          if(i!=j){
            DW = W*r1r2*(-1.0/h/h);
            //D2W = W*( r1r2*r1r2/pow(h,4.0) - r1r2/rr/h/h );
          }
          else{
            DW = 0.0;
          }
          dmudx[i] += mb/rho[j]*DW*mu0*pow(T[j]/(*gas).T0,0.5);

          ddx[i+j*N] = mb/rho[j]*DW;
          d2dx2[i+j*N] = mb/rho[j]*D2W;

          drhodx[i] += mb*DW;
          dudx[i]   += mb*U2[j]/rho[j]*DW;
          dTdx[i]   += mb*T[j]/rho[j]*DW;
          d2Tdx2[i] += mb*T[j]/rho[j]*D2W;

    		}
    	}
    }
  }



  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      if(i==j){
        A[i+j*N] = 1.0;
      }
      A[i+j*N] += -(*gas).delta_t*(4.0/3.0)/rho[i]*dmudx[i]*ddx[i+j*N];
      A[i+j*N] += -(*gas).delta_t*(4.0/3.0)/rho[i]*mu*d2dx2[i+j*N];
    }
    rA[i] = U2[i];
  }

  int nrhs = 1;
  int lda=N;
  int ipiv[N];
  int ldb = N;
  int info;

  double nrA;
  /*
  nrA = 0.0;
  for(i=0; i<N; i++){
      for(j=0; j<N; j++){
          A0[i*N+j] = A[i+j*N];
        }
      }
  for(i=0; i<N; i++){
    rA0[i] = rA[i];
    nrA += rA[i];
  }
  */
  //dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);
  for(i=0; i<N; i++){
    vs[i] = rA[i];
  }
  /*
  res = 0.0;
  for(i=0; i<N; i++){
      res += rA0[i];
      for(j=0; j<N; j++){
        res += -A0[i*N+j]*vs[j];
      }
  }
  */
  ////////////////////////
  ///////////////////////          solve p n+1
  ///////////////////////


  for(i=0; i<N; i++){
    dusdx[i] = 0.0;
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
          if(i!=j){
            DW = W*r1r2*(-1.0/h/h);
          }
          else{
            DW = 0.0;
          }
          dusdx[i]  += mb*vs[j]/rho[j]*DW;
    		}
    	}
    }
  }


  for(i=0; i<N; i++){
    p = rho[i]/R*T[i];
    mu = mu0*pow(T[i]/(*gas).T0,0.5);
    k = 15.0*(*gas).kb*mu/(4.0*(*gas).m);

    dkdx = k/mu*dmudx[i];
    dTdt = (gama-1.0)/rho[i]/R*( -p*dudx[i] + (4.0/3.0)*mu*dudx[i]*dudx[i] + dkdx*dTdx[i] + k*d2Tdx2[i] );
    rA [i] = -rho[i]/(*gas).delta_t;
    rA [i] += rho[i]*dusdx[i];
    rA [i] += -rho[i]/T[i]*dTdt;

    for(j=0; j<N; j++){
      A[i+j*N] = 0.0;
      if(i==j){
          A[i+j*N] = -1.0/( (*gas).delta_t*R*T[i] );
      }
      A[i+j*N] += -(*gas).delta_t/rho[i]*drhodx[i]*ddx[i+j*N];
      A[i+j*N] += (*gas).delta_t*d2dx2[i+j*N];
    }
  }

  //Gauss_Jordan(B, rB, N, pnp1);

/*
  nrA = 0.0;
  for(i=0; i<N; i++){
      for(j=0; j<N; j++){
          A0[i*N+j] = A[i+j*N];
        }
      }
  for(i=0; i<N; i++){
    rA0[i] = rA[i];
    nrA += rA[i];
  }
*/
  //for(i=0; i<N; i++)
  //  for(j=0; j<N; j++)
  //    a[i+j*N] = A[i*N+j];
  //dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);

  for(i=0; i<N; i++)
    pnp1[i] = rA[i];

/*
  res = 0.0;
  for(i=0; i<N; i++){
      res += rA0[i];
      for(j=0; j<N; j++){
        res += -A0[i*N+j]*pnp1[j];
      }
  }
*/


  // update v n+1
  for(i=0; i<(*gas).N; i++){
    dpdx[i] = 0.0;

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
          if(i!=j){
            DW = W*r1r2*(-1.0/h/h);
          }
          else{
            DW = 0.0;
          }

          dpdx[i]   += mb*pnp1[j]/rho[j]*DW;

        }
      }
    }
  }

  for(i=0; i<N; i++){
    U2[i]   = vs[i] + (*gas).delta_t/rho[i]*dpdx[i];
  }



  for(i=0; i<N; i++){
    p = rho[i]/R*T[i];
    mu = mu0*pow(T[i]/(*gas).T0,0.5);
    k = 15.0*(*gas).kb*mu/(4.0*(*gas).m);
    for(j=0; j<N; j++){
      dkdx = k/mu*dmudx[i];
      //dTdt = (gama-1.0)/rho[i]/R*( -p*dudx[i] + (4.0/3.0)*mu*dudx[i]*dudx[i] + dkdx*dTdx[i] + k*d2Tdx2[i] );

      A[i+j*N] = 0.0;
      if(i==j){
          A[i+j*N] = 1.0;
      }
      A[i+j*N] += -(*gas).delta_t/rho[i]/R*(gama-1.0)*dkdx*ddx[i+j*N];
      A[i+j*N] += -(*gas).delta_t/rho[i]/R*(gama-1.0)*k*d2dx2[i+j*N];
    }
    rA [i] = T[i] + (*gas).delta_t*(gama-1.0)/R/rho[i]*( -p*dudx[i] + 4.0/3.0*mu*dudx[i]*dudx[i] );
  }


  /*
  nrA = 0.0;
  for(i=0; i<N; i++){
      for(j=0; j<N; j++){
          A0[i*N+j] = A[i+j*N];
        }
      }
  for(i=0; i<N; i++){
    rA0[i] = rA[i];
    nrA += rA[i];
  }
*/
//  dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);

  for(i=0; i<N; i++)
    T[i] = rA[i];
/*
  res = 0.0;
  for(i=0; i<N; i++){
      res += rA0[i];
      for(j=0; j<N; j++){
        res += -A0[i*N+j]*T[j];
      }
  }
*/




//  for(i=0; i<N; i++)
//      T[i] = 0.0;
  //Gauss_Jordan(C, rC, N, T);
  // solve T n+1

  // update rho
  for(i=0; i<N; i++)
      rho[i] = pnp1[i]/R/T[i];
  //for(i=0;i<N;i++){
  //     delete [] A[i];
  //}
  //delete  [] A;
  free(rA); free(A);// free(a);
  free(vs); free(pnp1);
  free(rA0); free(dmudx);
  free(dusdx); free(dudx); free(drhodx);
  free(d2Tdx2); free(dTdx); free(dpdx);
  free(A0); free(ddx); free(d2dx2);

}
void Collision_ESMC_dcones(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::uniform_real_distribution<> dis_y(-1.0, 1.0);
  std::uniform_real_distribution<> dis_z(-1.0, 1.0);

  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::uniform_real_distribution<> dis_rr(0.0, 1.0);

  std::uniform_real_distribution<> dis_u(0.0, 1.0);
  std::uniform_real_distribution<> dis_v(0.0, 1.0);
  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[3], x_mid[3];
  int  ind[3], ind_mid[3];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double omega, r, Vc;

  double dummy;
  int num_act_col, cand_num;
  int num, Nc, id, num2;
  double M_cand;
  int in, out;

  double rr;
  int sgn;
  int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
  int Ny = (*box).Ny;
  int num_with_max_coll = 0;
  int dummy2=0;
  double X, Y, Xm, Ym;
  double theta, phi, u,v;
  for( num = 0; num<Nx*Ny; num++){
     omega = 0.0;

     Nc = cells[num].num_inside;
     Vc = cells[num].volume;
     M_cand = 0.5*cells[num].omega_max*( 1.0*Nc );
     num_act_col = 0;
     if( Nc>0){
	std::uniform_int_distribution<> disss1(0, cells[num].num_inside-1);
        for( cand_num=0; cand_num < floor(M_cand); cand_num++)
        {
	    out = 0;
 	    i = cells[num].indices_inside[ disss1(gen) ];
	    I = num;
            next = 0;
	    while(next ==0){
		    // picking points on surface of sphere from a uniform distribution
		    // http://mathworld.wolfram.com/SphericalCoordinates.html
		    u = dis_u(gen);
		    v = dis_v(gen);
		    theta = 2.0*pi*u;
		    phi = acos(2.0*v-1.0);
		    dir[0] = cos(theta)*sin(phi);
		    dir[1] = sin(theta)*sin(phi);
		    dir[2] = cos(phi);
/*
	            dir[0] = dis_x(gen);
        	    dir[1] = dis_y(gen);
        	    dir[2] = dis_z(gen);
        	    size = sqrt(dir[0]*dir[0] + dir[1]*dir[1]+dir[2]*dir[2]);
        	    for(ii=0; ii<3; ii++)
		          dir[ii] = dir[ii]/size;
*/
            	    // finding the cell at r + sigma*dir
            	    x[0] = x1[i] + (*gas).sigma*dir[0];
            	    x[1] = x2[i] + (*gas).sigma*dir[1];
		    x[2] = (*gas).sigma*dir[2];
		    X = x[0];
		    Y = sqrt(x[1]*x[1]+x[2]*x[2]);
            	    // if the other partilce is outside the domain, go for the next particle'
             	    J = dcone_index(X, Y, cells,  box, &ind[0], &ind[1]);

	    	   if(ind[0]>=0 && ind[0]<Nx && ind[1]>=0 && ind[1]<Ny && cells[J].num_inside>0){
			// finding the cell at r + sigma/2*dir
			x_mid[0] = x1[i] + (*gas).sigma*dir[0]/2.0;
			x_mid[1] = x2[i] + (*gas).sigma*dir[1]/2.0;
			x_mid[2] = 0.0 + (*gas).sigma*dir[2]/2.0;

			Xm = x_mid[0];
			Ym = sqrt(x_mid[1]*x_mid[1]+x_mid[2]*x_mid[2]);
			IJmid = dcone_index(Xm, Ym, cells,  box, &ind_mid[0], &ind_mid[1]);
			if( (ind_mid[0]>=0 && ind_mid[0]<Nx && ind_mid[1]>=0 && ind_mid[1]<Ny ) && cells[IJmid].num_inside>0)
				next = 1;
		   }
	   }
	   // finding the other particle
	   std::uniform_int_distribution<> disss2(0, cells[J].num_inside-1);
	   id = disss2(gen);
	   j = cells[J].indices_inside[ id ];
	   if( i!=j || (i==j && cells[J].num_inside>1) ){
	 	while( i == j ){
			id = disss2(gen);
			j = cells[J].indices_inside[ id ];
	   	}
	   	Ymid = increase_collision_rate(cells[IJmid].n, gas);

	   	g[0] = U1[i] - U1[j];
	   	g[1] = U2[i] - U2[j];
	   	g[2] = U3[i] - U3[j];
	   	g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

	   	dummy = 4.0*pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
	   	if( dummy - omega > 0.0)
	   		omega =  dummy;

	   	if (dummy/cells[num].omega_max> dis_r(gen) && g_dot_dir>0.0){

			U1[i] = U1[i] - g_dot_dir*dir[0];
			U2[i] = U2[i] - g_dot_dir*dir[1];
			U3[i] = U3[i] - g_dot_dir*dir[2];
			U1[j] = U1[j] + g_dot_dir*dir[0];
			U2[j] = U2[j] + g_dot_dir*dir[1];
			U3[j] = U3[j] + g_dot_dir*dir[2];
		        num_act_col++;
	   	}

	   }
	}
     }
     if(num_act_col > dummy2){
	dummy2 = num_act_col;
	num_with_max_coll = num;
     }
        if(Nc < 5 && Nc>0){
		printf("\n************            number of actual collisions = %d  at cell = %d step = %d  M_can=%e     *****************\n", num_act_col, num, (*box).step, M_cand);
		printf("Nc = %d, omega_max = %e, omega=%e",Nc, cells[num].omega_max,omega);
	}
        cells[num].omega_max = max( omega, cells[num].omega_max );
  }
  if((*box).step % (*box).every == 0)
	  printf("\n +++   cell %d had the most collisions %d with %d particles and n=%e  +++\n", num_with_max_coll, dummy2, cells[num_with_max_coll].num_inside, cells[num_with_max_coll].n);
}

void Collision_ESMC_new3(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis_r(0.0, 1.0);

  std::uniform_real_distribution<> dis_u(0.0, 1.0);
  std::uniform_real_distribution<> dis_v(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[2], x_mid[2];
  int ind[2], ind_mid[2];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double  r, Vc;

  double dummy;
  int  cand_num;
  int num, Nc, id, num2;
  int in, out;

double Vr;
int sgn;
double u, v, theta, phi;
int N_cell = (*box).N[0]*(*box).N[1]*(*box).N[2];
double *omega = (double *) malloc( N_cell * sizeof(double) );
int *sum_M = (int  *) malloc( N_cell * sizeof(int) );
int *num_act_col = (int  *) malloc( N_cell * sizeof(int) );
int m, z;
int found, M_cand;
int sw;
M_cand = 0;
int mtemp;
for(num=0; num<N_cell; num++){
  num_act_col[num] = 0;
	Nc = cells[num].num_inside;
  mtemp = floor(0.5*cells[num].omega_max*( 1.0*Nc )+dis_r(gen));
	sum_M[num] = M_cand + mtemp;
	M_cand = sum_M[num];
	omega[num]=0.0;//cells[num].omega_max;
  if(mtemp<2){
    printf("cell %d has only %d candidates\n", num, mtemp);
  }
}

for(m=0; m<M_cand; m++){
  if(sum_M[N_cell-1]>0){
	std::uniform_int_distribution<> disss1(0, sum_M[N_cell-1]);
	z = disss1(gen);
	num = 0;
	found = 0;
	while(num<N_cell && found == 0){
		if( z<=sum_M[num] ){
			I = num;
			found = 1;
			for(num2=I; num2<N_cell; num2++){
				sum_M[num2]--;
			}
		}
		num++;
	}
  if(found == 1 && cells[I].num_inside>0){
		std::uniform_int_distribution<> disss2(0, cells[I].num_inside-1);
		id = disss2(gen);
		i = cells[I].indices_inside[id];

		u = dis_u(gen);
		v = dis_v(gen);
		theta = 2.0*pi*u;
		phi = acos(2.0*v-1.0);
		dir[0] = cos(theta)*sin(phi);
		dir[1] = sin(theta)*sin(phi);
		dir[2] = cos(phi);

		// finding the cell at r + sigma*dir
		x[0] = x1[i] + (*gas).sigma*dir[0];
		x[1] = x2[i] + (*gas).sigma*dir[1];
		for(ii=0; ii<2; ii++)
			ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
		J =  ind[1];

		if(J<0)
			J = 0;
		else if(J>(*box).N[1]-1)
			J = (*box).N[1]-1;



		next = 0;
		if(J>=0 && J<(*box).N[1] )
			if( cells[J].num_inside>0 )
				next = 1;

		if(next ==1){
			// finding the other particle
		  	std::uniform_int_distribution<> disss3(0, cells[J].num_inside-1);
			  id = disss3(gen);
		  	j = cells[J].indices_inside[ id ];
        if( i == j )
				    next = 0;
			  if( next ==1){
	   			// finding the cell at r + sigma/2*dir
			  	x_mid[0] = (x1[i]+x1[j])/2.0;
			  	x_mid[1] = (x2[i]+x2[j])/2.0;
			  	//x_mid[2] = x3[i] + (*gas).sigma*dir[2]/2.0;

				  for(ii=0; ii<2; ii++)
					     ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
				  IJmid =  ind_mid[1];//*(*box).N[0] + ind_mid[0];
				  Ymid = increase_collision_rate(cells[IJmid].n, gas);

				  g[0] = U1[i] - U1[j];
				  g[1] = U2[i] - U2[j];
				  g[2] = U3[i] - U3[j];
				  g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];
          //Vr = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
	        dummy = 4.0*pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
				  if( dummy - omega[I] > 0.0)
	        	 omega[I] =  dummy;
	        if (dummy/cells[I].omega_max> dis_r(gen)) {
					       U1[i] = U1[i] - g_dot_dir*dir[0];
					       U2[i] = U2[i] - g_dot_dir*dir[1];
					       U3[i] = U3[i] - g_dot_dir*dir[2];

						     U1[j] = U1[j] + g_dot_dir*dir[0];
						     U2[j] = U2[j] + g_dot_dir*dir[1];
						     U3[j] = U3[j] + g_dot_dir*dir[2];
				         num_act_col[I]++;
				}
			}
		}
	}
 }
}
for(num=0; num<N_cell; num++){
	cells[num].omega_max = max(omega[num],cells[num].omega_max);
}

// end of function call
free(omega); free(sum_M); free(num_act_col);
}

void Collision_ESMC_new4(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  std::random_device rd;
  std::mt19937 gen(rd());


  std::uniform_real_distribution<> dis_r(0.0, 1.0);

  std::uniform_real_distribution<> dis_u(0.0, 1.0);
  std::uniform_real_distribution<> dis_v(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[2], x_mid[2];
  int ind[2], ind_mid[2];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double  r;

  double dummy;
  int num_act_col, cand_num;
  int Nc, id;
  double M_cand;

double rr;
int sgn;
double u, v, theta, phi;
int N_cell = (*box).N[0]*(*box).N[1]*(*box).N[2];
//double *omega = (double *) malloc( N_cell * sizeof(double) );
double omega;
int m;
int m_cand;
M_cand = 0.0;

//for(I=0; I<N_cell; I++){
for(I=N_cell-1; I>=0; I--){
  omega = 0.0;
  std::uniform_int_distribution<> diss(0, cells[I].num_inside-1);
  Nc = cells[I].num_inside;
  if(cells[I].num_inside >1){
    M_cand = 0.5*cells[I].omega_max*( (double) Nc );
    m_cand = floor(M_cand+dis_r(gen));
    num_act_col = 0;
    if( m_cand<2 )
      printf("m_cand in cell %d is %d\n", I, m_cand);
    for ( cand_num = 0; cand_num < m_cand; cand_num++)
    {
  	  i = cells[I].indices_inside[ diss(gen) ];
      next = 0;
      u = dis_u(gen);
    	v = dis_v(gen);
    	theta = 2.0*pi*u;
    	phi = acos(2.0*v-1.0);
    	dir[0] = cos(theta)*sin(phi);
    	dir[1] = sin(theta)*sin(phi);
    	dir[2] = cos(phi);

    	// finding the cell at r + sigma*dir
    	x[0] = x1[i] + (*gas).sigma*dir[0];
    	x[1] = x2[i] + (*gas).sigma*dir[1];
    	for(ii=0; ii<2; ii++){
    		ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
      }
    	J =  ind[1];
      if(J<0){
    		J = 0;
    	}
    	else if(J>(*box).N[1]-1){
    		J = (*box).N[1]-1;
    	}

      if( cells[J].num_inside>0 )
    		next = 1;
    	if(next ==1){
        std::uniform_int_distribution<> disss2(0, cells[J].num_inside-1);
  		  id = disss2(gen);
  	  	j = cells[J].indices_inside[ id ];
        if( i == j )
  			   next = 0;
  		  if( next ==1){
  	   		// finding the cell at r + sigma/2*dir
  		  	x_mid[0] = (x1[i]+x1[j])/2.0;// + (*gas).sigma*dir[0]/2.0;
  		  	x_mid[1] = (x2[i]+x2[j])/2.0;//x2[i] + (*gas).sigma*dir[1]/2.0;
          for(ii=0; ii<2; ii++)
    				ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
    			IJmid =  ind_mid[1];

          Ymid = increase_collision_rate(cells[IJmid].n, gas);
          g[0] = U1[i] - U1[j];
    			g[1] = U2[i] - U2[j];
    			g[2] = U3[i] - U3[j];
    			g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

    	    dummy = 4.0*pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
    			if( dummy - omega > 0.0)
    	        		omega =  dummy;
    	    if (dummy/cells[I].omega_max> dis_r(gen) && g_dot_dir>0.0) {
    				U1[i] = U1[i] - g_dot_dir*dir[0];
    				U2[i] = U2[i] - g_dot_dir*dir[1];
    				U3[i] = U3[i] - g_dot_dir*dir[2];

    				U1[j] = U1[j] + g_dot_dir*dir[0];
    				U2[j] = U2[j] + g_dot_dir*dir[1];
    				U3[j] = U3[j] + g_dot_dir*dir[2];
    			  num_act_col++;
    			}
        }
      }
    }
  }
  cells[I].omega_max = max( omega, cells[I].omega_max );
}

/// end of function
}


void Collision_ESMC_new(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::uniform_real_distribution<> dis_y(-1.0, 1.0);
  std::uniform_real_distribution<> dis_z(-1.0, 1.0);

  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::uniform_real_distribution<> dis_rr(0.0, 1.0);

  std::uniform_real_distribution<> dis_u(0.0, 1.0);
  std::uniform_real_distribution<> dis_v(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[2], x_mid[2];
  int ind[2], ind_mid[2];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double  r, Vc;

  double dummy;
  int num_act_col, cand_num;
  int num, Nc, id, num2;
  double M_cand;
  int in, out;

double rr;
int sgn;
double u, v, theta, phi;
int N_cell = (*box).N[0]*(*box).N[1]*(*box).N[2];
double *omega = (double *) malloc( N_cell * sizeof(double) );
int m;
int m_cand;
M_cand = 0.0;
for(num=0; num<N_cell; num++){
	Nc = cells[num].num_inside;
	Vc = cells[num].volume;
	M_cand += 0.5*cells[num].omega_max*( 1.0*Nc );
	omega[num]=0.0;
}
m_cand = floor(M_cand+dis_r(gen));
std::uniform_int_distribution<> disss1(0, (*gas).N-1);
num_act_col = 0;
for(m=0; m<m_cand; m++){
	i = disss1(gen);
	I = index[i];
	next = 0;
	// choose a random direction
	u = dis_u(gen);
	v = dis_v(gen);
	theta = 2.0*pi*u;
	phi = acos(2.0*v-1.0);
	dir[0] = cos(theta)*sin(phi);
	dir[1] = sin(theta)*sin(phi);
	dir[2] = cos(phi);

	// finding the cell at r + sigma*dir
	x[0] = x1[i] + (*gas).eqdist*dir[0];
	x[1] = x2[i] + (*gas).eqdist*dir[1];
	for(ii=0; ii<2; ii++)
		ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
	J =  ind[1];
	if(J<0){
		//printf("partner out from the bottom");
		out = J;
		J = 0;
	}
	else if(J>(*box).N[1]-1){
		//printf("partner out from the top");
		out = J-((*box).N[1]-1);
		J = (*box).N[1]-1;
	}
	if( cells[J].num_inside>0 )
		next = 1;
	if(next ==1){
		// finding the other particle
	  	std::uniform_int_distribution<> disss2(0, cells[J].num_inside-1);
		id = disss2(gen);
	  	j = cells[J].indices_inside[ id ];
                if( i == j )
			next = 0;
		if( next ==1){
	   		// finding the cell at r + sigma/2*dir
		  	x_mid[0] = x1[i] + (*gas).eqdist*dir[0]/2.0;
		  	x_mid[1] = x2[i] + (*gas).eqdist*dir[1]/2.0;
		  	//x_mid[2] = x3[i] + (*gas).sigma*dir[2]/2.0;

			for(ii=0; ii<2; ii++)
				ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
			IJmid =  ind_mid[1];//*(*box).N[0] + ind_mid[0];
			if( IJmid<0 )
				IJmid = 0;
			else if( IJmid>(*box).N[1]-1 )
				IJmid = (*box).N[1]-1;
			// calculating Y(ri, ri+sigma/2)
			Ymid = increase_collision_rate(cells[IJmid].n, gas);

			g[0] = U1[i] - U1[j];
			g[1] = U2[i] - U2[j];
			g[2] = U3[i] - U3[j];
			g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

	    dummy = 4.0*pi*(*gas).eqdist*(*gas).eqdist*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
			if( dummy - omega[I] > 0.0)
	        		omega[I] =  dummy;
	        	if (dummy/cells[I].omega_max> dis_r(gen) && g_dot_dir>0.0) {
				U1[i] = U1[i] - g_dot_dir*dir[0];
				U2[i] = U2[i] - g_dot_dir*dir[1];
				U3[i] = U3[i] - g_dot_dir*dir[2];

				U1[j] = U1[j] + g_dot_dir*dir[0];
				U2[j] = U2[j] + g_dot_dir*dir[1];
				U3[j] = U3[j] + g_dot_dir*dir[2];
			        num_act_col++;
			}
		}
	}
}
if (num_act_col<10){
  printf("num_act_col = %d\n", num_act_col);
}
for(num=0; num<N_cell; num++){
	cells[num].omega_max = max( omega[num], cells[num].omega_max );
}

// end of function call
}


void Collision_ESMC_new2(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::uniform_real_distribution<> dis_y(-1.0, 1.0);
  std::uniform_real_distribution<> dis_z(-1.0, 1.0);

  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::uniform_real_distribution<> dis_rr(0.0, 1.0);

  std::uniform_real_distribution<> dis_u(0.0, 1.0);
  std::uniform_real_distribution<> dis_v(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[2], x_mid[2];
  int ind[2], ind_mid[2];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double omega, r, Vc;

  double dummy;
  int num_act_col, cand_num;
  int num, Nc, id, num2;
  double M_cand;
  int in, out;

double rr;
int sgn;
double u, v, theta, phi;
  for( num = (*box).N[0]*(*box).N[1]*(*box).N[2]-1; num>=0; num--){
     omega = 0.0;
     std::uniform_int_distribution<> disss1(0, cells[num].num_inside-1);

     Nc = cells[num].num_inside;
     Vc = cells[num].volume;
     M_cand = 0.5*cells[num].omega_max*( 1.0*Nc )+dis_r(gen);
     num_act_col = 0;
     if( Nc>0){
        for( cand_num=0; cand_num < floor(M_cand); cand_num++)
        {
	    out = 0;
 	    i = cells[num].indices_inside[ disss1(gen) ];
	    I = num;
            next = 0;
           // while (next ==0) {
	  	// choose a random direction
		    u = dis_u(gen);
		    v = dis_v(gen);
		    theta = 2.0*pi*u;
		    phi = acos(2.0*v-1.0);
		    dir[0] = cos(theta)*sin(phi);
		    dir[1] = sin(theta)*sin(phi);
		    dir[2] = cos(phi);

	           // finding the cell at r + sigma*dir
	           x[0] = x1[i] + (*gas).eqdist*dir[0];
	           x[1] = x2[i] + (*gas).eqdist*dir[1];
  	              // if the other partilce is outside the domain, go for the next particle

  	  	    for(ii=0; ii<2; ii++)
			ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
		    if(ind[1]>0 && ind[1]<(*box).N[1]-1){
			     if(fabs(x[1]-(*box).delta_dim[1]*ind[1])<1e-16 || fabs(x[1]-(*box).delta_dim[1]*(ind[1]+1))<1e-16){
					rr = dis_rr(gen);
					if(rr>0.5)
						sgn = 1;
					else
						sgn = -1;
					ind[1] = floor( (x[1]+1e-15*sgn) /(*box).delta_dim[1]);
			     }
		    }

	  	    J =  ind[1];//*(*box).N[0] + ind[0];
		    if(J<0){
			printf("partner out from the bottom");
			out = J;
			J = 0;
		    }
		    else if(J>(*box).N[1]-1){
			//printf("partner out from the top");
			out = J-((*box).N[1]-1);
			J = (*box).N[1]-1;
		    }
	            if( cells[J].num_inside>0 )
		            next = 1;
            // }
	  // continue the loop if only the other cell is still in the domain

	  if(next ==1){
	  	// finding the other particle
	  	std::uniform_int_distribution<> disss2(0, cells[J].num_inside-1);
		id = disss2(gen);
	  	j = cells[J].indices_inside[ id ];
                if( i == j && out == 0){
			next = 0;
	   	}
		if( next ==1){
	   	  	// finding the cell at r + sigma/2*dir
		  	x_mid[0] = x1[i] + (*gas).eqdist*dir[0]/2.0;
		  	x_mid[1] = x2[i] + (*gas).eqdist*dir[1]/2.0;
		  	//x_mid[2] = x3[i] + (*gas).sigma*dir[2]/2.0;

		  	for(ii=0; ii<2; ii++)
				ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
			if(ind_mid[1]>0 && ind_mid[1]<(*box).N[1]-1){
			        if(fabs(x_mid[1]-(*box).delta_dim[1]*ind_mid[1])<1e-16 || fabs(x_mid[1]-(*box).delta_dim[1]*(ind_mid[1]+1))<1e-16){
					rr = dis_rr(gen);
					if(rr>0.5)
						sgn = 1;
					else
						sgn = -1;
					ind_mid[1] = floor( (x_mid[1]+1e-15*sgn) /(*box).delta_dim[1]);
			        }
			}
		  	IJmid =  ind_mid[1];//*(*box).N[0] + ind_mid[0];
			if( IJmid<0 )
				IJmid = 0;
			else if( IJmid>(*box).N[1]-1 )
				IJmid = (*box).N[1]-1;
		  	// calculating Y(ri, ri+sigma/2)
		  	Ymid = increase_collision_rate(cells[IJmid].n, gas);

			g[0] = U1[i] - U1[j];
			g[1] = U2[i] - U2[j];
			g[2] = U3[i] - U3[j];
			g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

	          	dummy = 4.0*pi*(*gas).eqdist*(*gas).eqdist*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
		  	if( dummy - omega > 0.0)
	               		omega =  dummy;

	         	if (dummy/cells[num].omega_max> dis_r(gen) && g_dot_dir>0.0) {

				U1[i] = U1[i] - g_dot_dir*dir[0];
				U2[i] = U2[i] - g_dot_dir*dir[1];
				U3[i] = U3[i] - g_dot_dir*dir[2];
				//if(out == 0){
					U1[j] = U1[j] + g_dot_dir*dir[0];
					U2[j] = U2[j] + g_dot_dir*dir[1];
					U3[j] = U3[j] + g_dot_dir*dir[2];
				//}
		                num_act_col++;
			}
		}
	  }
        }

        //if(num_act_col < 1)
	//	printf("\n************            number of actual collisions = %d          *****************\n", num_act_col);
        cells[num].omega_max = max( omega, cells[num].omega_max );
    }
  }
}


void MD_force(double *dis, double *F1,double *F2,double *F3, struct GAS *gas){
      double r = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
      *F1 += -4.0*(*gas).epsilon*6.0*pow((*gas).sigma/r,5.0)*dis[0]/r*(*gas).sigma/r/r;
      *F2 += -4.0*(*gas).epsilon*6.0*pow((*gas).sigma/r,5.0)*dis[1]/r*(*gas).sigma/r/r;
      *F3 += -4.0*(*gas).epsilon*6.0*pow((*gas).sigma/r,5.0)*dis[2]/r*(*gas).sigma/r/r;
}
void MD_potential(double *dis, double *phi, struct GAS *gas){
      double r = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
      *phi += 4.0*(*gas).epsilon*pow((*gas).sigma/r,6.0);
}
void MD_8th(double *U1, double *U2, double *U3, double *x1, double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *F1, double *F2, double *F3,
                  double *F1_old, double *F2_old, double *F3_old, double *F1_m1, double *F1_m2, double *F1_m3, double *F1_m4, double *F1_m5, double *F1_m6, double *F1_m7,
                                                                  double *F2_m1, double *F2_m2, double *F2_m3, double *F2_m4, double *F2_m5, double *F2_m6, double *F2_m7,
                                                                  double *F3_m1, double *F3_m2, double *F3_m3, double *F3_m4, double *F3_m5, double *F3_m6, double *F3_m7,
                                                                  double *x1_m1, double *x2_m1, double *x3_m1, double *x1_old, double *x2_old, double *x3_old, double *phi)
{
      int i, id,id2, num,num1, num2, neigh;
      double cut_off = 3.2*(*gas).sigma;
      double epsilon = 1e-15;
      double dis[3];
      double r;
      double rx[7];
      double ry[7];
      double rz[7];
      int ii,jj,kk, count, nn;
      double dt;
      neigh = 0.0;//floor((cut_off+epsilon)/(*box).delta_dim[1]);




      //#pragma omp parallel for private(num, num1, num2, neigh, i, id2, dis, id )

      for(id = 0; id<(*gas).N; id++){
            F1[id] = 0.0; F2[id] = 0.0; F3[id] = 0.0;
            phi[id] = 0.0;
            num = index[id];
            for(num1=num-neigh; num1<num+neigh+1; num1++){
                  if(num1 <0)
                        num2 = 0;
                  else if(num1 > (*box).N[1]-1)
                        num2 = (*box).N[1]-1;
                  else
                        num2 = num1;
                  for(i=0; i < cells[num2].num_inside; i++){
                        id2 = cells[num2].indices_inside[i];
                        if(id2 != id){
                              count = 0;
                              rx[count] = -x1[id]+x1[id2];
                              ry[count] = -x2[id]+x2[id2];
                              rz[count] = -x3[id]+x3[id2];
                              count ++;
                              for(ii=-1;ii<2;ii++){
                                    if(ii!=0){
                                          rx[count] = -x1[id]+x1[id2] + ii*(*box).delta_dim[0];
                                          ry[count] = -x2[id]+x2[id2];
                                          rz[count] = -x3[id]+x3[id2];
                                          count ++;
                                    }
                              }
                              for(jj=-1;jj<2;jj++){
                                    if(jj!=0){
                                          rx[count] = -x1[id]+x1[id2];
                                          ry[count] = -x2[id]+x2[id2] + jj*(*box).delta_dim[1];
                                          rz[count] = -x3[id]+x3[id2];
                                          count ++;
                                    }
                              }
                              for(kk=-1;kk<2;kk++){
                                    if(kk!=0){
                                          rx[count] = -x1[id]+x1[id2];
                                          ry[count] = -x2[id]+x2[id2];
                                          rz[count] = -x3[id]+x3[id2] + kk*(*box).delta_dim[2];
                                          count ++;
                                    }
                              }
                              /*
                              for(ii=-1;ii<2;ii++){
                                    for(jj=-1;jj<2;jj++){
                                          for(kk=-1;kk<2;kk++){
                                                rx[count] = -x1[id]+x1[id2] + ii*(*box).delta_dim[0];
                                                ry[count] = -x2[id]+x2[id2] + jj*(*box).delta_dim[1];
                                                rz[count] = -x3[id]+x3[id2] + kk*(*box).delta_dim[2];
                                                count ++;
                                          }
                                    }
                              }
                              */
                              min_rij(7, rx, ry, rz, &nn, &r);
                              dis[0] = rx[nn];
                              dis[1] = ry[nn];
                              dis[2] = rz[nn];
//                              dis[0] = -x1[id]+x1[id2];
//                              dis[1] = -x2[id]+x2[id2];//+(num1-num2)*(*box).delta_dim[1]);
//                              dis[2] = -x3[id]+x3[id2];
//                              r = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
                              if(r-(*gas).sigma>1e-13){
                                    MD_force(dis, &F1[id],&F2[id],&F3[id],gas);
                                    MD_potential(dis, &phi[id], gas);
                              }
                              if(F1[id] != F1[id]){
                                    printf("F1[%d] = %e in step = %d\n",id, F1[id] ,(*box).step);
                                    printf("dis = %e  %e  %e\n",dis[0], dis[1], dis[2] );
                                    printf("id = %d and id2 = %d\n", id, id2);
                                    getchar();
                              }
                              if(F2[id] != F2[id]){
                                    printf("F2[%d] = %e in step = %d\n",id, F2[id] ,(*box).step);
                                    getchar();
                              }
                              if(F3[id] != F3[id]){
                                    printf("F3[%d] = %e in step = %d\n",id, F3[id] ,(*box).step);
                                    getchar();
                              }
                        }
                  }
            }
      }
      /*
      if((*box).step == 1){
         for(id = 0; id<(*gas).N; id++){
            F1_m1[id] = F1[id];
            F1_m2[id] = F1[id];
            F1_m3[id] = F1[id];
            F1_m4[id] = F1[id];
            F1_m5[id] = F1[id];
            F1_m6[id] = F1[id];
            F1_m7[id] = F1[id];

            F2_m1[id] = F2[id];
            F2_m2[id] = F2[id];
            F2_m3[id] = F2[id];
            F2_m4[id] = F2[id];
            F2_m5[id] = F2[id];
            F2_m6[id] = F2[id];
            F2_m7[id] = F2[id];

            F3_m1[id] = F3[id];
            F3_m2[id] = F3[id];
            F3_m3[id] = F3[id];
            F3_m4[id] = F3[id];
            F3_m5[id] = F3[id];
            F3_m6[id] = F3[id];
            F3_m7[id] = F3[id];

         }
      }
      */
      //#pragma omp parallel for private(id)


if((*box).step ==1){

      dt = (*gas).delta_t;
      //#pragma omp parallel for private(id)
      for(id = 0; id<(*gas).N; id++){
            U1[id] = U1[id] + dt*(F1[id]+F1_old[id])/2.0/(*gas).m;
            U2[id] = U2[id] + dt*(F2[id]+F2_old[id])/2.0/(*gas).m;
            U3[id] = U3[id] + dt*(F3[id]+F3_old[id])/2.0/(*gas).m;
      }

      dt = (*gas).delta_t;
      if( (*box).direction[0] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
               x1[id] =  x1[id] + dt*U1[id];
      }
      if( (*box).direction[1] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  x2[id] + dt*U2[id];
       }
       if( (*box).direction[2] == 1 ){
            // #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  x3[id] + dt*U3[id];
       }
       //BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
       BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);

}
else{
      dt = (*gas).delta_t;
      if( (*box).direction[0] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
               x1[id] =  2.0*x1[id]-x1_m1[id] + dt*dt/(60480.0)*( 88324.0*F1_old[id]-121797.0*F1_m1[id]+245598.0*F1_m2[id]-300227.0*F1_m3[id]+236568.0*F1_m4[id]-117051.0*F1_m5[id]+33190.0*F1_m6[id]-4125.0*F1_m7[id] )/(*gas).m;
      }
      if( (*box).direction[1] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  2.0*x2[id]-x2_m1[id] + dt*dt/(60480.0)*( 88324.0*F2_old[id]-121797.0*F2_m1[id]+245598.0*F2_m2[id]-300227.0*F2_m3[id]+236568.0*F2_m4[id]-117051.0*F2_m5[id]+33190.0*F2_m6[id]-4125.0*F2_m7[id] )/(*gas).m;
        }
        if( (*box).direction[2] == 1 ){
            // #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  2.0*x3[id]-x3_m1[id] + dt*dt/(60480.0)*( 88324.0*F3_old[id]-121797.0*F3_m1[id]+245598.0*F3_m2[id]-300227.0*F3_m3[id]+236568.0*F3_m4[id]-117051.0*F3_m5[id]+33190.0*F3_m6[id]-4125.0*F3_m7[id] )/(*gas).m;
        }

        for(id = 0; id<(*gas).N; id++){
             U1[id] = 1.0/(2.0*dt)*(x1[id] - x1_m1[id]) + dt/3.0*(F1[id]+2.0*F1_old[id])/(*gas).m;
             U2[id] = 1.0/(2.0*dt)*(x2[id] - x2_m1[id]) + dt/3.0*(F2[id]+2.0*F2_old[id])/(*gas).m;
             U3[id] = 1.0/(2.0*dt)*(x3[id] - x3_m1[id]) + dt/3.0*(F3[id]+2.0*F3_old[id])/(*gas).m;
        }
}

        for(id=0; id<(*gas).N; id++){
             F1_old[id] = F1[id];
             F2_old[id] = F2[id];
             F3_old[id] = F3[id];

             F1_m1[id] = F1_old[id];
             F1_m2[id] = F1_m1[id];
             F1_m3[id] = F1_m2[id];
             F1_m4[id] = F1_m3[id];
             F1_m5[id] = F1_m4[id];
             F1_m6[id] = F1_m5[id];
             F1_m7[id] = F1_m6[id];

             F2_m1[id] = F2_old[id];
             F2_m2[id] = F2_m1[id];
             F2_m3[id] = F2_m2[id];
             F2_m4[id] = F2_m3[id];
             F2_m5[id] = F2_m4[id];
             F2_m6[id] = F2_m5[id];
             F2_m7[id] = F2_m6[id];

             F3_m1[id] = F3_old[id];
             F3_m2[id] = F3_m1[id];
             F3_m3[id] = F3_m2[id];
             F3_m4[id] = F3_m3[id];
             F3_m5[id] = F3_m4[id];
             F3_m6[id] = F3_m5[id];
             F3_m7[id] = F3_m6[id];

             x1_m1[id] = x1_old[id];
             x2_m1[id] = x2_old[id];
             x3_m1[id] = x3_old[id];
        }

        //BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
        BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);
}


void MD_new(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index,
         double *F1, double *F2, double *F3){


      int i, j, i2, id,id2, num,num1, num2, neigh;
      double cut_off = 3.0*(*gas).sigma;
      double epsilon = 1e-15;
      double dis[3];
      double r;
      double rx;
      double ry;
      double rz;
      int ii,jj,kk, count, nn;
      double dt;
      double neighU, neighL;

      dt = (*gas).delta_t/2.0;
      if( (*box).direction[0] == 1 ){
          for(id=0; id<(*gas).N; id++)
         // #pragma omp parallel for
               x1[id] =  x1[id] + dt*U1[id];
       }
       if( (*box).direction[1] == 1 ){
      //       #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  x2[id] + dt*U2[id];
        }
        if( (*box).direction[2] == 1 ){
      //       #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  x3[id] + dt*U3[id];
        }
      //BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);

	if( (*box).problem == evaporation || (*box).problem == vacuum){
	        BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);
	}
	else if( (*box).problem == inverted ){
		BC_specular_periodic_reflection(U1, U2, U3, x1, x2,x3, gas, box, cells);
	}
	// should be here cell_update(x1, x2, U1, U2, U3, gas, cells, box, index);
      //#pragma omp parallel for private(num, num1, num2, neigh, i, id2, dis, id )

//#pragma omp parallel private(neighU, neighL, j, num2, i2, id, ii, kk, rx, ry, rz, r)
	for(i=0; i<(*gas).N; i++){
		F1[i] = 0.0; F2[i] = 0.0; F3[i] = 0.0;
		neighU =  ceil( (x2[i]+cut_off)/(*box).delta_dim[1] );
		neighL = floor( (x2[i]-cut_off)/(*box).delta_dim[1] );
		for(j=neighL; j<= neighU; j++){
			if(j<0)
				num2 = 0;
			else if(j>(*box).N[1]-1)
				num2 = (*box).N[1]-1;
			else
				num2 = j;
			for(i2=0; i2<cells[num2].num_inside; i2++){
				id2 = cells[num2].indices_inside[i2];
//				for(ii=0; ii<1; ii++){
//					for(kk=0; kk<1; kk++){
				for(ii=-1; ii<=1; ii++){
					for(kk=-1; kk<=1; kk++){
						rx = x1[id2]-x1[i]+ii*(*box).delta_dim[0];
						ry = x2[id2]-x2[i]+(j-num2)*(*box).delta_dim[1];
						rz = x3[id2]-x3[i]+kk*(*box).delta_dim[2];
						r = sqrt( rx*rx+ry*ry+rz*rz );
						if(r<cut_off && r>0.1*(*gas).sigma){
							F1[i] += -4.0*(*gas).epsilon*6.0*( pow((*gas).sigma/r,5.0)-2.0*pow((*gas).sigma/r,11.0) )*rx/r*(*gas).sigma/r/r;
							F2[i] += -4.0*(*gas).epsilon*6.0*( pow((*gas).sigma/r,5.0)-2.0*pow((*gas).sigma/r,11.0) )*ry/r*(*gas).sigma/r/r;
							F3[i] += -4.0*(*gas).epsilon*6.0*( pow((*gas).sigma/r,5.0)-2.0*pow((*gas).sigma/r,11.0) )*rz/r*(*gas).sigma/r/r;
						}
					}
				}
			}
		}
	}


      dt = (*gas).delta_t;
      //#pragma omp parallel for private(id)
      for(id = 0; id<(*gas).N; id++){
            U1[id] = U1[id] - dt*(F1[id])/(*gas).m;
            U2[id] = U2[id] - dt*(F2[id])/(*gas).m;
            U3[id] = U3[id] - dt*(F3[id])/(*gas).m;
      }


//      free(F1); free(F2); free(F3);
//      free(F1_old); free(F2_old); free(F3_old);
      dt = (*gas).delta_t/2.0;
      if( (*box).direction[0] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
               x1[id] =  x1[id] + dt*U1[id];
      }
      if( (*box).direction[1] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  x2[id] + dt*U2[id];
        }
        if( (*box).direction[2] == 1 ){
            // #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  x3[id] + dt*U3[id];
        }
	if( (*box).problem == evaporation || (*box).problem == vacuum){
	        BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);
	}
	else if( (*box).problem== inverted){
		BC_specular_periodic_reflection(U1, U2, U3, x1, x2,x3, gas, box, cells);
	}
	// should be here cell_update(x1, x2, U1, U2, U3, gas, cells, box, index);
}
void MD(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index,
         double *F1, double *F2, double *F3, double *F1_old, double *F2_old, double *F3_old){


      int i, id,id2, num,num1, num2, neigh;
      double cut_off = 2.0*(*gas).sigma;
      double epsilon = 1e-15;
      double dis[3];
      double r;
      double rx[27];
      double ry[27];
      double rz[27];
      int ii,jj,kk, count, nn;
      double dt;
      neigh = 0.0;//floor((cut_off+epsilon)/(*box).delta_dim[1]);


      dt = (*gas).delta_t/2.0;
      if( (*box).direction[0] == 1 ){
          for(id=0; id<(*gas).N; id++)
         // #pragma omp parallel for
               x1[id] =  x1[id] + dt*U1[id];//dt*(U1[id]+0.5*F1_old[id]/(*gas).m*dt);
       }
       if( (*box).direction[1] == 1 ){
      //       #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  x2[id] + dt*U2[id];//dt*(U2[id]+0.5*F2_old[id]/(*gas).m*dt);
        }
        if( (*box).direction[2] == 1 ){
      //       #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  x3[id] + dt*U3[id];//dt*(U3[id]+0.5*F3_old[id]/(*gas).m*dt);
        }
      //BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
      BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);

      //#pragma omp parallel for private(num, num1, num2, neigh, i, id2, dis, id )

      for(id = 0; id<(*gas).N; id++){
            F1[id] = 0.0; F2[id] = 0.0; F3[id] = 0.0;
            num = index[id];
            for(num1=num-neigh; num1<num+neigh+1; num1++){
                  if(num1 <0)
                        num2 = 0;
                  else if(num1 > (*box).N[1]-1)
                        num2 = (*box).N[1]-1;
                  else
                        num2 = num1;
                  for(i=0; i < cells[num2].num_inside; i++){
                        id2 = cells[num2].indices_inside[i];
                        if(id2 != id){
                              count = 0;
                              for(ii=-1;ii<2;ii++){
                                    for(jj=-1;jj<2;jj++){
                                          for(kk=-1;kk<2;kk++){
                                                rx[count] = -x1[id]+x1[id2] + ii*(*box).delta_dim[0];
                                                ry[count] = -x2[id]+x2[id2] + jj*(*box).delta_dim[1];
                                                rz[count] = -x3[id]+x3[id2] + kk*(*box).delta_dim[2];
                                                count ++;
                                          }
                                    }
                              }
                              min_rij(27, rx, ry, rz, &nn, &r);
                              dis[0] = rx[nn];
                              dis[1] = ry[nn];
                              dis[2] = rz[nn];
//                              dis[0] = -x1[id]+x1[id2];
//                              dis[1] = -x2[id]+x2[id2];//+(num1-num2)*(*box).delta_dim[1]);
//                              dis[2] = -x3[id]+x3[id2];
//                              r = sqrt(dis[0]*dis[0] + dis[1]*dis[1] + dis[2]*dis[2]);
                              if(r-(*gas).sigma>1e-14)
                                    MD_force(dis, &F1[id],&F2[id],&F3[id],gas);
                              if(F1[id] != F1[id]){
                                    printf("F1[%d] = %e in step = %d\n",id, F1[id] ,(*box).step);
                                    printf("dis = %e  %e  %e\n",dis[0], dis[1], dis[2] );
                                    printf("id = %d and id2 = %d\n", id, id2);
                                    getchar();
                              }
                              if(F2[id] != F2[id]){
                                    printf("F2[%d] = %e in step = %d\n",id, F2[id] ,(*box).step);
                                    getchar();
                              }
                              if(F3[id] != F3[id]){
                                    printf("F3[%d] = %e in step = %d\n",id, F3[id] ,(*box).step);
                                    getchar();
                              }
                        }
                  }
            }
      }
      dt = (*gas).delta_t;
      //#pragma omp parallel for private(id)
      for(id = 0; id<(*gas).N; id++){
            U1[id] = U1[id] + dt*(F1[id]+F1_old[id])/2.0/(*gas).m;
            U2[id] = U2[id] + dt*(F2[id]+F2_old[id])/2.0/(*gas).m;
            U3[id] = U3[id] + dt*(F3[id]+F3_old[id])/2.0/(*gas).m;
      }
      for(id=0; id<(*gas).N; id++){
            F1_old[id] = F1[id];
            F2_old[id] = F2[id];
            F3_old[id] = F3[id];
      }


//      free(F1); free(F2); free(F3);
//      free(F1_old); free(F2_old); free(F3_old);
      dt = (*gas).delta_t/2.0;
      if( (*box).direction[0] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
               x1[id] =  x1[id] + dt*U1[id];
      }
      if( (*box).direction[1] == 1 ){
            //#pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                x2[id] =  x2[id] + dt*U2[id];
        }
        if( (*box).direction[2] == 1 ){
            // #pragma omp parallel for
          for(id=0; id<(*gas).N; id++)
                 x3[id] =  x3[id] + dt*U3[id];
        }
        //BC_specular_reflection(U1, U2, U3, x1, x2,x3, gas, box);
        BC_periodic(U1, U2, U3, x1, x2,x3, gas, box, cells);

}
void Collision_ESMC_1D(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *ESMC, int *index){
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::mt19937 gen2(rd());
  std::mt19937 gen3(rd());
  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::uniform_real_distribution<> dis_y(-1.0, 1.0);
  std::uniform_real_distribution<> dis_z(-1.0, 1.0);

  std::uniform_real_distribution<> dis_r(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[3], x_mid[3], ind[3], ind_mid[3];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double omega, r, Vc, g_magn;

  (*ESMC).pi = acos(-1.0);
  double dummy;
  double check;
  int num_act_col;
  int sw;
  //printf("omega_max = %lf\n", cells[1].omega_max);

  for( (*ESMC).num = 0; (*ESMC).num<(*box).N[0]*(*box).N[1]*(*box).N[2]; (*ESMC).num++){
	sw = 0;
	//if(cells[(*ESMC).num].cell_center[1]<4.0*(*gas).sigma || cells[(*ESMC).num].cell_center[1]> (*box).Len[1]-4.0*(*gas).sigma)
	//	sw = 1;

     omega = 0.0;
     std::uniform_int_distribution<> diss(0, cells[(*ESMC).num].num_inside-1);

     (*ESMC).Nc = cells[(*ESMC).num].num_inside;
    // if(cells[(*ESMC).num].num_inside <3){
//	printf("in cell %d, number of particles are %d\n", (*ESMC).num, cells[(*ESMC).num].num_inside);
  //   }
if(cells[(*ESMC).num].num_inside >1){
     (*ESMC).M_cand = 0.5*cells[(*ESMC).num].omega_max*( (double)(*ESMC).Nc );
     num_act_col = 0;
     if( floor( (*ESMC).M_cand )<2 )
        printf("M_cand in cell %d is %f\n", (*ESMC).num, floor((*ESMC).M_cand) );
     for ( (*ESMC).cand_num = 0; (*ESMC).cand_num < floor((*ESMC).M_cand+dis_x(gen)); (*ESMC).cand_num++)
     {
	  (*ESMC).pair[0] = cells[(*ESMC).num].indices_inside[ diss(gen) ];
 	  i = (*ESMC).pair[0];
	  I = (*ESMC).num;
        next = 0;
       // while (next ==0) {
	//if(x2[i]>(*box).ghost && x2[i]<(*box).Len[1]-( (*box).ghost)){
//	if(x2[i]>(*gas).sigma && x2[i]<(*box).Len[1]-( (*gas).sigma)){
	//if(x2[i]>(*gas).sigma && x2[i]<(*box).Len[1]-( (*gas).sigma)){
	  	// choose a random direction
	           dir[0] = dis_x(gen);
	           dir[1] = dis_y(gen);
	           dir[2] = dis_z(gen);
	           size = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	           for(ii=0; ii<3; ii++)
		          dir[ii] = dir[ii]/size;

	           // finding the cell at r + sigma*dir
	           x[0] = x1[i] + (*gas).sigma*dir[0];
	           x[1] = x2[i] + (*gas).sigma*dir[1];
	           x[2] = x3[i] + (*gas).sigma*dir[2];
                 //if( x[0] - (*box).Len[0] < 0.0 && x[0] > 0.0 && x[1] - (*box).Len[1] < 0.0 && x[1] > 0.0 && x[2] - (*box).Len[2] < 0.0 && x[2] > 0.0){
		  	for(ii=0; ii<3; ii++)
				ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
		  	J = ind[2]*(*box).N[0]*(*box).N[1] + ind[1]*(*box).N[0] + ind[0];
			//J = I;
			//if(sw==1)
			//	J = I;
      if(J>=0 && J<(*box).N[1] ){
  			if( cells[J].num_inside>0 ){
  				next = 1;
  				//if( out == 0){
  					//for(num2=J; num2<N_cell; num2++){
  					//	sum_M[num2]--;
  					//}
  				//}
  			}
  		}
		 //}
	//}
      //  }
	  // continue the loop if only the other cell is still in the domain

	  if(next ==1){

	  	// finding the other particle
	  	std::uniform_int_distribution<> disss(0, cells[J].num_inside-1);
		int id = disss(gen2);
	  	j = cells[J].indices_inside[ id ];


           if( i == j ){
                 j = cells[(*ESMC).num].indices_inside[ disss(gen2) ];
           }
	   if( i == j || cells[J].num_inside < 1){
		next =0;
	   }
/*
           // correct the direction of x1+sigma*dir
           dir[0] = x1[j]-x1[i];
           dir[1] = x2[j]-x2[i];
           dir[2] = x3[j]-x3[i];
           size = sqrt(dir[0]*dir[0] + dir[1]*dir[1]+dir[2]*dir[2]);
           for(ii=0; ii<3; ii++)
                 dir[ii] = dir[ii]/size;
*/
	   if(next == 1){
   	  	// finding the cell at r + sigma/2*dir
	  	x_mid[0] = (x1[i] + x1[j])/2.0;//(*gas).sigma*dir[0]/2.0;
	  	x_mid[1] = (x2[i] + x2[j])/2.0;//x2[i] + (*gas).sigma*dir[1]/2.0;
	  	x_mid[2] = (x3[i] + x3[j])/2.0;//x3[i] + (*gas).sigma*dir[2]/2.0;

	  	for(ii=0; ii<3; ii++)
			ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
	  	IJmid = ind_mid[2]*(*box).N[0]*(*box).N[1] + ind_mid[1]*(*box).N[0] + ind_mid[0];

	  	// calculating Y(ri, ri+sigma/2)
	  	Ymid = increase_collision_rate(cells[IJmid].n, gas);


		g[0] = U1[i] - U1[j];
		g[1] = U2[i] - U2[j];
		g[2] = U3[i] - U3[j];
		g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

                (*ESMC).Vr = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
              //sigmat = pow( (*DSMC).Vr,(2.0-2.0*(*gas).visp) )*pow( (*gas).crref,2.0*(*gas).visp-1.0 );
            //if( sigmat - crm > 0.0){
            //       crm =  sigmat;
            //}
                  dummy = 4.0*(*ESMC).pi*(*gas).sigma*(*gas).sigma*(*ESMC).Vr*Ymid*cells[J].n*(*gas).delta_t;
          	//dummy = 4.0*(*ESMC).pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
	  	if( dummy - omega > 0.0)
               		omega =  dummy;
                //check = 4.0*(*ESMC).pi*(*gas).sigma*(*gas).sigma*cells[J].n*Ymid*g_dot_dir*(*gas).delta_t;
         	if (dummy/cells[(*ESMC).num].omega_max> dis_r(gen3) && g_dot_dir>0.0) {

			U1[i] = U1[i] - g_dot_dir*dir[0];
			U2[i] = U2[i] - g_dot_dir*dir[1];
			U3[i] = U3[i] - g_dot_dir*dir[2];

			U1[j] = U1[j] + g_dot_dir*dir[0];
			U2[j] = U2[j] + g_dot_dir*dir[1];
			U3[j] = U3[j] + g_dot_dir*dir[2];
            		num_act_col++;
		}
	     }
	  }

     }

   //if(num_act_col < 1)
	  // printf("\n************            number of actual collisions = %d          *****************\n", num_act_col);
   cells[(*ESMC).num].omega_max = max( omega, cells[(*ESMC).num].omega_max );
}
  }

/*

  for(i=0; i<(*gas).N; i++){
	// don't go next loop before doing the current one
	next = 0;

	I = index[i];

	// choose a random direction
	dir[0] = dis_x(gen);
	dir[1] = dis_y(gen);
	dir[2] = dis_z(gen);
	size = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	for(ii=0; ii<3; ii++)
		dir[ii] = dir[ii]/size;

	// finding the cell at r + sigma*dir
	x[0] = x1[i] + (*gas).sigma*dir[0];
	x[1] = x2[i] + (*gas).sigma*dir[1];
	x[2] = x3[i] + (*gas).sigma*dir[2];

	// if the other partilce is outside the domain, go for the next particle
	if( x[1] - (*box).Len[1] > 0.0 || x[1] < 0.0)
		next = 1;
	// continue the loop if only the other cell is still in the domain
	if(next ==0){
		for(ii=0; ii<3; ii++)
			ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
		J = ind[2]*(*box).N[0]*(*box).N[1] + ind[1]*(*box).N[0] + ind[0];

		// finding the other particle
		std::uniform_int_distribution<> disss(0, cells[J].num_inside-1);
		j = cells[J].indices_inside[ disss(gen) ];

		// finding the cell at r + sigma/2*dir
		x_mid[0] = x1[i] + (*gas).sigma*dir[0]/2.0;
		x_mid[1] = x2[i] + (*gas).sigma*dir[1]/2.0;
		x_mid[2] = x3[i] + (*gas).sigma*dir[2]/2.0;

		for(ii=0; ii<3; ii++)
			ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
		IJmid = ind_mid[2]*(*box).N[0]*(*box).N[1] + ind_mid[1]*(*box).N[0] + ind_mid[0];

		// calculating Y(ri, ri+sigma/2)
		Ymid = increase_collision_rate(cells[IJmid].n, gas);

		g[0] = U1[i] - U1[j];
		g[1] = U2[i] - U2[j];
		g[2] = U3[i] - U3[j];
		g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];
		g_magn = sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
		Vc = ( cells[I].dim[1]-cells[I].dim[0] )*( cells[I].dim[3]-cells[I].dim[2] )*( cells[I].dim[5]-cells[I].dim[4] );

          	dummy = 4.0*(*ESMC).pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
	  	if( dummy - cells[I].omega_max > 0.0)
               		cells[I].omega_max =  dummy;

		// check whether the collision happens or not
		r = dis_r(gen);
		if( dummy/cells[I].omega_max> r && g_dot_dir > 0.0){
			U1[i] = U1[i] - g_dot_dir*dir[0];
			U2[i] = U2[i] - g_dot_dir*dir[1];
			U3[i] = U3[i] - g_dot_dir*dir[2];

			U1[j] = U1[j] + g_dot_dir*dir[0];
			U2[j] = U2[j] + g_dot_dir*dir[1];
			U3[j] = U3[j] + g_dot_dir*dir[2];
		}
	}
  }
*/

}


void Collision_ESMC(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *ESMC, int *index){
 // omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::uniform_real_distribution<> dis_y(-1.0, 1.0);
  std::uniform_real_distribution<> dis_z(-1.0, 1.0);

  std::uniform_real_distribution<> dis_r(0.0, 1.0);

  double pi = acos(-1.0);
  int i,j, I, J, IJmid;
  int ii,jj;
  double dir[3], x[2], x_mid[2], ind[2], ind_mid[2];
  double g[3], Ymid;
  double size, g_dot_dir;
  int next = 0;
  double omega, r, Vc, g_magn;

  (*ESMC).pi = acos(-1.0);
  double dummy;
  int num_act_col;

  //printf("omega_max = %lf\n", cells[1].omega_max);

  for( (*ESMC).num = 0; (*ESMC).num<(*box).N[0]*(*box).N[1]*(*box).N[2]; (*ESMC).num++){
     omega = 0.0;
     std::uniform_int_distribution<> disss(0, cells[(*ESMC).num].num_inside-1);

     (*ESMC).Nc = cells[(*ESMC).num].num_inside;
     (*ESMC).Vc = ( cells[(*ESMC).num].dim[1]-cells[(*ESMC).num].dim[0] )*( cells[(*ESMC).num].dim[3]-cells[(*ESMC).num].dim[2] )*( cells[(*ESMC).num].dim[5]-cells[(*ESMC).num].dim[4] );
     //(*ESMC).M_cand = (*ESMC).pi*(*gas).sigma*(*gas).sigma*( (double ) (*ESMC).Nc*cells[(*ESMC).num].n )*cells[(*ESMC).num].crm*(*gas).delta_t/2.0;
     //(*ESMC).Y = increase_collision_rate(cells[(*ESMC).num].n, gas);

//     (*ESMC).M_cand = (*ESMC).Y*(*ESMC).M_cand;
     (*ESMC).M_cand = 0.5*cells[(*ESMC).num].omega_max*( (double)(*ESMC).Nc );
     num_act_col = 0;

     for ( (*ESMC).cand_num = 0; (*ESMC).cand_num < floor((*ESMC).M_cand+dis_x(gen)); (*ESMC).cand_num++)
     {
	  (*ESMC).pair[0] = cells[(*ESMC).num].indices_inside[ disss(gen) ];
 	  i = (*ESMC).pair[0];
	  I = (*ESMC).num;
        next = 0;
        while (next ==0) {
	  	// choose a random direction
	           dir[0] = dis_x(gen);
	           dir[1] = dis_y(gen);
	           dir[2] = dis_z(gen);
	           size = sqrt(dir[0]*dir[0] + dir[1]*dir[1]+dir[2]*dir[2]);
	           for(ii=0; ii<3; ii++)
		          dir[ii] = dir[ii]/size;

	           // finding the cell at r + sigma*dir
	           x[0] = x1[i] + (*gas).sigma*dir[0];
	           x[1] = x2[i] + (*gas).sigma*dir[1];
  	              // if the other partilce is outside the domain, go for the next particle
	           if( x[0] - (*box).Len[0] < 0.0 && x[0] > 0.0 && x[1] - (*box).Len[1] < 0.0 && x[1] > 0.0 )
		            next = 1;
         }
	  // continue the loop if only the other cell is still in the domain

	  if(next ==1){

	  	for(ii=0; ii<2; ii++)
			ind[ii] = floor(x[ii]/(*box).delta_dim[ii]);
	  	J =  ind[1]*(*box).N[0] + ind[0];

	  	// finding the other particle
	  	std::uniform_int_distribution<> disss(0, cells[J].num_inside-1);
		int id = disss(gen);
	  	j = cells[J].indices_inside[ id ];
            (*ESMC).pair[1] = j;
            while( (*ESMC).pair[0] == (*ESMC).pair[1] ){
			(*ESMC).pair[1] = cells[(*ESMC).num].indices_inside[ disss(gen) ];
			j = (*ESMC).pair[1];
	   	}

		// correct the direction of x1+sigma*dir
	//	dir[0] = x1[j]-x1[i];
	//	dir[1] = x2[j]-x2[i];
      //          size = sqrt(dir[0]*dir[0] + dir[1]*dir[1]+dir[2]*dir[2]);
      //          for(ii=0; ii<3; ii++)
      //               dir[ii] = dir[ii]/size;

   	  	// finding the cell at r + sigma/2*dir
	  	x_mid[0] = x1[i] + (*gas).sigma*dir[0]/2.0;
	  	x_mid[1] = x2[i] + (*gas).sigma*dir[1]/2.0;
	  	//x_mid[2] = x3[i] + (*gas).sigma*dir[2]/2.0;

	  	for(ii=0; ii<2; ii++)
			ind_mid[ii] = floor(x_mid[ii]/(*box).delta_dim[ii]);
	  	IJmid =  ind_mid[1]*(*box).N[0] + ind_mid[0];

	  	// calculating Y(ri, ri+sigma/2)
	  	Ymid = increase_collision_rate(cells[IJmid].n, gas);

		g[0] = U1[i] - U1[j];
		g[1] = U2[i] - U2[j];
		g[2] = U3[i] - U3[j];
		g_dot_dir = g[0]*dir[0] + g[1]*dir[1] + g[2]*dir[2];

          	dummy = 4.0*(*ESMC).pi*(*gas).sigma*(*gas).sigma*g_dot_dir*Ymid*cells[J].n*(*gas).delta_t;
	  	if( dummy - omega > 0.0)
               		omega =  dummy;

         	if (dummy/cells[(*ESMC).num].omega_max> dis_r(gen) && g_dot_dir>0.0) {
            //if (dummy/cells[(*ESMC).num].omega_max> dis_r(gen)) {

			U1[i] = U1[i] - g_dot_dir*dir[0];
			U2[i] = U2[i] - g_dot_dir*dir[1];
			U3[i] = U3[i] - g_dot_dir*dir[2];

			U1[j] = U1[j] + g_dot_dir*dir[0];
			U2[j] = U2[j] + g_dot_dir*dir[1];
			U3[j] = U3[j] + g_dot_dir*dir[2];
                  num_act_col++;
		}
	  }

     }

   if(num_act_col < 1)
	printf("\n************            number of actual collisions = %d          *****************\n", num_act_col);
   cells[(*ESMC).num].omega_max = max( omega, cells[(*ESMC).num].omega_max );
  }

}



void Collision_DSMC_HS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC){
 // omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());

  (*DSMC).pi = acos(-1.0);
  double sigmat;
  double crm;
  int num_act_col;
  for( (*DSMC).num = 0; (*DSMC).num<(*box).N[0]*(*box).N[1]*(*box).N[2]; (*DSMC).num++){
     crm = 0.0;
     std::uniform_real_distribution<> dis_x(0, 1);
     std::uniform_int_distribution<> disss(0, cells[(*DSMC).num].num_inside-1);
	 //(*DSMC).Vr_max = sqrt((*DSMC).pi* (*gas).kb*cells[(*DSMC).num].T/(*gas).m);
         //(*DSMC).Vr_max = cr_max [ (*DSMC).num ];
     (*DSMC).Nc = cells[(*DSMC).num].num_inside;
     (*DSMC).Vc = ( cells[(*DSMC).num].dim[1]-cells[(*DSMC).num].dim[0] )*( cells[(*DSMC).num].dim[3]-cells[(*DSMC).num].dim[2] )*( cells[(*DSMC).num].dim[5]-cells[(*DSMC).num].dim[4] );
	 //(*DSMC).M_cand = (*gas).Fn*(*DSMC).Nc*(*DSMC).Nc*(*DSMC).pi*(*gas).sigma*(*gas).sigma*(*DSMC).Vr_max*(*gas).delta_t/(2.0*(*DSMC).Vc);
         //(*DSMC).M_cand = ( (double ) (*DSMC).Nc*cells[(*DSMC).num].n )*(*DSMC).pi*(*gas).sigma*(*gas).sigma*(*DSMC).Vr_max*(*gas).delta_t/2.0;
         //(*DSMC).M_cand = ( (double ) (*DSMC).Nc*cells[(*DSMC).num].n )*(*DSMC).pi*(*gas).sigma*(*gas).sigma*thermal*pow( (*gas).crref/thermal,2.0*(*gas).visp-1.0 )*(*gas).delta_t/2.0;
         //(*DSMC).M_cand = ( (double ) (*DSMC).Nc*cells[(*DSMC).num].n )*cells[ (*DSMC).num ].crm*(*gas).delta_t/2.0;
     (*DSMC).M_cand = (*DSMC).pi*(*gas).sigma*(*gas).sigma*( (double ) (*DSMC).Nc*cells[(*DSMC).num].n )*cells[(*DSMC).num].crm*(*gas).delta_t/2.0;
     num_act_col = 0;
     for ( (*DSMC).cand_num = 0; (*DSMC).cand_num < floor((*DSMC).M_cand+dis_x(gen)); (*DSMC).cand_num++)
     {
	  (*DSMC).pair[0] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
	  (*DSMC).pair[1] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
	  while( (*DSMC).pair[0] == (*DSMC).pair[1] ){
		(*DSMC).pair[1] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
	   }
	  (*DSMC).Vr = sqrt( 	  (U1[(*DSMC).pair[0]]-U1[(*DSMC).pair[1]])*(U1[(*DSMC).pair[0]]-U1[(*DSMC).pair[1]])
				+ (U2[(*DSMC).pair[0]]-U2[(*DSMC).pair[1]])*(U2[(*DSMC).pair[0]]-U2[(*DSMC).pair[1]])
				+ (U3[(*DSMC).pair[0]]-U3[(*DSMC).pair[1]])*(U3[(*DSMC).pair[0]]-U3[(*DSMC).pair[1]]) );

//            sigmat = (*DSMC).pi*(*gas).sigma*(*gas).sigma*pow( (*DSMC).Vr,(2.0-2.0*(*gas).visp) )*pow( (*gas).crref,2.0*(*gas).visp-1.0 );
          sigmat = (*DSMC).Vr;
	  if( sigmat - crm > 0.0)
               crm =  sigmat;

         if (sigmat/cells[(*DSMC).num].crm> dis_x(gen)) {
//	    (*DSMC).vr_pre_coll[0] = U1[(*DSMC).pair[0]]-U1[(*DSMC).pair[1]];
//	    (*DSMC).vr_pre_coll[1] = U2[(*DSMC).pair[0]]-U2[(*DSMC).pair[1]];
//	    (*DSMC).vr_pre_coll[2] = U3[(*DSMC).pair[0]]-U3[(*DSMC).pair[1]];
	    (*DSMC).alpha1 = dis_x(gen);//dis_x(gen);
	    (*DSMC).alpha2 = dis_x(gen);
	    (*DSMC).theta = acos( 2.0*(*DSMC).alpha1-1.0 );
	    (*DSMC).phi = 2.0*(*DSMC).pi*(*DSMC).alpha2;

	    (*DSMC).v_r_prime[0] = (*DSMC).Vr*cos((*DSMC).theta);
	    (*DSMC).v_r_prime[1] = (*DSMC).Vr*sin((*DSMC).theta)*cos((*DSMC).phi);
	    (*DSMC).v_r_prime[2] = (*DSMC).Vr*sin((*DSMC).theta)*sin((*DSMC).phi);

	    (*DSMC).v_c[0] = 0.5*( U1[(*DSMC).pair[0]]+ U1[(*DSMC).pair[1]] );
	    (*DSMC).v_c[1] = 0.5*( U2[(*DSMC).pair[0]]+ U2[(*DSMC).pair[1]] );
	    (*DSMC).v_c[2] = 0.5*( U3[(*DSMC).pair[0]]+ U3[(*DSMC).pair[1]] );

	    U1[(*DSMC).pair[0]] = (*DSMC).v_c[0] + 0.5*(*DSMC).v_r_prime[0];
	    U2[(*DSMC).pair[0]] = (*DSMC).v_c[1] + 0.5*(*DSMC).v_r_prime[1];
	    U3[(*DSMC).pair[0]] = (*DSMC).v_c[2] + 0.5*(*DSMC).v_r_prime[2];

	    U1[(*DSMC).pair[1]] = (*DSMC).v_c[0] - 0.5*(*DSMC).v_r_prime[0];
	    U2[(*DSMC).pair[1]] = (*DSMC).v_c[1] - 0.5*(*DSMC).v_r_prime[1];
	    U3[(*DSMC).pair[1]] = (*DSMC).v_c[2] - 0.5*(*DSMC).v_r_prime[2];
            num_act_col++;
	}
     }
   if(num_act_col < 1)
	printf("\n************            number of actual collisions = %d          *****************\n", num_act_col);
   cells[(*DSMC).num].crm = max( 1.1*crm, 0.9*cells[(*DSMC).num].crm );
  }
}

void Collision_DSMC_VHS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC){
 // omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());

  (*DSMC).pi = acos(-1.0);
  double sigmat;

  double crm;
  int num_act_col;
  int i,j;
  int N_cells;
  if((*box).problem == dcones)
	N_cells = (*box).Ny*( (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2] );
  else
	N_cells = (*box).N[0]*(*box).N[1]*(*box).N[2];
  for( (*DSMC).num = 0; (*DSMC).num<N_cells; (*DSMC).num++){
      crm = 0.0;
      std::uniform_real_distribution<> dis_x(0, 1);
      std::uniform_int_distribution<> disss(0, cells[(*DSMC).num].num_inside-1);

      (*DSMC).Nc = cells[(*DSMC).num].num_inside;
      if((*DSMC).Nc>1){
	 (*DSMC).Vc = ( cells[(*DSMC).num].dim[1]-cells[(*DSMC).num].dim[0] )*( cells[(*DSMC).num].dim[3]-cells[(*DSMC).num].dim[2] )*( cells[(*DSMC).num].dim[5]-cells[(*DSMC).num].dim[4] );
         (*DSMC).M_cand = (*DSMC).pi*(*gas).sigma*(*gas).sigma*( (double ) (*DSMC).Nc*cells[(*DSMC).num].n )*cells[(*DSMC).num].crm*(*gas).delta_t/2.0;
     	 num_act_col = 0;
	 //write_relaxation(cells,  box, num_act_col);
	for ( (*DSMC).cand_num = 0; (*DSMC).cand_num < floor((*DSMC).M_cand+dis_x(gen)); (*DSMC).cand_num++)
	{
		  (*DSMC).pair[0] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
		  (*DSMC).pair[1] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
		  while( (*DSMC).pair[0] == (*DSMC).pair[1] ){
			(*DSMC).pair[1] = cells[(*DSMC).num].indices_inside[ disss(gen) ];
		   }
		  (*DSMC).Vr = sqrt(    (U1[(*DSMC).pair[0]]-U1[(*DSMC).pair[1]])*(U1[(*DSMC).pair[0]]-U1[(*DSMC).pair[1]])
				      + (U2[(*DSMC).pair[0]]-U2[(*DSMC).pair[1]])*(U2[(*DSMC).pair[0]]-U2[(*DSMC).pair[1]])
				      + (U3[(*DSMC).pair[0]]-U3[(*DSMC).pair[1]])*(U3[(*DSMC).pair[0]]-U3[(*DSMC).pair[1]]) );
        	  sigmat = pow( (*DSMC).Vr,(2.0-2.0*(*gas).visp) )*pow( (*gas).crref,2.0*(*gas).visp-1.0 );
		 if( sigmat - crm > 0.0){
        	       crm =  sigmat;
		 }
         	if (sigmat/cells[(*DSMC).num].crm> dis_x(gen)) {
	    	(*DSMC).alpha1 = dis_x(gen);
	    	(*DSMC).alpha2 = dis_x(gen);
	    	(*DSMC).theta = acos( 2.0*(*DSMC).alpha1-1.0 );
	    	(*DSMC).phi = 2.0*(*DSMC).pi*(*DSMC).alpha2;

	    	(*DSMC).v_r_prime[0] = (*DSMC).Vr*cos((*DSMC).theta);
	    	(*DSMC).v_r_prime[1] = (*DSMC).Vr*sin((*DSMC).theta)*cos((*DSMC).phi);
	    	(*DSMC).v_r_prime[2] = (*DSMC).Vr*sin((*DSMC).theta)*sin((*DSMC).phi);

	    	(*DSMC).v_c[0] = 0.5*( U1[(*DSMC).pair[0]]+ U1[(*DSMC).pair[1]] );
	    	(*DSMC).v_c[1] = 0.5*( U2[(*DSMC).pair[0]]+ U2[(*DSMC).pair[1]] );
	    	(*DSMC).v_c[2] = 0.5*( U3[(*DSMC).pair[0]]+ U3[(*DSMC).pair[1]] );

	    	i = (*DSMC).pair[0]; j = (*DSMC).pair[1];

	    	U1[(*DSMC).pair[0]] = (*DSMC).v_c[0] + 0.5*(*DSMC).v_r_prime[0];
	    	U2[(*DSMC).pair[0]] = (*DSMC).v_c[1] + 0.5*(*DSMC).v_r_prime[1];
	    	U3[(*DSMC).pair[0]] = (*DSMC).v_c[2] + 0.5*(*DSMC).v_r_prime[2];

	    	U1[(*DSMC).pair[1]] = (*DSMC).v_c[0] - 0.5*(*DSMC).v_r_prime[0];
	    	U2[(*DSMC).pair[1]] = (*DSMC).v_c[1] - 0.5*(*DSMC).v_r_prime[1];
	    	U3[(*DSMC).pair[1]] = (*DSMC).v_c[2] - 0.5*(*DSMC).v_r_prime[2];
            	num_act_col++;
	  	}
	 }
         if((*DSMC).Nc < 5){
		printf("\n************            number of actual collisions = %d  at cell = %d step = %d  M_can=%e     *****************\n", num_act_col, (*DSMC).num, (*box).step, (*DSMC).M_cand);
		printf("Nc = %d, omega_max = %e, omega=%e",(*DSMC).Nc,cells[(*DSMC).num].crm,crm);
	}
   	cells[(*DSMC).num].crm = max( 1.1*crm, 0.9*cells[(*DSMC).num].crm );
      }
//      if((*DSMC).M_cand < 1.0)
//	printf("M_can=%e, Nc = %d, num = %d \n",(*DSMC).M_cand, (*DSMC).Nc,(*DSMC).num);
   }
}


void Collision_CBA_VHS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *CBA){
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());


  (*CBA).pi = acos(-1.0);
  double sigmat;
  double crm;
  int num_act_col;
  for( (*CBA).num = 0; (*CBA).num<(*box).N[0]*(*box).N[1]*(*box).N[2]; (*CBA).num++){
     crm = 0.0;
	std::uniform_real_distribution<> dis_x(0, 1);
	std::uniform_int_distribution<> disss(0, cells[(*CBA).num].num_inside-1);

	 //(*CBA).Vr_max = sqrt((*CBA).pi* (*gas).kb*cells[(*CBA).num].T/(*gas).m);
         //(*CBA).Vr_max = cr_max [ (*CBA).num ];
	 (*CBA).Nc = cells[(*CBA).num].num_inside;
	 (*CBA).Vc = ( cells[(*CBA).num].dim[1]-cells[(*CBA).num].dim[0] )*( cells[(*CBA).num].dim[3]-cells[(*CBA).num].dim[2] )*( cells[(*CBA).num].dim[5]-cells[(*CBA).num].dim[4] );
	 //(*CBA).M_cand = (*gas).Fn*(*CBA).Nc*(*CBA).Nc*(*CBA).pi*(*gas).sigma*(*gas).sigma*(*CBA).Vr_max*(*gas).delta_t/(2.0*(*CBA).Vc);
         //(*CBA).M_cand = ( (double ) (*CBA).Nc*cells[(*CBA).num].n )*(*CBA).pi*(*gas).sigma*(*gas).sigma*(*CBA).Vr_max*(*gas).delta_t/2.0;
         //(*CBA).M_cand = ( (double ) (*CBA).Nc*cells[(*CBA).num].n )*(*CBA).pi*(*gas).sigma*(*gas).sigma*thermal*pow( (*gas).crref/thermal,2.0*(*gas).visp-1.0 )*(*gas).delta_t/2.0;
         //(*CBA).M_cand = ( (double ) (*CBA).Nc*cells[(*CBA).num].n )*cells[ (*CBA).num ].crm*(*gas).delta_t/2.0;
         (*CBA).M_cand = (*CBA).pi*(*gas).sigma*(*gas).sigma*( (double ) (*CBA).Nc*cells[(*CBA).num].n )*cells[(*CBA).num].crm*(*gas).delta_t/2.0;

	 (*CBA).Y = increase_collision_rate(cells[(*CBA).num].n, gas);//(1-(*gas).b*rho/8.0)/ (pow((1.0 - (*gas).b*rho/4.0),3));
 	 (*CBA).M_cand = (*CBA).M_cand*(*CBA).Y;

     num_act_col = 0;
for ( (*CBA).cand_num = 0; (*CBA).cand_num < floor((*CBA).M_cand+dis_x(gen)); (*CBA).cand_num++)
	{
	  (*CBA).pair[0] = cells[(*CBA).num].indices_inside[ disss(gen) ];
	  (*CBA).pair[1] = cells[(*CBA).num].indices_inside[ disss(gen) ];
	  while( (*CBA).pair[0] == (*CBA).pair[1] ){
		(*CBA).pair[1] = cells[(*CBA).num].indices_inside[ disss(gen) ];
	   }
	  (*CBA).Vr = sqrt(    (U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]])*(U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]])
			      + (U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]])*(U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]])
			      + (U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]])*(U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]]) );
          sigmat = pow( (*CBA).Vr,(2.0-2.0*(*gas).visp) )*pow( (*gas).crref,2.0*(*gas).visp-1.0 );
	 if( sigmat - crm > 0.0){
               crm =  sigmat;
	 }
         if (sigmat/cells[(*CBA).num].crm> dis_x(gen)) {
	    (*CBA).vr_pre_coll[0] = U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]];
	    (*CBA).vr_pre_coll[1] = U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]];
	    (*CBA).vr_pre_coll[2] = U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]];

	    (*CBA).alpha1 = dis_x(gen);
	    (*CBA).alpha2 = dis_x(gen);
	    (*CBA).theta = acos( 2.0*(*CBA).alpha1-1.0 );
	    (*CBA).phi = 2.0*(*CBA).pi*(*CBA).alpha2;

//	    (*CBA).q = 1.0-2*(*CBA).alpha1;
//	    (*CBA).cos_teta = (*CBA).q;
//           (*CBA).sin_teta = sqrt(1.0-(*CBA).q*(*CBA).q);

	    (*CBA).v_r_prime[0] = (*CBA).Vr*cos((*CBA).theta);
	    (*CBA).v_r_prime[1] = (*CBA).Vr*sin((*CBA).theta)*cos((*CBA).phi);
	    (*CBA).v_r_prime[2] = (*CBA).Vr*sin((*CBA).theta)*sin((*CBA).phi);

	    (*CBA).v_c[0] = 0.5*( U1[(*CBA).pair[0]]+ U1[(*CBA).pair[1]] );
	    (*CBA).v_c[1] = 0.5*( U2[(*CBA).pair[0]]+ U2[(*CBA).pair[1]] );
	    (*CBA).v_c[2] = 0.5*( U3[(*CBA).pair[0]]+ U3[(*CBA).pair[1]] );

	    U1[(*CBA).pair[0]] = (*CBA).v_c[0] + 0.5*(*CBA).v_r_prime[0];
	    U2[(*CBA).pair[0]] = (*CBA).v_c[1] + 0.5*(*CBA).v_r_prime[1];
	    U3[(*CBA).pair[0]] = (*CBA).v_c[2] + 0.5*(*CBA).v_r_prime[2];

	    U1[(*CBA).pair[1]] = (*CBA).v_c[0] - 0.5*(*CBA).v_r_prime[0];
	    U2[(*CBA).pair[1]] = (*CBA).v_c[1] - 0.5*(*CBA).v_r_prime[1];
	    U3[(*CBA).pair[1]] = (*CBA).v_c[2] - 0.5*(*CBA).v_r_prime[2];


	    // Doing CBA correction
	    (*CBA).vr_post_coll[0] = U1[(*CBA).pair[0]] - U1[(*CBA).pair[1]];
	    (*CBA).vr_post_coll[1] = U2[(*CBA).pair[0]] - U2[(*CBA).pair[1]];
	    (*CBA).vr_post_coll[2] = U3[(*CBA).pair[0]] - U3[(*CBA).pair[1]];
//
	    (*CBA).norm_post_pre = sqrt( pow((*CBA).vr_post_coll[0]-(*CBA).vr_pre_coll[0],2) +
					 pow((*CBA).vr_post_coll[1]-(*CBA).vr_pre_coll[1],2) +
					 pow((*CBA).vr_post_coll[2]-(*CBA).vr_pre_coll[2],2));
	    (*CBA).d[0] = ((*CBA).vr_post_coll[0] - (*CBA).vr_pre_coll[0])*( (*gas).sigma/(*CBA).norm_post_pre);
	    (*CBA).d[1] = ((*CBA).vr_post_coll[1] - (*CBA).vr_pre_coll[1])*( (*gas).sigma/(*CBA).norm_post_pre);
	    (*CBA).d[2] = ((*CBA).vr_post_coll[2] - (*CBA).vr_pre_coll[2])*( (*gas).sigma/(*CBA).norm_post_pre);


	    x1[(*CBA).pair[0]] = x1[(*CBA).pair[0]] + (*CBA).d[0];
	    x2[(*CBA).pair[0]] = x2[(*CBA).pair[0]] + (*CBA).d[1];
	//    x3[(*CBA).pair[0]] = x3[(*CBA).pair[0]] + (*CBA).d[2]*factor;

	    x1[(*CBA).pair[1]] = x1[(*CBA).pair[1]] - (*CBA).d[0];
	    x2[(*CBA).pair[1]] = x2[(*CBA).pair[1]] - (*CBA).d[1];
	//    x3[(*CBA).pair[1]] = x3[(*CBA).pair[1]] - (*CBA).d[2]*factor;


            num_act_col++;
	  }
	}
   if(num_act_col < 1)
	printf("\n************            number of actual collisions = %d          *****************\n", num_act_col);
   cells[(*CBA).num].crm = max( 1.1*crm, 0.9*cells[(*CBA).num].crm );
   }
}

void Collision_CBA_HS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *CBA){
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams

//omp_set_num_threads( (*box).num_thr);
  std::random_device rd;
  std::mt19937 gen(rd());

  double factor;
  double rho = (*gas).n*(*gas).m;

  for( (*CBA).num = 0; (*CBA).num<(*box).N[0]*(*box).N[1]*(*box).N[2]; (*CBA).num++){

	std::uniform_real_distribution<> dis_x(0, 1);
	(*CBA).pi = acos(-1.0);

	  std::uniform_int_distribution<> disss(0, cells[(*CBA).num].num_inside-1);

         rho = cells[(*CBA).num].n*(*gas).m;
	 (*CBA).Vr_max = sqrt((*CBA).pi*(*gas).kb*cells[(*CBA).num].T/(*gas).m);
	 (*CBA).Nc = cells[(*CBA).num].num_inside;
	 (*CBA).Vc = ( cells[(*CBA).num].dim[1]-cells[(*CBA).num].dim[0] )*( cells[(*CBA).num].dim[3]-cells[(*CBA).num].dim[2] )*( cells[(*CBA).num].dim[5]-cells[(*CBA).num].dim[4] );
	 (*CBA).M_cand = (*gas).Fn*(*CBA).Nc*(*CBA).Nc*(*CBA).pi*(*gas).sigma*(*gas).sigma*(*CBA).Vr_max*(*gas).delta_t/(2.0*(*CBA).Vc);

	 (*CBA).Y = (1-(*gas).b*rho/8.0)/ (pow((1.0 - (*gas).b*rho/4.0),3));
 	 (*CBA).M_cand = (*CBA).M_cand*(*CBA).Y;

	// p_aim = (*gas).n*(*gas).kb*(*gas).T*(1.0 + (*CBA).Y*(*gas).b*(*gas).n);
	 factor = 1.0;
// 	 dummy = 0;
for ( (*CBA).cand_num = 0; (*CBA).cand_num < floor((*CBA).M_cand); (*CBA).cand_num++)
//	while( cand_num <  M_cand) //         check whether we have reached required number of collisions or not.
	{

	  (*CBA).pair[0] = cells[(*CBA).num].indices_inside[ disss(gen) ];
	  (*CBA).pair[1] = cells[(*CBA).num].indices_inside[ disss(gen) ];

	  (*CBA).r = dis_x(gen);
	  (*CBA).Vr = sqrt( (U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]])*(U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]]) + (U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]])*(U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]]) + (U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]])*(U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]]) );
	  if ((*CBA).Vr/(*CBA).Vr_max>(*CBA).r) {

	    if((*CBA).Vr - (*CBA).Vr_max > 1e-14){
		(*CBA).Vr_max	=  (*CBA).Vr;
		(*CBA).M_cand = (*gas).Fn*(*CBA).Nc*(*CBA).Nc*(*CBA).pi*(*gas).sigma*(*gas).sigma*(*CBA).Vr_max*(*gas).delta_t/(2.0*(*CBA).Vc);
		(*CBA).M_cand = (*CBA).M_cand*(*CBA).Y;
	    }

// 	    dummy++;
	    (*CBA).vr_pre_coll[0] = U1[(*CBA).pair[0]]-U1[(*CBA).pair[1]];
	    (*CBA).vr_pre_coll[1] = U2[(*CBA).pair[0]]-U2[(*CBA).pair[1]];
	    (*CBA).vr_pre_coll[2] = U3[(*CBA).pair[0]]-U3[(*CBA).pair[1]];
	    (*CBA).alpha1 = (*CBA).r;//dis_x(gen);
	    (*CBA).alpha2 = dis_x(gen);
	    (*CBA).phi = 2.0*(*CBA).pi*(*CBA).alpha2;
	    (*CBA).q = 1.0-2*(*CBA).alpha1;
	    (*CBA).cos_teta = (*CBA).q;
            (*CBA).sin_teta = sqrt(1.0-(*CBA).q*(*CBA).q);

	    (*CBA).v_r_prime[2] = (*CBA).Vr*(*CBA).sin_teta*cos((*CBA).phi);
	    (*CBA).v_r_prime[1] = (*CBA).Vr*(*CBA).sin_teta*sin((*CBA).phi);
	    (*CBA).v_r_prime[0] = (*CBA).Vr*(*CBA).cos_teta;

	    (*CBA).v_c[0] = 0.5*( U1[(*CBA).pair[0]]+ U1[(*CBA).pair[1]] );
	    (*CBA).v_c[1] = 0.5*( U2[(*CBA).pair[0]]+ U2[(*CBA).pair[1]] );
	    (*CBA).v_c[2] = 0.5*( U3[(*CBA).pair[0]]+ U3[(*CBA).pair[1]] );

	    U1[(*CBA).pair[0]] = (*CBA).v_c[0] + 0.5*(*CBA).v_r_prime[0];
	    U2[(*CBA).pair[0]] = (*CBA).v_c[1] + 0.5*(*CBA).v_r_prime[1];
	    U3[(*CBA).pair[0]] = (*CBA).v_c[2] + 0.5*(*CBA).v_r_prime[2];

	    U1[(*CBA).pair[1]] = (*CBA).v_c[0] - 0.5*(*CBA).v_r_prime[0];
	    U2[(*CBA).pair[1]] = (*CBA).v_c[1] - 0.5*(*CBA).v_r_prime[1];
	    U3[(*CBA).pair[1]] = (*CBA).v_c[2] - 0.5*(*CBA).v_r_prime[2];

	    // Doing CBA correction
	    (*CBA).vr_post_coll[0] = U1[(*CBA).pair[0]] - U1[(*CBA).pair[1]];
	    (*CBA).vr_post_coll[1] = U2[(*CBA).pair[0]] - U2[(*CBA).pair[1]];
	    (*CBA).vr_post_coll[2] = U3[(*CBA).pair[0]] - U3[(*CBA).pair[1]];
//
	    (*CBA).norm_post_pre = sqrt( pow((*CBA).vr_post_coll[0]-(*CBA).vr_pre_coll[0],2) + pow((*CBA).vr_post_coll[1]-(*CBA).vr_pre_coll[1],2) + pow((*CBA).vr_post_coll[2]-(*CBA).vr_pre_coll[2],2));
	    (*CBA).d[0] = ((*CBA).vr_post_coll[0] - (*CBA).vr_pre_coll[0])*( (*gas).sigma/(*CBA).norm_post_pre);
	    (*CBA).d[1] = ((*CBA).vr_post_coll[1] - (*CBA).vr_pre_coll[1])*( (*gas).sigma/(*CBA).norm_post_pre);
	    (*CBA).d[2] = ((*CBA).vr_post_coll[2] - (*CBA).vr_pre_coll[2])*( (*gas).sigma/(*CBA).norm_post_pre);


	    x1[(*CBA).pair[0]] = x1[(*CBA).pair[0]] + (*CBA).d[0];
	    x2[(*CBA).pair[0]] = x2[(*CBA).pair[0]] + (*CBA).d[1];
	//    x3[(*CBA).pair[0]] = x3[(*CBA).pair[0]] + (*CBA).d[2]*factor;

	    x1[(*CBA).pair[1]] = x1[(*CBA).pair[1]] - (*CBA).d[0];
	    x2[(*CBA).pair[1]] = x2[(*CBA).pair[1]] - (*CBA).d[1];
	//    x3[(*CBA).pair[1]] = x3[(*CBA).pair[1]] - (*CBA).d[2]*factor;

	  }
	}
       }
}
