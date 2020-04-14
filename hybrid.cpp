#include "cpp_headers.h"

double compute_DKL_MED_MW(double *lambda){
  double a = 10.0;
  double pi = acos(-1.0);

  return (exp(-a*a/2.0)*(sqrt(2.0)*sqrt(pi)*exp(a*a/2.0)*erf(a/sqrt(2.0))-2.0*a)*(2*lambda[1]-1))/(pow(2.0,3/2)*sqrt(pi));
}
double compute_fisher_metric(double *lambda){
  // exp(-x^2/2)/(sqrt(2*pi))*(-x+l1+2*l2*x+3*l3*x^2)^2
  double a = 10.0;
  double pi = acos(-1.0);
  double metric = 0.0;
  double l_1 = lambda[0];
  double l_2 = lambda[1];
  double l_3 = lambda[2];

  metric = -( exp(-a*a/2.0)*( (9.0*sqrt(2.0)*pow(a,3)+27.0*sqrt(2.0)*a)*pow(l_3,2.0)+3.0*pow(2.0,3.0/2.0)*a*l_1*l_3
                +pow(2.0,5.0/2.0)*a*pow(l_2,2.0)-pow(2.0,5.0/2.0)*a*l_2+sqrt(2.0)*a) )/sqrt(pi);
  metric += +27.0*erf(a/sqrt(2.0))*pow(l_3,2.0)+6.0*erf(a/sqrt(2.0))*l_1*l_3+4.0*erf(a/sqrt(2.0))*pow(l_2,2.0)-4.0*erf(a/sqrt(2.0))*l_2+erf(a/sqrt(2.0))*pow(l_1,2)+erf(a/sqrt(2.0));
  return metric;
}


void evaluation_GPR(struct GPR *gpr, double x, double *y){

  double *output = (double *) malloc( dimY * sizeof(double) );
  int i,j;
  int Ni,Nj;
  double input, r, dummy;
  /*
  def transf_to(xx,m,v,i1,i2):
    xxnew = xx.copy();
    for i in range(i1,i2):
        xxnew[:,i] = (xx[:,i] - m[i] ) / v[i] ** 0.5+1.0
    return xxnew

  def transf_back(xx,m,v,i1,i2):
    xxnew = xx.copy();
    for i in range(i1,i2):
        xxnew[:, i] = (xx[:,i]-1.0)* v[i] ** 0.5 + m[i]
    return xxnew
    */
  input = (x - (*gpr).mean_input[0]) / pow( (*gpr).variance_input[0], 0.5) + 1.0 ;
  Ni = Ndata;
  for( i=0; i<Ni; i++){
    r = fabs((*gpr).X[i] - input) / (*gpr).k_lsc;
    (*gpr).kxX[i] =  (*gpr).k_var * exp(-r*r/2.0);
  }
  // np.matmul(kxX, kst)
  Ni = Ndata; Nj = dimY;
  for( j=0; j<Nj; j++){
    dummy = 0.0;
    for(i=0; i<Ni; i++){
      dummy += (*gpr).kxX[i]*(*gpr).kst[i*Nj+j];
    }
    output[j] = dummy;
  }

  for(j=0; j<dimY; j++){
    y[j] = (output[j]-1.0)*pow( (*gpr).variance_output[j], 0.5) + (*gpr).mean_output[j];
  }
  //printf("y computed\n");
  free(output);
}

void evaluation(struct NerualNetwork *NN, double x, double *y){
  double *h1 = (double *) malloc( dimX*L11 * sizeof(double) );
  double *h2 = (double *) malloc( dimX*L22 * sizeof(double) );
  double *h3 = (double *) malloc( dimX*L33 * sizeof(double) );
  double *output = (double *) malloc( dimY * sizeof(double) );
  int i,j;
  int Ni,Nj;
  double input;

  input = (x - (*NN).mean_input[0]) / (*NN).dev_input[0];

  Ni = dimX; Nj = L11;
  for(i=0; i<Ni; i++){
    for(j=0; j<Nj; j++){
      h1[i*Nj+j] = tanh( input*(*NN).W1[i*Nj+j] + (*NN).b1[i*Nj+j]);
    }
  }
  //printf("h1 computed\n");

  double sum;
  for(j=0; j<L22; j++){
    sum = 0.0;
    for(i=0; i<L11; i++){
        sum+= h1[i]*(*NN).W2[i*L22+j];
    }
    h2[j] = tanh( sum + (*NN).b2[j] );
  }
  //printf("h2 computed\n");

  for(j=0; j<L33; j++){
    sum = 0.0;
    for(i=0; i<L22; i++){
        sum+= h2[i]*(*NN).W3[i*L33+j];
    }
    h3[j] = tanh( sum + (*NN).b3[j] );
  }
  //printf("h3 computed\n");

  for(j=0; j<dimY; j++){
    sum = 0.0;
    for(i=0; i<L33; i++){
        sum+= h3[i]*(*NN).Wo[i*dimY+j];
    }
    output[j] = sum + (*NN).bo[j];
  }


  for(j=0; j<dimY; j++){
    y[j] = output[j]*(*NN).dev_output[j] + (*NN).mean_output[j];
  }
  //printf("y computed\n");

  free(h1); free(h2); free(h3);
  free(output);
}


double MED(double x, double *lambda, int N, double Z){
  int i;
  double sum= 0.0;
  for(i=0; i<N; i++){
    sum += -pow(x,i+1)*lambda[i];
  }
  sum = exp(sum);
  return sum;

}

double MEDx(double x, double *lambda, int N, double Z, double SC){
  int i;
  double sum= 0.0;
  for(i=0; i<N; i++){
    sum += -pow(x,i+1)*lambda[i];
  }
  sum = (x-SC)*exp(sum);
  return sum;

}


double ChEn(double x, double rho, double RT, double q, double tau){
  double phi;
  phi = 2.0/5.0*q*x/(rho*pow(RT,2))*( x*x/2.0/RT - 5.0/2.0 );
  phi = phi - tau*x*x/rho/pow(RT,2.0);
  //printf("phi = %e    ", phi);
  return  exp( -(x*x)/2.0 )/sqrt(2.0*acos(-1.0))*(1.0+phi);
}

double Normal(double x){
  return exp( -(x*x)/2.0 )/sqrt(2.0*acos(-1.0));
}
double Normalm(double x,double m){
  return exp( -(x-m)*(x-m)/2.0 )/sqrt(2.0*acos(-1.0));
}

double compute_Mom(double *lambda, int N, int pth){
  // for Z, pth should be 0

 int nn = 50;
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
  int nn = 20;
  double xi[20] = {-0.9931286,  -0.96397193, -0.91223443, -0.83911697, -0.74633191, -0.63605368,
 -0.510867,   -0.37370609, -0.22778585, -0.07652652,  0.07652652,  0.22778585,
  0.37370609,  0.510867,    0.63605368,  0.74633191,  0.83911697,  0.91223443,
  0.96397193,  0.9931286};
  double wi[20] = {0.01761401, 0.04060143, 0.06267205, 0.08327674, 0.10193012, 0.11819453,
 0.13168864, 0.14209611, 0.14917299, 0.15275339, 0.15275339, 0.14917299,
 0.14209611, 0.13168864, 0.11819453, 0.10193012, 0.08327674, 0.06267205,
 0.04060143, 0.01761401};
*/
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

 double a, b;

b = 10.0;
a = -10.0;
 double sum = 0.0, w, x;
 int i;
 for(i=0; i<nn; i++){
   x = (b-a)/2.0*xi[i]+(a+b)/2.0;
   w = wi[i]*(b-a)/2.0;
   sum += w*MED(x, lambda, N, 1.0)*pow(x,pth);
 }
 //free(xi);
 //free(wi);
 return sum;
}

double compute_mom_flux(double *lambda, int N, double mu, double sigma, double Z){
  // SC is mu/sigma
  // for Z, pth should be 0

  int nn = 50;
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
  int nn = 20;
  double xi[20] = {-0.9931286,  -0.96397193, -0.91223443, -0.83911697, -0.74633191, -0.63605368,
 -0.510867,   -0.37370609, -0.22778585, -0.07652652,  0.07652652,  0.22778585,
  0.37370609,  0.510867,    0.63605368,  0.74633191,  0.83911697,  0.91223443,
  0.96397193,  0.9931286};
  double wi[20] = {0.01761401, 0.04060143, 0.06267205, 0.08327674, 0.10193012, 0.11819453,
 0.13168864, 0.14209611, 0.14917299, 0.15275339, 0.15275339, 0.14917299,
 0.14209611, 0.13168864, 0.11819453, 0.10193012, 0.08327674, 0.06267205,
 0.04060143, 0.01761401};
 */
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

 double SC = mu/sigma;
 double b = 10.0;
 double a = -SC;
 int pth = 1;

 double sum = 0.0, w, x;
 int i;
 for(i=0; i<nn; i++){
   x = (b-a)/2.0*xi[i]+(a+b)/2.0;
   w = wi[i]*(b-a)/2.0;
   sum += w*(MED(x, lambda, N, 1.0)*pow(x+SC,pth));
 }
 sum = sum*sigma/Z;
// printf("SC=%e, Z=%e, sigma=%e, a=%e\n", SC, Z, sigma, a);
 //free(xi);
 //free(wi);
 return sum;
}

void MomentsofSamples(int nSamples, double *x, double *m, int Nm){
  int i,j;
  for(i=0; i<Nm; i++){
    m[i] = 0.0;
  }
  for(j=0; j<nSamples; j++){
    for(i=0; i<Nm; i++){
        m[i] += pow(x[j],i+1);
    }
  }
  for(i=0; i<Nm; i++){
    m[i] = m[i]/(1.0*nSamples);
  }
}

void MetroHastingChEn(int nSamples, double rho, double RT, double q, double tau, double *x){
  //*nSamples = max(10000,*nSamples);
  //double *x = (double *) malloc( (nSamples) * sizeof(double) );
  double xstar, c, alpha, u;
  int t, i;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::normal_distribution<> dis_norm(0.0,1.0);
  for(i=0; i<nSamples; i++){
    x[i] = 0.0;
  }
  int burn = 1, ib = 0;
  double xx, x_prev;
  t = 0;
  x_prev = 0.0;
  xx = x_prev;
  while(t<nSamples){
    xstar = dis_norm(gen);
    c = Normal(x_prev)/Normal(xstar);
    alpha = min(1.0, ChEn(xstar, rho, RT, q, tau)/ChEn(x_prev, rho, RT, q, tau)*c);
    //alpha = min(1.0, MED(xstar, lambda, N, Z)/MED(x_prev, lambda, N, Z)*c);
    // accept?
    u = dis_r(gen);
    if(u < alpha){
        //x[t] = xstar;
        x_prev = xx;
        xx = xstar;
    }
    else{
        //x[t] = x[t-1];
        xx = x_prev;
    }
    if(ib>burn){
      x[t] = xx;
      t++;
    }
    ib++;
  }
  //return x;
}


void MetroHastingFace(int nSamples, double *lambda, int N, double *x, double SC){
  //*nSamples = max(10000,*nSamples);
  //double *x = (double *) malloc( (nSamples) * sizeof(double) );
  double xstar, c, alpha, u;
  int t, i;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::normal_distribution<> dis_norm(0.0,1.0);
  double Z = 1.0;//compute_Mom(lambda, N, 0);
  for(i=0; i<nSamples; i++){
    x[i] = 0.0;
  }
  int burn = 1000, ib = 0;
  double xx, x_prev;
  t = 0;
  x_prev = 0.0;
  xx = x_prev;
  while(t<nSamples){
    xstar = dis_norm(gen);
    c = Normal(x_prev)/Normal(xstar);
    alpha = min(1, MEDx(xstar, lambda, N, Z, SC)/MEDx(x_prev, lambda, N, Z, SC)*c);
    // accept?
    u = dis_r(gen);
    if(u < alpha ){
        //x[t] = xstar;
        x_prev = xx;
        xx = xstar;
    }
    else{
        //x[t] = x[t-1];
        xx = x_prev;
    }
    if(ib>burn){
      x[t] = xx;
      t++;
    }
    ib++;
  }
  //return x;
}

void MetroHasting(int nSamples, double *lambda, int N, double *x){
  //*nSamples = max(10000,*nSamples);
  //double *x = (double *) malloc( (nSamples) * sizeof(double) );
  double xstar, c, alpha, u;
  int t, i;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_r(0.0, 1.0);
  std::normal_distribution<> dis_norm(0.0,1.0);
  double Z = 1.0;//compute_Mom(lambda, N, 0);
  for(i=0; i<nSamples; i++){
    x[i] = 0.0;
  }
  int burn = 1000, ib = 0;
  double xx, x_prev;
  t = 0;
  x_prev = 0.0;
  xx = x_prev;
  while(t<nSamples){
    xstar = dis_norm(gen);
    c = Normal(x_prev)/Normal(xstar);
    alpha = min(1, MED(xstar, lambda, N, Z)/MED(x_prev, lambda, N, Z)*c);
    // accept?
    u = dis_r(gen);
    if(u < alpha){
        //x[t] = xstar;
        x_prev = xx;
        xx = xstar;
    }
    else{
        //x[t] = x[t-1];
        xx = x_prev;
    }
    if(ib>burn){
      x[t] = xx;
      t++;
    }
    ib++;
  }
  //return x;
}

void SPH_velocity_temperature_update_hybrid(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){
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
  mu0 = (*gas).mu;
  gama = 1.4;
  r_cut = (*gas).r_cut;
  h = (*gas).h;
  mb = (*gas).mb;
  search = floor(r_cut/(*box).delta_dim[1])+1;
  sqrt_pi = sqrt(2.0*acos(-1.0));

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
          A[i*N+j] = 0.0;
          ddx[i*N+j] = 0.0;
          d2dx2[i*N+j] = 0.0;
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

          if(i!=j){
            DW = W*r1r2*(-1.0/h/h);
          }
          else{
            DW = 0.0;
          }
          dmudx[i] += mb/rho[j]*DW*mu0*pow(T[j]/(*gas).T0,0.5);

          ddx[i*N+j] = mb/rho[j]*DW;
          d2dx2[i*N+j] = mb/rho[j]*D2W;

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
      A[i+j*N] += -(*gas).delta_t*(4.0/3.0)/rho[i]*dmudx[i]*ddx[i*N+j];
      A[i+j*N] += -(*gas).delta_t*(4.0/3.0)/rho[i]*mu*d2dx2[i*N+j];
    }
    rA[i] = U2[i];
  }

  int nrhs = 1;
  int lda=N;
  int ipiv[N];
  int ldb = N;
  int info;

  //dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);
  for(i=0; i<N; i++)
    vs[i] = rA[i];

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
    rA [i] = - rho[i]/(*gas).delta_t;
    rA [i] += rho[i]*dusdx[i];
    //rA [i] += - rho[i]/T[i]*dTdt;

    for(j=0; j<N; j++){
      A[i+j*N] = 0.0;
      if(i==j){
          A[i+j*N] = -1.0/( (*gas).delta_t*R*T[i] );
      }
      A[i+j*N] += -(*gas).delta_t/rho[i]*drhodx[i]*ddx[i*N+j];
      A[i+j*N] += (*gas).delta_t*d2dx2[i*N+j];
    }
  }

  //Gauss_Jordan(B, rB, N, pnp1);

  double nrA = 0.0;
  for(i=0; i<N; i++)
      for(j=0; j<N; j++)
          A0[i*N+j] = A[i+j*N];
  for(i=0; i<N; i++){
    rA0[i] = rA[i];
    nrA += rA[i];
  }

  //for(i=0; i<N; i++)
  //  for(j=0; j<N; j++)
  //    a[i+j*N] = A[i*N+j];
  //dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);

  for(i=0; i<N; i++)
    pnp1[i] = rA[i];

  res = 0.0;
  for(i=0; i<N; i++){
      res += rA0[i];
      for(j=0; j<N; j++){
        res += -A0[i*N+j]*pnp1[j];
      }
  }



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
    U2[i]   = vs[i] - (*gas).delta_t/rho[i]*dpdx[i];
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
      A[i+j*N] += -(*gas).delta_t/rho[i]/R*(gama-1.0)*dkdx*ddx[i*N+j];
      A[i+j*N] += (*gas).delta_t/rho[i]/R*(gama-1.0)*k*d2dx2[i*N+j];
    }
    rA [i] = T[i] + (*gas).delta_t*(gama-1.0)/R/rho[i]*( -p*dudx[i] + 4.0/3.0*mu/rho[i]*dudx[i]*dudx[i] );
  }



  nrA = 0.0;
  for(i=0; i<N; i++)
      for(j=0; j<N; j++)
          A0[i*N+j] = A[i+j*N];
  for(i=0; i<N; i++){
    rA0[i] = rA[i];
    nrA += rA[i];
  }

  //dgesv(&N, &nrhs, A, &lda, ipiv, rA, &ldb, &info);

  for(i=0; i<N; i++)
    T[i] = rA[i];

  res = 0.0;
  for(i=0; i<N; i++){
      res += rA0[i];
      for(j=0; j<N; j++){
        res += -A0[i*N+j]*T[j];
      }
  }

  // update rho
  for(i=0; i<N; i++)
      rho[i] = pnp1[i]/R/T[i];
  free(rA); free(A);// free(a);
  free(vs); free(pnp1);
  free(rA0); free(dmudx);
  free(dusdx); free(dudx); free(drhodx);
  free(d2Tdx2); free(dTdx); free(dpdx);
  free(A0); free(ddx); free(d2dx2);

}

void hybrid_surface_flux(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho,  struct GPR *gpr, struct NerualNetwork *NN){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j,Nnew1, Nnew3;
int Nnew1_SPH, Nnew3_SPH;
double tol = 1e-15;
double Cratio, Hratio;
int NnewC,NnewH,N_new, Ntot;
int num;
double dt_remain,yw;
int Nout = 0;
double pi = acos(-1.0);
double d1, d3;
int ind_h1 = -1, ind_h2 = -1;
for(num=0; num<(*box).N[1]; num++){
  if(cells[num].model == hybrid && ind_h1 == -1){
    ind_h1 = num;
  }
  else if(cells[num].model == hybrid && ind_h1 != -1){
    ind_h2 = num;
  }
}
if(ind_h1 != -1 && ind_h2 == -1){
  ind_h2 = ind_h1;
}
if(ind_h1 != -1){
//  ind_h1 = ind_h1 - 1;
//  ind_h2 = ind_h2 + 1;
double std1 = sqrt((*gas).kb*cells[ind_h1-1].T/(*gas).m);
double std3 = sqrt((*gas).kb*cells[ind_h2+1].T/(*gas).m);

// Now, how many new particles enter domain
double beta, dummy, Z;
double xx, m, stdT, VMP;
double sign;
double *U_temp;
int N;
double lambda[3], mm[3];
/*
beta = 1.0/(std1*sqrt(2.0));
//beta = 1.0/(stdc);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*cells[ind_h1].n/(beta*2.0*sqrt(pi));
Nnew1 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;
*/
xx = 0.5*cells[ind_h1-1].Q[1]/( cells[ind_h1-1].n*(*gas).m*pow( (*gas).kb*cells[ind_h1-1].T/(*gas).m, 1.5 ) );
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
Z = compute_Mom(lambda, 3, 0);
dummy = compute_mom_flux(lambda, 3, cells[ind_h1-1].U_space[1], std1, Z);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*cells[ind_h1-1].n*dummy;
Nnew1 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;

/*
beta = 1.0/(std3*sqrt(2.0));
//beta = 1.0/(stdh);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*cells[ind_h2].n/(beta*2.0*sqrt(pi));
Nnew3 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;
*/

xx = 0.5*cells[ind_h2+1].Q[1]/( cells[ind_h2+1].n*(*gas).m*pow( (*gas).kb*cells[ind_h2+1].T/(*gas).m, 1.5 ) );
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
Z = compute_Mom(lambda, 3, 0);
dummy = compute_mom_flux(lambda, 3, -cells[ind_h2+1].U_space[1], std3, Z);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*cells[ind_h2+1].n*dummy;
Nnew3 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;

Nnew1_SPH = 0;
Nnew3_SPH = 0;
/*
d1 = cells[ind_h1].U_space[1]*(*gas).delta_t;
if(d1<0.0){
  Nnew1_SPH = floor(cells[ind_h1].n*fabs(d1)/( (*gas).mb/(*gas).m ));
  //printf("d1=%e, Nnew1_SPH=%d\n", d1, Nnew1_SPH);
}
d3 = cells[ind_h2].U_space[1]*(*gas).delta_t;
if(d3>0.0){
  Nnew3_SPH = floor(cells[ind_h2].n*d3/( (*gas).mb/(*gas).m ));
  //printf("d3=%e, Nnew3_SPH=%d\n", d3, Nnew3_SPH);
}
*/

N_new = Nnew1 + Nnew3 + Nnew1_SPH + Nnew3_SPH;
if(N_new>0){
  //printf("**   N_new1=%d and Nnew3=%d\n",Nnew1, Nnew3);
  Ntot = N_new + (*gas).N;


double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew, *Tnew, *rhonew;
int *colornew;
//if(Ntot != (*gas).N){
  x1new = (double *) malloc( Ntot * sizeof(double) );
  x2new = (double *) malloc( Ntot * sizeof(double) );
  U1new = (double *) malloc( Ntot * sizeof(double) );
  U2new = (double *) malloc( Ntot * sizeof(double) );
  U3new = (double *) malloc( Ntot * sizeof(double) );
  x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
  x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
  colornew = (int *) malloc( Ntot * sizeof(int) );

  Tnew = (double *) malloc( Ntot * sizeof(double) );
  rhonew = (double *) malloc( Ntot * sizeof(double) );
//}
// store the surviving partices
j=0;
for(i=0; i<(*gas).N; i++){
   x1new[j] = (*x1)[i];
   x2new[j] = (*x2)[i];
   U1new[j] = (*U1)[i];
   U2new[j] = (*U2)[i];
   U3new[j] = (*U3)[i];
   colornew[j] = (*color)[i];
   x1_oldnew[j] = (*x1_old)[i];
   x2_oldnew[j] = (*x2_old)[i];
   Tnew[j] = (*T)[i];
   rhonew[j] = (*rho)[i];
   j++;
}
N = Nnew1;
if(Nnew3>N){
  N = Nnew3;
}
U_temp = (double *) malloc( N * sizeof(double) );

stdT = sqrt( (*gas).kb*cells[ind_h1-1].T/(*gas).m );
VMP = stdT;//sqrt(2.0*(*gas).kb*cells[ind_h1-1].T/(*gas).m);
m = cells[ind_h1-1].U_space[1];
xx = 0.5*cells[ind_h1-1].Q[1]/( cells[ind_h1-1].n*(*gas).m*pow( (*gas).kb*cells[ind_h1-1].T/(*gas).m, 1.5 ) );
//evaluation(NN, xx, lambda);
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
sign = 1.0;
MetroHastingFace(Nnew1, lambda, 3, U_temp, -m/stdT);
//MomentsofSamples(Nnew1, U_temp, mm, 3);
//printf("mm before %lf %lf %lf\n", mm[0], mm[1], mm[2]);
for(i=0; i<Nnew1; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U2new[j] = U_temp[i]*VMP;//+m;
    //if(U2new[j] < 0.0){
    //  printf("Oops! U2new = %e < 0.0 for ind_h1\n", U2new[j]);
    //}
    x2new[j] = cells[ind_h1].dim[2] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h1-1].T;
    rhonew[j] = cells[ind_h1-1].n*(*gas).m;
		j++;
}


stdT = sqrt( (*gas).kb*cells[ind_h2+1].T/(*gas).m );
VMP = stdT;//sqrt(2.0*(*gas).kb*cells[ind_h2+1].T/(*gas).m);
m = cells[ind_h2+1].U_space[1];
xx = 0.5*cells[ind_h2+1].Q[1]/( cells[ind_h2+1].n*(*gas).m*pow( (*gas).kb*cells[ind_h2+1].T/(*gas).m, 1.5 ) );
//evaluation(NN, xx, lambda);
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
sign = -1.0;
MetroHastingFace(Nnew3, lambda, 3, U_temp, m/stdT);
for(i=0; i<Nnew3; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U2new[j] = -U_temp[i]*VMP;//-m;
    //if(U2new[j] > 0.0){
    //  printf("Oops! U2new = %e > 0.0 for ind_h2\n", U2new[j]);
    //}
    x2new[j] = cells[ind_h2].dim[3] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h2+1].T;
    rhonew[j] = cells[ind_h2+1].n*(*gas).m;
		j++;
}
free(U_temp);

/*
double stdT;
double VMP, SC, QA, FS1, FS2, U, pratio,FS1m;
VMP = sqrt(2.0*(*gas).kb*cells[ind_h1-1].T/(*gas).m);
SC = cells[ind_h1-1].U_space[1]/VMP;
QA = 3.0;
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));

stdT = sqrt( (*gas).kb*cells[ind_h1-1].T/(*gas).m );
for(i=0; i<Nnew1; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		}
		U2new[j] = (U+SC)*VMP;
    x2new[j] = cells[ind_h1].dim[2] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h1-1].T;
    rhonew[j] = cells[ind_h1-1].n*(*gas).m;
		j++;
}

VMP = sqrt(2.0*(*gas).kb*cells[ind_h2+1].T/(*gas).m);
SC = cells[ind_h2+1].U_space[1]/VMP;
//VMP = sqrt(2.0*(*gas).kb*cells[(*box).N[1]-1].T/(*gas).m);
//stdT = sqrt( (*gas).kb*cells[(*box).N[1]-1].T/(*gas).m );
stdT = sqrt( (*gas).kb*cells[ind_h2+1].T/(*gas).m );
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));
for(i=0; i<Nnew3; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1m-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1m-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		}
		U2new[j] = (U+SC)*VMP;
		x2new[j] = cells[ind_h2].dim[3] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h2+1].T;
    rhonew[j] =cells[ind_h2+1].n*(*gas).m;
		j++;
}
*/

for(i=0; i<Nnew1_SPH; i++){
  x2new[j] = cells[ind_h1].dim[2] + (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew1_SPH);
  U1new[j] = cells[ind_h1].U_space[0];
  U2new[j] = cells[ind_h1].U_space[1];
  U3new[j] = cells[ind_h1].U_space[2];
  x1new[j] = (*x1)[0];
  x1_oldnew[j] = x1new[j];
  x2_oldnew[j] = x2new[j];
  colornew[j] = sph;
  Tnew[j] = cells[ind_h1].T;
  rhonew[j] = cells[ind_h1].n*(*gas).m;
  j++;
}

for(i=0; i<Nnew3_SPH; i++){
  x2new[j] = cells[ind_h2].dim[3] + (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew3_SPH);
  U1new[j] = cells[ind_h2].U_space[0];
  U2new[j] = cells[ind_h2].U_space[1];
  U3new[j] = cells[ind_h2].U_space[2];
  x1new[j] = (*x1)[0];
  x1_oldnew[j] = x1new[j];
  x2_oldnew[j] = x2new[j];
  colornew[j] = sph;
  Tnew[j] = cells[ind_h2].T;
  rhonew[j] = cells[ind_h2].n*(*gas).m;
  j++;
}




if(j!=Ntot){
  printf("---- in addrem_BC_hybrid j=%d Ntot=%d\n",j,Ntot);
}
/////  reallocating memory
free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
free(*color); free(*x1); free(*x1_old); free(*T); free(*rho);
(*gas).N = Ntot;
*x1 = x1new;
*x2 = x2new;
*T = Tnew;
*rho = rhonew;
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
}
}
}



void add_ghost_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho,  struct GPR *gpr, struct NerualNetwork *NN){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j;
int Nnew1_SPH, Nnew3_SPH;
int Nnew1_DSMC, Nnew3_DSMC;
double tol = 1e-15;
int N_new, Ntot;
int num;
double dt_remain,yw;
int Nout = 0;
double pi = acos(-1.0);
double d1, d3;
int ind_h1 = -1, ind_h2 = -1;
for(num=0; num<(*box).N[1]; num++){
  if(cells[num].model == hybrid && ind_h1 == -1){
    ind_h1 = num;
  }
  else if(cells[num].model == hybrid && ind_h1 != -1 && ind_h2 == -1){
    ind_h2 = num;
  }
}
//if(ind_h1 != -1 && ind_h2 == -1){
//  ind_h2 = ind_h1;
//}
if(ind_h1 != -1 ){
double std1 = sqrt((*gas).kb*cells[ind_h1-1].T/(*gas).m);
double std3;// = sqrt((*gas).kb*cells[ind_h2+1].T/(*gas).m);

// Now, how many new particles enter domain
Nnew1_DSMC = floor( cells[ind_h1-1].n_smooth*cells[ind_h1-1].volume/(*gas).Fn+ dis_x(gen));
//printf("ind_h1=%d & ind_h2 = %d\n", ind_h1, ind_h2);
if( ind_h2!=-1 && ind_h2+1 < (*box).N[1]){
  Nnew3_DSMC = floor( cells[ind_h2+1].n_smooth*cells[ind_h2+1].volume/(*gas).Fn+ dis_x(gen));
  std3 = sqrt((*gas).kb*cells[ind_h2+1].T/(*gas).m);
}
else{
  Nnew3_DSMC = 0;
}
//Nnew3_DSMC = 0;
Nnew1_SPH = 0;//floor(cells[ind_h1].n*cells[ind_h1].volume/( (*gas).mb/(*gas).m ));
Nnew3_SPH = 0;//floor(cells[ind_h2].n*cells[ind_h2].volume/( (*gas).mb/(*gas).m ));

N_new = Nnew1_DSMC + Nnew3_DSMC + Nnew1_SPH + Nnew3_SPH;
if(N_new>0){
  //printf("**   N_new=%d and Nout=%d\n", N_new, Nout);
  Ntot = N_new + (*gas).N;


double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew, *Tnew, *rhonew;
int *colornew;
//if(Ntot != (*gas).N){
  x1new = (double *) malloc( Ntot * sizeof(double) );
  x2new = (double *) malloc( Ntot * sizeof(double) );
  U1new = (double *) malloc( Ntot * sizeof(double) );
  U2new = (double *) malloc( Ntot * sizeof(double) );
  U3new = (double *) malloc( Ntot * sizeof(double) );
  x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
  x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
  colornew = (int *) malloc( Ntot * sizeof(int) );

  Tnew = (double *) malloc( Ntot * sizeof(double) );
  rhonew = (double *) malloc( Ntot * sizeof(double) );
//}
// store the surviving partices
j=0;
for(i=0; i<(*gas).N; i++){
   x1new[j] = (*x1)[i];
   x2new[j] = (*x2)[i];
   U1new[j] = (*U1)[i];
   U2new[j] = (*U2)[i];
   U3new[j] = (*U3)[i];
   colornew[j] = (*color)[i];
   x1_oldnew[j] = (*x1_old)[i];
   x2_oldnew[j] = (*x2_old)[i];
   Tnew[j] = (*T)[i];
   rhonew[j] = (*rho)[i];
   j++;
}
double xx, m;
double sign;
double *U_temp;
int N = Nnew1_DSMC;
double lambda[3], mm[3];
if(Nnew3_DSMC>N){
  N = Nnew3_DSMC;
}
U_temp = (double *) malloc( N * sizeof(double) );

//VMP = stdT;//sqrt(2.0*(*gas).kb*cells[ind_h1-1].T/(*gas).m);
m = cells[ind_h1-1].U_space[1];
xx = 0.5*cells[ind_h1-1].Q[1]/( cells[ind_h1-1].n*(*gas).m*pow( (*gas).kb*cells[ind_h1-1].T/(*gas).m, 1.5 ) );
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
MetroHasting(Nnew1_DSMC, lambda, 3, U_temp);
//MomentsofSamples(Nnew1, U_temp, mm, 3);
//printf("mm before %lf %lf %lf\n", mm[0], mm[1], mm[2]);
for(i=0; i<Nnew1_DSMC; i++){
	  U1new[j] = dis_norm(gen)*std1;
		U3new[j] = dis_norm(gen)*std1;
		U2new[j] = U_temp[i]*std1+m;
    x2new[j] = cells[ind_h1-1].dim[2] + (*box).delta_dim[1]*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h1-1].T;
    rhonew[j] = cells[ind_h1-1].n*(*gas).m;
		j++;
}


m = cells[ind_h2+1].U_space[1];
xx = 0.5*cells[ind_h2+1].Q[1]/( cells[ind_h2+1].n*(*gas).m*pow( (*gas).kb*cells[ind_h2+1].T/(*gas).m, 1.5 ) );
if( (*gas).gp_nn == GP ){
  evaluation_GPR(gpr, xx, lambda);
}
else{
  evaluation(NN, xx, lambda);
}
//evaluation_GPR(gpr, xx, lambda);
MetroHasting(Nnew3_DSMC, lambda, 3, U_temp);
//MomentsofSamples(Nnew1, U_temp, mm, 3);
//printf("mm before %lf %lf %lf\n", mm[0], mm[1], mm[2]);
for(i=0; i<Nnew3_DSMC; i++){
	  U1new[j] = dis_norm(gen)*std3;
		U3new[j] = dis_norm(gen)*std3;
		U2new[j] = U_temp[i]*std3+m;
    x2new[j] = cells[ind_h2+1].dim[2] + (*box).delta_dim[1]*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h2+1].T;
    rhonew[j] = cells[ind_h2+1].n*(*gas).m;
		j++;
}

free(U_temp);

/*
double stdT;
double VMP, SC, QA, FS1, FS2, U, pratio,FS1m;
VMP = sqrt(2.0*(*gas).kb*cells[ind_h1-1].T/(*gas).m);
SC = cells[ind_h1-1].U_space[1]/VMP;
QA = 3.0;
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));

stdT = sqrt( (*gas).kb*cells[ind_h1-1].T/(*gas).m );
for(i=0; i<Nnew1; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		}
		U2new[j] = (U+SC)*VMP;
    x2new[j] = cells[ind_h1].dim[2] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h1-1].T;
    rhonew[j] = cells[ind_h1-1].n*(*gas).m;
		j++;
}

VMP = sqrt(2.0*(*gas).kb*cells[ind_h2+1].T/(*gas).m);
SC = cells[ind_h2+1].U_space[1]/VMP;
//VMP = sqrt(2.0*(*gas).kb*cells[(*box).N[1]-1].T/(*gas).m);
//stdT = sqrt( (*gas).kb*cells[(*box).N[1]-1].T/(*gas).m );
stdT = sqrt( (*gas).kb*cells[ind_h2+1].T/(*gas).m );
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));
for(i=0; i<Nnew3; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1m-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1m-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		}
		U2new[j] = (U+SC)*VMP;
		x2new[j] = cells[ind_h2].dim[3] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = cells[ind_h2+1].T;
    rhonew[j] =cells[ind_h2+1].n*(*gas).m;
		j++;
}
*/

for(i=0; i<Nnew1_SPH; i++){
  x2new[j] = cells[ind_h1].dim[2] + (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew1_SPH);
  U1new[j] = cells[ind_h1].U_space[0];
  U2new[j] = cells[ind_h1].U_space[1];
  U3new[j] = cells[ind_h1].U_space[2];
  x1new[j] = (*x1)[0];
  x1_oldnew[j] = x1new[j];
  x2_oldnew[j] = x2new[j];
  colornew[j] = sph;
  Tnew[j] = cells[ind_h1].T;
  rhonew[j] = cells[ind_h1].n*(*gas).m;
  j++;
}

for(i=0; i<Nnew3_SPH; i++){
  x2new[j] = cells[ind_h2].dim[3] - (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew3_SPH);
  U1new[j] = cells[ind_h2].U_space[0];
  U2new[j] = cells[ind_h2].U_space[1];
  U3new[j] = cells[ind_h2].U_space[2];
  x1new[j] = (*x1)[0];
  x1_oldnew[j] = x1new[j];
  x2_oldnew[j] = x2new[j];
  colornew[j] = sph;
  Tnew[j] = cells[ind_h2].T;
  rhonew[j] = cells[ind_h2].n*(*gas).m;
  j++;
}




if(j!=Ntot){
  printf("---- in addrem_BC_hybrid j=%d Ntot=%d\n",j,Ntot);
}
/////  reallocating memory
free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
free(*color); free(*x1); free(*x1_old); free(*T); free(*rho);
(*gas).N = Ntot;
*x1 = x1new;
*x2 = x2new;
*T = Tnew;
*rho = rhonew;
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
}
}
}


void addrem_BC_hybrid(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho){
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::normal_distribution<> dis_norm(0.0,1.0);
//double sige = 0.9, alpha = 0.1;

int i,j,Nnew1, Nnew3;
double tol = 1e-15;
double Cratio, Hratio;
int NnewC,NnewH,N_new, Ntot;
int num;
double std1 = sqrt((*gas).kb*(*gas).Tv1/(*gas).m);
double std3 = sqrt((*gas).kb*(*gas).Tv3/(*gas).m);
double dt_remain,yw;
int Nout = 0;
double pi = acos(-1.0);
// particles that leave the domain needed to be removed
for(i=0; i<(*gas).N; i++){
	if( (*x2)[i] < tol ){
    if( (*color)[i] == sph ){
      (*x2)[i] = 2.0*tol;//(*x2)[i] + 2.0*fabs( (*x2)[i]  ) + tol;
      //(*T)[i] = (*gas).Tv1;
      (*U2)[i] = fabs( (*U2)[i] );
      (*flag)[i] = 1;
    }
    else{
      (*flag)[i] = -1;
      Nout++;
    }
  }
  else if( (*x2)[i] > (*box).Len[1]-tol){
    if( (*color)[i] == sph ){
      (*x2)[i] = (*box).Len[1] - 2.0*tol;
      //(*x2)[i] = (*x2)[i] - 2.0*fabs( (*x2)[i] - (*box).Len[1] ) - tol;
      //(*T)[i] = (*gas).Tv3;
      (*U2)[i] = -fabs( (*U2)[i] );
      (*flag)[i] = 1;
    }
    else{
      (*flag)[i] = -1;
      Nout++;
    }
  }
  else{
    (*flag)[i] = 1;
  }
}
// Now, how many new particles enter domain
double beta, dummy;
beta = 1.0/(std1*sqrt(2.0));
//beta = 1.0/(stdc);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nv1/(beta*2.0*sqrt(pi));
if(cells[0].model == dsmc)
  Nnew1 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;
else
  Nnew1 = 0;

beta = 1.0/(std3*sqrt(2.0));
//beta = 1.0/(stdh);
dummy = (*gas).delta_t*(*box).Len[0]*(*box).Len[2]*(*gas).nv3/(beta*2.0*sqrt(pi));
if(cells[(*box).N[1]-1].model == dsmc)
  Nnew3 = floor(dummy/(*gas).Fn+ dis_x(gen));//-count_reem_c;
else
  Nnew3 = 0;
//if( NnewC<2){
//	printf("\n*** NnewC<1, = %d\n\n",NnewC);
//}
//if( NnewH<2){
//	printf("\n*** NnewH<1, = %d\n\n",NnewH);
//}
N_new = Nnew1 + Nnew3;
if(N_new>0 || Nout>0){
//printf("**   N_new=%d and Nout=%d\n", N_new, Nout);
Ntot = N_new + (*gas).N - Nout;


double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew, *Tnew, *rhonew;
int *colornew;
//if(Ntot != (*gas).N){
  x1new = (double *) malloc( Ntot * sizeof(double) );
  x2new = (double *) malloc( Ntot * sizeof(double) );
  U1new = (double *) malloc( Ntot * sizeof(double) );
  U2new = (double *) malloc( Ntot * sizeof(double) );
  U3new = (double *) malloc( Ntot * sizeof(double) );
  x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
  x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
  colornew = (int *) malloc( Ntot * sizeof(int) );

  Tnew = (double *) malloc( Ntot * sizeof(double) );
  rhonew = (double *) malloc( Ntot * sizeof(double) );
//}
// store the surviving partices
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
   Tnew[j] = (*T)[i];
   rhonew[j] = (*rho)[i];
   j++;
 }
}

double stdT;
double VMP, SC, QA, FS1, FS2, U, pratio,FS1m;
VMP = sqrt(2.0*(*gas).kb*cells[0].T/(*gas).m);
SC = 0.0;//(*box).U_wall_5/VMP;
QA = 3.0;
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));

stdT = sqrt( (*gas).kb*cells[0].T/(*gas).m );
for(i=0; i<Nnew1; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
		}
		U2new[j] = (U+SC)*VMP;
    x2new[j] = U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = (*gas).Tv1;
    rhonew[j] = (*gas).nv1*(*gas).m;
		j++;
}

VMP = sqrt(2.0*(*gas).kb*cells[(*box).N[1]-1].T/(*gas).m);
stdT = sqrt( (*gas).kb*cells[(*box).N[1]-1].T/(*gas).m );
if(SC<-3.0)
	QA = fabs(SC)+1.0;
FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
FS1m = SC-sqrt(SC*SC+2.0);
FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));
for(i=0; i<Nnew3; i++){
	  U1new[j] = dis_norm(gen)*stdT;
		U3new[j] = dis_norm(gen)*stdT;
		U = (10.0*dis_x(gen)*FS1m-SC);
		pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		while(pratio<dis_x(gen)){
			U = (10.0*dis_x(gen)*FS1m-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1-U*U )/FS1m;
		}
		U2new[j] = (U+SC)*VMP;
		x2new[j] = (*box).Len[1] + U2new[j]*(*gas).delta_t*dis_x(gen);
    x1new[j] = (*x1)[0];
    x1_oldnew[j] = x1new[j];
    x2_oldnew[j] = x2new[j];
    colornew[j] = dsmc;
    Tnew[j] = (*gas).Tv3;
    rhonew[j] = (*gas).nv3*(*gas).m;
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

if(j!=Ntot){
  printf("---- in addrem_BC_hybrid j=%d Ntot=%d\n",j,Ntot);
}
/////  reallocating memory
free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
free(*color); free(*x1); free(*x1_old); free(*T); free(*rho);
(*gas).N = Ntot;
*x1 = x1new;
*x2 = x2new;
*T = Tnew;
*rho = rhonew;
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
}
}

void SPH_simple_hybrid(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *color){
  int N = (*gas).N, i, j, id;
  int search, kk, kkk, num;
	double W, sqrt_pi, mb, rr, r1r2, r2, r1, h, r_cut;
  double R = (*gas).kb/(*gas).m;
  double D2W, DW, k, mui, gama;
  double  mu0;
  double cv;
  cv = 3.0*(*gas).kb/(*gas).m/2.0;
  mu0 = (*gas).mu;// 5.0/( 16.0*(*gas).sigma*(*gas).sigma )*sqrt( (*gas).m*(*gas).kb*(*gas).T0/acos(-1.0) );
  gama = 1.4;
  r_cut = (*gas).r_cut;
  h = (*gas).h;
  mb = (*gas).mb;
  search = floor(r_cut/(*box).delta_dim[1])+1;
  sqrt_pi = sqrt(2.0*acos(-1.0));

  double *F = (double *) malloc( N * sizeof(double) );
  double *G = (double *) malloc( N * sizeof(double) );
  double *H = (double *) malloc( N * sizeof(double) );

  double *eps = (double *) malloc( N * sizeof(double) );
  //double *dmudx = (double *) malloc( N * sizeof(double) );
  //double *d2udx2 = (double *) malloc( N * sizeof(double) );
  //double *dkdx = (double *) malloc( N * sizeof(double) );
  //double *dTdx = (double *) malloc( N * sizeof(double) );
  //double *p = (double *) malloc( N * sizeof(double) );
  //double *dpdx = (double *) malloc( N * sizeof(double) );
  //double *d2Tdx2 = (double *) malloc( N * sizeof(double) );

  double  f_shear, f_pr, ki, kj;
	double c12, rho12, mu12, alpha,beta,eta2, v1v2, pi12;
	double muj, kij;
	alpha = 1.0;
	beta = 2.0;
	eta2 = pow( h*0.1, 2 );
  int N_dummy;
  double rhoj, U2j, Tj, x2j, epsj;
  for(i=0; i<N; i++){
    eps[i] = 0.0;
    //dmudx[i] = 0.0;
    //d2udx2[i] = 0.0;
    //dkdx[i] = 0.0;
    //dTdx[i] = 0.0;
    //d2Tdx2[i] = 0.0;
    //dpdx[i] = 0.0;
    num = index[i];
    if( cells[num].model == sph ){
      if( fabs( rho[i] ) < 1e-14 ){
        printf("rho[i=%d] = %e\n", i, rho[i]);
        exit(0);
      }
      r1 = x2[i];
      mui = mu0*pow(T[i]/(*gas).T0,0.5);
      ki = 15.0*(*gas).kb*mui/(4.0*(*gas).m);
      for(kk = num-search; kk<= num+search; kk++){
        if(kk<0)
           kkk = 0;
        else if(kk>(*box).N[1]-1)
           kkk = (*box).N[1]-1;
        else
           kkk = kk;
        if(cells[kkk].model == sph){
          for(id=0; id<cells[kkk].num_inside; id++){
            j = cells[kkk].indices_inside[id];
            if(color[j] == sph){
               if( fabs( rho[j] ) < 1e-14 ){
                 printf("rho[j=%d] = %e\n", j, rho[j]);
                 exit(0);
                }
                r2 = x2[j]+(kk-kkk)*(*box).delta_dim[1];
                r1r2 = r2-r1;
                rr = fabs(r2-r1);
                if(rr<r_cut){
                  W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                  D2W = W*( -1.0+rr*rr/h/h )/h/h;
                  if(i!=j){
                    DW = W*r1r2*(-1.0/h/h);
                  }
                  else{
                    DW = 0.0;// D2W = 0.0;
                  }
                  muj = mu0*pow(T[j]/(*gas).T0,0.5);
                  kj = 15.0*(*gas).kb*muj/(4.0*(*gas).m);

                  eps[i] += 4.0/3.0*mb/rho[j]*(U2[i]-U2[j])*DW;
                  //dudx[i] += mb*U2[j]/rho[j]*DW;
                  //dTdx[i] += mb*T[j]/rho[j]*DW;
                  //dmudx[i] += mb*muj/rho[j]*DW;
                  //dkdx[i] += mb*kj/rho[j]*DW;
                  //d2udx2[i] += mb*U2[j]/rho[j]*D2W;
                  //d2Tdx2[i] += mb*T[j]/rho[j]*D2W;
                  //dpdx[i] += mb*R*T[j]*DW;
              }
            }
          }
        }

        else{
          N_dummy = floor(cells[kkk].n*cells[kkk].volume/( (*gas).mb/(*gas).m ));//integer
          rhoj = cells[kkk].n*(*gas).m;
          if( fabs( rhoj ) < 1e-14 ){
            printf("rhoj = %e\n",  rhoj);
            exit(0);
           }
          U2j = cells[kkk].U_space[1];
          Tj = cells[kkk].T;
          for(j=0; j<N_dummy; j++){
            x2j = cells[kkk].cell_center[1] -(*box).delta_dim[1]/2.0 + (j+0.5)*(*box).delta_dim[1]/(1.0*N_dummy);
            r2 = x2j+(kk-kkk)*(*box).delta_dim[1];
            r1r2 = r2-r1;
            rr = fabs(r2-r1);
            if(rr<r_cut){
                  W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                  if(rr>1e-15){
                    DW = W*r1r2*(-1.0/h/h);
                  }
                  else{
                    DW = 0.0;
                  }
                  muj = mu0*pow(Tj/(*gas).T0,0.5);
                  //muj = mu0;
                  kj = 15.0*(*gas).kb*muj/(4.0*(*gas).m);

                  eps[i] += 4.0/3.0*mb/rhoj*(U2[i]-U2j)*DW;
            }
          }
        }

        if( eps[i] != eps[i] ){
          printf("eps[i] = %e\n",  eps[i]);
          exit(0);
         }
      }
    }
  }

  for(i=0; i<N; i++){
    F[i] = 0.0;
    G[i] = 0.0;
    H[i] = 0.0;
    num = index[i];
    if( cells[num].model == sph ){
      r1 = x2[i];
      mui = mu0*pow(T[i]/(*gas).T0,0.5);
      ki = 15.0*(*gas).kb*mui/(4.0*(*gas).m);
      for(kk = num-search; kk<= num+search; kk++){
        if(kk<0)
    	     kkk = 0;
    	  else if(kk>(*box).N[1]-1)
    	     kkk = (*box).N[1]-1;
    	  else
    	     kkk = kk;
        if(cells[kkk].model == sph){
    	    for(id=0; id<cells[kkk].num_inside; id++){
    		    j = cells[kkk].indices_inside[id];
            if(color[j] == sph){
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
                       DW = 0.0;// D2W = 0.0;
                     }
                     if(W!=W || DW!=DW || U2[j]!=U2[j] || T[j]!=T[j] || rho[j] != rho[j] || rho[i] != rho[i]){
                       num = index[i];
                       printf("W!=W || DW!=DW\n");
                       printf("rho[i]=%e and rho[j] = %e\n", rho[i], rho[j]);
                       exit(0);
                     }

								     muj = mu0*pow(T[j]/(*gas).T0,0.5);
                     kj = 15.0*(*gas).kb*muj/(4.0*(*gas).m);
								     v1v2 = U2[j]-U2[i];


								             //if( v1v2*r1r2 < 0.0 ){
								         //	c12 = 0.5*( sqrt(gama*(*gas).kb*T[i]/(*gas).m) + sqrt(gama*(*gas).kb*T[j]/(*gas).m) );
								        //	rho12 = 0.5*(rho[i]+rho[j]);
								         //	mu12 = h*v1v2*r1r2/(r1r2*r1r2+eta2);
                        //  mu12 = -16.0*mui*muj/(mui+muj)*v1v2*r1r2/(r1r2*r1r2+eta2);
                        //  //mu12 = mui*muj/(mui+muj);
								        //	//pi12 = (-alpha*c12*mu12+beta*mu12*mu12)/rho12;
                        //    pi12 = mu12/rho12;
							               //	}
						               //		else{
						                //			pi12 = 0.0;
						                      //		}

								                          //pi12=0.0;
                                          //F[i] +=  mb*(R*T[i]/rho[i]+R*T[j]/rho[j] + pi12)*DW;
                                          //F[i] +=  mb*mu*U2[j]/rho[j]*D2W;
                                          //G[i] +=  0.5*mb*R*( T[i]/rho[i] + T[j]/rho[j] )*( U2[j]-U2[i] )*DW;
								                                  //G[i] +=  0.5*mb*( R*T[i]/rho[i] + R*T[j]/rho[j]  + pi12)*( U2[j]-U2[i] )*DW;

								       f_pr = ( rho[i]*R*T[i]+rho[j]*R*T[j] )/(rho[i]*rho[j]);//+R*T[j]/rho[j]);//*DW;
                       //f_pr =  R*T[j]/rho[i]*DW;
								               //f_shear = pi12*DW;//0.5*mu0*(U2[j]-U2[i])/rho[j]*D2W/rho[i];
								       //f_shear = muj*(U2[j]-U2[i])/rho[j]*D2W/rho[i];
								               //f_shear = 0.0;//4.0/3.0*mui*(U2[j]-U2[i])/rho[j]*D2W/rho[i];
                               //f_shear = 4.0/3.0*( mui*d2udx2[i] + dmudx[i]*dudx[i])/rho[i];
                        f_shear = -(mui*eps[i]+muj*eps[j])/(rho[i]*rho[j]);
							                 //if(fabs(U2[j])>1e-1){
								                       //	printf("here\n");
								                               //}
                        if(f_pr!=f_pr || f_shear!=f_shear){
                            printf("f_pr!=f_pr || f_shear!=f_shear\n");
                            printf("rho[i]=%e and rho[j] = %e\n", rho[i], rho[j]);
                            printf("T[i]=%e and T[j] = %e\n", T[i], T[j]);
                            printf("eps[i]=%e and eps[j] = %e\n", eps[i], eps[j]);
                            exit(0);
                        }
                        //kij = 4.0*ki*kj/(ki+kj);
                        kij = 0.5*(ki+kj);
								        F[i] +=  mb*(f_pr+f_shear)*DW;
								        G[i] +=  -0.5*mb*( f_pr*(U2[i]-U2[j]) )*DW + mb*kij*(T[i]-T[j])*DW/rho[i]/rho[j];
                        //G[i] +=  + mb*kij*(T[i]-T[j])*DW/rho[i]/rho[j];
                        H[i] += mb*(U2[i]-U2[j])*DW/rho[j]*rho[i];
                 }
             }
    		   }
          }

          else{
            N_dummy = floor(cells[kkk].n*cells[kkk].volume/( (*gas).mb/(*gas).m ));//integer
            rhoj = cells[kkk].n*(*gas).m;
            U2j = cells[kkk].U_space[1];
            Tj = cells[kkk].T;
            if(Tj<1e-16){
              printf("Tj = %e\n", Tj);
              printf("rhoj = %e\n", rhoj);
              printf("U2j = %e\n", U2j);
              printf("N_dummy = %d\n", N_dummy);

              exit(0);
            }
            muj = mu0*pow(Tj/(*gas).T0,0.5);
            //muj = mu0;
            kj = 15.0*(*gas).kb*muj/(4.0*(*gas).m);
            epsj = -(cells[kkk+1].U_space[1]-cells[kkk-1].U_space[1])/(*box).delta_dim[1]/2.0*4.0/3.0;
            //epsj =  cells[kkk].PIJ[3] /muj*mb/rhoj*0.5;//*mb/rhoj;//cells[kkk].n*(*gas).m*0.5;
            for(j=0; j<N_dummy; j++){
              x2j = cells[kkk].cell_center[1] -(*box).delta_dim[1]/2.0 + (j+0.5)*(*box).delta_dim[1]/(1.0*N_dummy);
              r2 = x2j+(kk-kkk)*(*box).delta_dim[1];
              r1r2 = r2-r1;
              rr = fabs(r2-r1);
              if(rr<r_cut){
                    W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                    if(rr>1e-15){
                      DW = W*r1r2*(-1.0/h/h);
                    }
                    else{
                      DW = 0.0;
                    }
                    if(W!=W || DW!=DW){
                      num = index[i];
                      printf("W!=W || DW!=DW\n");
                      printf("rho[i]=%e and rhoj = %e\n", rho[i], rhoj);
                      exit(0);
                    }

                    v1v2 = U2j-U2[i];
                    f_pr = ( rho[i]*R*T[i]+rhoj*R*Tj )/(rho[i]*rhoj);
                    f_shear = -(mui*eps[i]+muj*epsj)/(rho[i]*rhoj);
                    //f_shear = 0.0;//-(mui*eps[i]+muj*epsj)/(rho[i]*rhoj);
                    //f_shear = -(mui*eps[i])/(rho[i]);
                    //kij = 4.0*ki*kj/(ki+kj);
                    kij = 0.5*(ki+kj);
                    F[i] +=  mb*(f_pr+f_shear)*DW;
                    G[i] +=  -0.5*mb*( f_pr*(U2[i]-U2j) )*DW + mb*kij*(T[i]-Tj)*DW/rho[i]/rhoj;
                    H[i] += mb*(U2[i]-U2j)*DW/rhoj*rho[i];
              }
            }
          }

        }
      //f_shear = 4.0/3.0*( mui*d2udx2[i] + dmudx[i]*dudx[i]);
      //F[i] =  (-dpdx[i]+4.0/3.0*( mui*d2udx2[i] + dmudx[i]*dudx[i]))/rho[i];
      //F[i] =  (-dpdx[i] +4.0/3.0*( mui*d2udx2[i]) )/rho[i];
      //G[i] =  (-rho[i]*R*T[i]*dudx[i] + 4.0/3.0*mui*d2udx2[i] + ki*d2Tdx2[i]+dkdx[i]*dTdx[i])/rho[i];



      G[i] += mui*eps[i]*eps[i]/2.0/rho[i];
    }
  }

  for(i=0; i<N; i++){
      U2[i] = U2[i] + F[i]*(*gas).delta_t;
      T[i]  = T[i]  + G[i]*(*gas).delta_t/cv*0.5;
      if(T[i]<1e-16){
        printf("T[i] = %e\n", T[i]);
        printf("num = %d\n", index[i]);
        exit(0);
      }
      //rho[i] = rho[i] + H[i]*(*gas).delta_t;
      if(U2[i]!=U2[i] || T[i]!=T[i]){
        num = index[i];
        printf("U2[%d]=%e, T[%d]=%e with model[%d]=%d while  model[%d]=%d\n", i, U2[i], i, T[i], num, cells[num].model, num+1,cells[num+1].model);
      }
  }



  //free(dudx); free(dmudx); free(d2udx2); free(dkdx); free(dTdx); free(d2Tdx2);  free(dpdx);
  free(eps);
  free(F); free(G); free(H);
}



void Collision_DSMC_VHS_hybrid(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC, int *color){
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
	N_cells = (*box).N[1];
  for( (*DSMC).num = 0; (*DSMC).num<N_cells; (*DSMC).num++){
    if( cells[ (*DSMC).num ].model != sph ){
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
              if( color[ (*DSMC).pair[0] ] == sph || color[ (*DSMC).pair[1] ] == sph){
                printf("pair0: color[%d]=%d\n", (*DSMC).pair[0], color[ (*DSMC).pair[0] ]);
                printf("pair1: color[%d]=%d\n", (*DSMC).pair[1], color[ (*DSMC).pair[1] ]);
              }
              if( color[ (*DSMC).pair[0] ] == dsmc &&  color[ (*DSMC).pair[1] ] == dsmc)
              {
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
         	      if (sigmat/cells[(*DSMC).num].crm> dis_x(gen))
                {
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
	       }
         //if(num_act_col < 5){
		     //    printf("\n************            number of actual collisions = %d  at cell = %d step = %d  M_can=%e     *****************\n", num_act_col, (*DSMC).num, (*box).step, (*DSMC).M_cand);
		     //    printf("Nc = %d, omega_max = %e, omega=%e",(*DSMC).Nc,cells[(*DSMC).num].crm,crm);
	       //}
   	     cells[(*DSMC).num].crm = max( crm, cells[(*DSMC).num].crm );
      }
   }
  }
}



void cell_update_hybrid(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double * T, double * rho, int *color)
{
  int Ncells = (*box).N[2]*(*box).N[1]*(*box).N[0];
  int num, i, j, id;
  double Mp1, Mp2, Mp3, energy;
  int *numbering = (int *) malloc(  Ncells* sizeof(int) );
  double trace;
  // initialize particles and cells to be sph for zeroth time step
  if( (*box).step == 0 ){
    for(i=0; i<(*gas).N; i++){
      if((*gas).model=="SPH" || (*gas).model=="Hybrid"){
        color[i] = sph;
      }
      else if( (*gas).model=="DSMC_VHS" ){
        color[i] = dsmc;
      }
      if(x2[i]<0.5*(*box).Len[1]){
        //rho[i] = (*gas).m*(*gas).nv1;
        T[i] = (*gas).Tv1;
      }
      else{
        //rho[i] = (*gas).m*(*gas).nv3;
        T[i] = (*gas).Tv3;
      }
    }
    for(num=0; num<(*box).N[1]; num++){
      if((*gas).model=="SPH" || (*gas).model=="Hybrid"){
        cells[num].model = sph;
      }
      else if( (*gas).model=="DSMC_VHS" ){
        cells[num].model = dsmc;
      }
    }
  }

  double *ws = (double *) malloc(  Ncells* sizeof(double) );

  // save old values
  double *Told = (double *) malloc(  Ncells* sizeof(double) );
  double *nold = (double *) malloc(  Ncells* sizeof(double) );
  double *U1old = (double *) malloc(  Ncells* sizeof(double) );
  double *U2old = (double *) malloc(  Ncells* sizeof(double) );
  double *U3old = (double *) malloc(  Ncells* sizeof(double) );
  for(num=0; num<Ncells; num++){
    Told[num] = cells[num].T;
    nold[num] = cells[num].n;
    U1old[num] = cells[num].U_space[0];
    U2old[num] = cells[num].U_space[1];
    U3old[num] = cells[num].U_space[2];
  }


  for (num=0; num < Ncells; num++){
    cells[num].num_inside = 0;
    numbering[num] = 0;
  }

  for (i=0; i < (*gas).N; i++){
      num = floor(x2[i]/(*box).delta_dim[1]);
      if(num<0){
          printf("OOPS, particle is out from down, x2=%e x2<0, color=%d\n",x2[i], color[i]);
          exit(1);
      }
      if(num>(*box).N[1]-1){
          printf("OOPS, particle is out from up, x2=%e x2>L, color=%d\n",x2[i], color[i]);
          exit(1);
      }
      index[i] = num;
      //if(cells[num].model != hybrid || color[i] != sph){
        cells[num].num_inside = cells[num].num_inside+1;
      //}
  }
  for (num=0; num<Ncells; num++){
   free(cells[num].indices_inside);
   cells[num].indices_inside = (int *) malloc( cells[num].num_inside * sizeof(int) );
   }
   for (i=0; i < (*gas).N; i++){
     num = index[i];
     //if(cells[num].model != hybrid || color[i] != sph){
       cells[num].indices_inside[ numbering[num] ] = i;
       numbering[num]++;
     //}
   }



   for(num=0; num< Ncells; num++){
     cells[num].n_smooth = 0.0;
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

   double w = 1.0;


   for(num=0; num< Ncells; num++){
     if(cells[num].model != sph){
       for(i=0; i<cells[num].num_inside; i++){
         id = cells[num].indices_inside[i];
         if(color[id]!=sph){
           cells[num].U_space[0] = cells[num].U_space[0]+U1[id]*w;
           cells[num].U_space[1] = cells[num].U_space[1]+U2[id]*w;
           cells[num].U_space[2] = cells[num].U_space[2]+U3[id]*w;
           cells[num].weight += w;
         }
       }
     }
   }
   for(num=0; num<Ncells; num++){
     if(cells[num].model != sph){
       if(cells[num].num_inside >0){
 	    	  cells[num].U_space[0] = cells[num].U_space[0]/cells[num].weight;
     		  cells[num].U_space[1] = cells[num].U_space[1]/cells[num].weight;
     		  cells[num].U_space[2] = cells[num].U_space[2]/cells[num].weight;

          U2old[num] = cells[num].U_space[1];
        }
     }
   }
   //#pragma omp parallel for  private(i, id, energy)
   for(num=0; num< Ncells; num++){
     if(cells[num].model != sph){
       for(i=0; i<cells[num].num_inside; i++){
         id = cells[num].indices_inside[i];
         if(color[id]!=sph){
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
     }
   }
   for(num=0; num<Ncells; num++){
     if(cells[num].model != sph){
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
       cells[num].n = (*gas).Fn*cells[num].weight/cells[num].volume;
       nold[num] = cells[num].n;
       trace = cells[num].PIJ[0]+cells[num].PIJ[3]+cells[num].PIJ[5];
       for(j=0; j<3; j++){
         cells[num].Q[j] = cells[num].Q[j]*cells[num].n*(*gas).m/2.0;

         cells[num].PIJ[j] = cells[num].PIJ[j] - (1.0/3.0)*trace;
         cells[num].PIJ[j] = cells[num].PIJ[j]*cells[num].n*(*gas).m;
       }
     }
   }

   int search, kk, kkk;
 	 double W, sqrt_pi, mb, rr, r1r2, r2, r1, h, r_cut, s2h2;
   double mu,k,DW,dudx,dTdx;
 	 double cv = 3.0*(*gas).kb/(*gas).m/2.0;
   //double gama = 1.4;
   //r_cut = 2.0*(*box).delta_dim[1];// working for 1e-5
   //h = r_cut/30.0;
   r_cut = 2.0*(*box).delta_dim[1];
   h = r_cut/30.0;

   (*gas).r_cut = r_cut; (*gas).h = h;
   search = floor(r_cut/(*box).delta_dim[1])+1;
   mb = (*gas).m*(*gas).Lv1/(*gas).Nv1*(*gas).nv1;
   (*gas).mb = mb;
   sqrt_pi = sqrt(2.0*acos(-1.0));
   //if((*box).step==0){
 	   for(i=0; i<(*gas).N; i++){
       if( color[i] == sph)
 		     rho[i] = 0.0;
 	   }
   //}



   for(num=0; num<Ncells; num++){
 	   if( cells[num].model == sph ){// if it was SPH
 		    cells[num].PIJ[3] = 0.0;
 			  cells[num].Q[1] = 0.0;
        cells[num].n = 0.0;
        cells[num].T = 0.0;
        for(i=0; i<3; i++){
          cells[num].U_space[i] = 0.0;
        }
 		 }
 	}

  // compute rho_i for particles in the sph cells
  double Tj, rhoj, U1j, U2j, U3j, x2j;
  int N_dummy;
//if((*box).step==0){
  for(i=0; i<(*gas).N; i++){
		r1 = x2[i];
		num = index[i];
    if(cells[num].model == sph){
		    for(kk = num-search; kk<= num+search; kk++){
			     if(kk<0)
				       kkk = 0;
			     else if(kk>(*box).N[1]-1)
				       kkk = (*box).N[1]-1;
			     else
				       kkk = kk;
           if(cells[kkk].model == sph){
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
           else{
             N_dummy = floor(nold[kkk]*cells[kkk].volume/( (*gas).mb/(*gas).m ));//integer
             rhoj = nold[kkk]*(*gas).m;
             for(j=0; j<N_dummy; j++){
               x2j = cells[kkk].cell_center[1] -(*box).delta_dim[1]/2.0 + (j+0.5)*(*box).delta_dim[1]/(1.0*N_dummy);
               r2 = x2j+(kk-kkk)*(*box).delta_dim[1];
               r1r2 = r2-r1;
               rr = fabs(r2-r1);
               if(rr<r_cut){
                  W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                  rho[i] += mb*W;
               }
             }
           }
			 }
		}
	}
//}
  // fix rho, T of dummy sph particles in hybrid cells to mean values
  /*
  for(i=0; i<(*gas).N; i++){
    num = index[i];
    if(cells[num].model == hybrid && color[i] == sph){
      rho[i] = (*gas).m*cells[num].n;
      T[i] = cells[num].T;
    }
  }
  */
  // compute cell values of tau and q, n, T for sph cell

/*
  for(i=0; i<(*gas).N; i++){
    num = index[i];
    if(cells[num].model == sph && color[i] == sph){
      cells[num].T += T[i];
      //cells[num].U_space[0] += U1[i];
      //cells[num].U_space[1] += U2[i];
      //cells[num].U_space[2] += U3[i];
      cells[num].n  += rho[i];
    }
  }
  for(num=0; num<(*box).N[1]; num++){
      if(cells[num].model == sph){
        cells[num].T = cells[num].T/(1.0*cells[num].num_inside);
        cells[num].n = cells[num].n/(1.0*cells[num].num_inside)/(*gas).m;
        //for(j=0; j<3; j++){
        //  cells[num].U_space[j] = cells[num].U_space[j]/(1.0*cells[num].num_inside);
        //}
      }
  }
*/


  for(num=0; num<(*box).N[1]; num++){
     if(cells[num].model == sph){
	       r1 = cells[num].cell_center[1];
	       for(kk = num-search; kk<= num+search; kk++){
		         if(kk<0)
			          kkk = 0;
		         else if(kk>(*box).N[1]-1)
			          kkk = (*box).N[1]-1;
		         else
			          kkk = kk;
             if(cells[kkk].model == sph){
		           for(id=0; id<cells[kkk].num_inside; id++){
			            j = cells[kkk].indices_inside[id];
			            r2 = x2[j]+(kk-kkk)*(*box).delta_dim[1];
			            r1r2 = r2-r1;
			            rr = fabs(r2-r1);
			            if(rr<r_cut){
				                W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                        if(rr>1e-15){
                          DW = W*r1r2*(-1.0/h/h);
                        }
                        else{
                          DW = 0.0;
                        }
                        mu = (*gas).mu*pow(T[j]/(*gas).T0,0.5);
                        k = (*gas).k/(*gas).mu*mu;
                        cells[num].n += mb*W/(*gas).m;
                        dudx   = mb*U2[j]/rho[j]*DW;
                      	dTdx   = mb*T[j]/rho[j]*DW;
            						cells[num].PIJ[3] += 4.0/3.0*mu*dudx;
                      	cells[num].Q[1] += -k*dTdx;
                        cells[num].T += mb*T[j]*W/rho[j];
                        cells[num].U_space[0] += 0.0;//mb*U1[j]*W/rho[j];
                        cells[num].U_space[1] += mb*U2[j]*W/rho[j];
                        cells[num].U_space[2] += 0.0;//mb*U3[j]*W/rho[j];
                  }
               }
             }
             else{
               N_dummy = floor(nold[kkk]*cells[kkk].volume/( (*gas).mb/(*gas).m ));//integer
               Tj = Told[kkk];
               rhoj = nold[kkk]*(*gas).m;
               U1j = U1old[kkk];
               U2j = U2old[kkk];
               U3j = U3old[kkk];
               //printf("Tj = %e\n", Tj);
               //printf("rhoj = %e\n", rhoj);
               //printf("U1j = %e\n", U1j);
               //printf("U2j = %e\n", U2j);
               //printf("U3j = %e\n", U3j);
               //printf("N_dummy = %d\n", N_dummy);
               for(j=0; j<N_dummy; j++){
                 x2j = cells[kkk].cell_center[1] -(*box).delta_dim[1]/2.0 + (j+0.5)*(*box).delta_dim[1]/(1.0*N_dummy);
                 r2 = x2j+(kk-kkk)*(*box).delta_dim[1];
                 r1r2 = r2-r1;
                 rr = fabs(r2-r1);
                 if(rr<r_cut){
                       W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
                       if(rr>1e-15){
                         DW = W*r1r2*(-1.0/h/h);
                       }
                       else{
                         DW = 0.0;
                       }
                       mu = (*gas).mu*pow(Tj/(*gas).T0,0.5);
                       k = (*gas).k/(*gas).mu*mu;
                       cells[num].n += mb*W/(*gas).m;
                       dudx   = mb*U2j/rhoj*DW;
                       dTdx   = mb*Tj/rhoj*DW;
                       cells[num].PIJ[3] += 4.0/3.0*mu*dudx;
                       cells[num].Q[1] += -k*dTdx;
                       cells[num].T += mb*Tj*W/rhoj;
                       cells[num].U_space[0] += mb*U1j*W/rhoj;
                       cells[num].U_space[1] += mb*U2j*W/rhoj;
                       cells[num].U_space[2] += mb*U3j*W/rhoj;
                 }
               }
             }
			   }
		 }
	}
/*
for(num=0; num<(*box).N[1]; num++){
  ws[num] = 0.0;
  if(cells[num].model == dsmc){
    cells[num].T = 0.0;
  }
}
*/
int num2;
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
        if(cells[num2].model == sph){
          mb = (*gas).mb/(*gas).m;
        }
        else{
          mb = (*gas).Fn;
        }
        for(i=0; i<cells[num2].num_inside; i++){
          id = cells[num2].indices_inside[i];
          r2 = x2[id]+(kk-kkk)*(*box).delta_dim[1];
          r1r2 = r2-r1;
          rr = fabs(r2-r1);
          if(rr<r_cut){
            W = 1.0/sqrt_pi/h*exp(-rr*rr/h/h/2.0);
            cells[num].n_smooth += mb*W;
            /*
            if(cells[num].model == dsmc){
              ws[num] += W;
              Mp1 = U1[i]-cells[num2].U_space[0];
              Mp2 = U2[i]-cells[num2].U_space[1];
              Mp3 = U3[i]-cells[num2].U_space[2];
              energy = Mp1*Mp1 + Mp2*Mp2 + Mp3*Mp3;
              cells[num].T += (*gas).m*energy/(3.0*(*gas).kb)*W;
              //printf("T = %e\n", (*gas).m*energy/(3.0*(*gas).kb)*W*mb);
            }
            */
          }
        }
      }
    }
  }

/*
  for(num=0; num<(*box).N[1]; num++){
    cells[num].n = cells[num].n_smooth;
    //if(cells[num].model == dsmc){
    //  cells[num].T = cells[num].T/ws[num];
    //  printf("T = %e\n", cells[num].T);
    //}
  }
*/
  free(numbering);
  free(Told); free(U1old); free(U2old); free(U3old); free(nold);
  /*
  /// checking dummy sph particles in hybrid cells
  int count_all_dummy = 0, count_comp_part = 0, count_dsmc_hybrid=0;
  for(i=0; i<(*gas).N; i++){
    num = floor(x2[i]/(*box).delta_dim[1]);
    if(cells[num].model==hybrid && color[i]==sph){
      count_all_dummy++;
    }
  }
  for(num=0; num<(*box).N[1]; num++){
      count_comp_part += cells[num].num_inside;
  }
  printf("count_all_dummy=%d, and gas.N-sum_num_inside=%d\n", count_all_dummy, (*gas).N-count_comp_part);

  // printing index of dummy particles and DSMC particles in the hybrid cells
  printf("dummy sph particles in hybrid cells\n" );
  for(i=0; i<(*gas).N; i++){
    num = floor(x2[i]/(*box).delta_dim[1]);
    if(cells[num].model == hybrid && color[i]==sph){
      printf("%d ", i);
    }
  }
  printf("dsmc particles in hybrid cells\n" );
  for(num=0; num<(*box).N[1]; num++){
    if(cells[num].model == hybrid){
      for(id = 0; id<cells[num].num_inside; id++){
        i = cells[num].indices_inside[id];
        printf("%d ", i);
        count_dsmc_hybrid++;
      }
    }
  }
  printf("\n **  dsmc in hybrid = %d while sph in hybrid is %d\n",count_dsmc_hybrid, count_all_dummy);
  */
}
void fix_mean_variance(double *U_temp, int N){
  double  variance=0.0;
  double mean = 0.0;
  double lamb;
  int i;
  for(i=0; i<N;++i)
  {
      mean += U_temp[i];
  }
  mean=mean/(1.0*N);
  for(i=0; i<N;++i){
    variance+= pow(U_temp[i]-mean,2.0);
  }
  variance = variance/(1.0*N);
  for(i=0; i<N; i++){
    U_temp[i] = (U_temp[i]-mean)*sqrt(1.0/variance);
  }
}
void fix_mean(double *U_temp, int N){
  double mean = 0.0;
  double lamb;
  int i;
  for(i=0; i<N;++i)
  {
      mean += U_temp[i];
  }
  mean=mean/(1.0*N);
  for(i=0; i<N; i++){
    U_temp[i] = (U_temp[i]-mean);
  }
}

void add_remove_hybrid(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho, struct NerualNetwork *NN,  struct GPR *gpr){
  double R = (*gas).kb/(*gas).m;
  double q,phi,tau,trace,epsilon;
  int model;
  double lambda[3];

  // working epsilon
  epsilon = 8e-5;
  //epsilon = 0.001;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_x(-1.0, 1.0);
  std::normal_distribution<> dis_norm(0.0,1.0);

  // add remove particles after each time step
  // flag -1 the particles in the wrong cells
  // also flag -1 the particles in the DSMC-SPH or SPH-DSMC
  // models: 1. SPH, 2. SPH-DSMC, 3. DSMC-SPH, 4. DSMC
  int NremSPH = 0, NremDSMC = 0;
  int Ntot;
  int i, j, num;
  int *Nnew_SPH = (int *) malloc( (*box).N[1] * sizeof(int) );
  int *Nnew_Hybrid = (int *) malloc( (*box).N[1] * sizeof(int) );
  int *Nnew_DSMC = (int *) malloc( (*box).N[1] * sizeof(int) );
  for(num=0; num<(*box).N[1]; num++){
    Nnew_SPH[num] = 0; Nnew_DSMC[num] = 0; Nnew_Hybrid[num] = 0;
  }
  int ind_h1 = -1, ind_h2 = -1;
  for(num=0; num<(*box).N[1]; num++){
    if(cells[num].model == hybrid && ind_h1 == -1){
      ind_h1 = num;
    }
    else if(cells[num].model == hybrid && ind_h1 != -1){
      ind_h2 = num;
    }
  }
  if(ind_h1 != -1 && ind_h2 == -1){
    ind_h2 = ind_h1;
  }
  //if(ind_h1!=-1 || ind_h2!=-1){
  //  printf("ind_h1=%d ind_h2=%d\n", ind_h1, ind_h2);
  //}
  // before here, cells are either dsmc, sph, or hybrid(dsmc)
  for(num=0; num<(*box).N[1]; num++){
    if(cells[num].model == hybrid){
      cells[num].model = dsmc;
    }
  }
  for(i=0; i<(*gas).N; i++){
    (*flag)[i] = 1;
  }

  // check cell models which must change
  double RT, xx;
  double Z, p1, p2;
  for(num=0; num<(*box).N[1]; num++){
    RT = R*cells[num].T;
    q = cells[num].Q[1];
    tau = cells[num].PIJ[3];
    //phi = 1.0/(cells[num].n*(*gas).kb*cells[num].T);
    //phi = phi*sqrt( 2.0/5.0*q*q/R/cells[num].T + 0.5*tau*tau );
    //phi = pow(1.0/RT,1.5)*sqrt(q*q+tau*tau*RT)/(cells[num].n*(*gas).m);
    //printf("phi=%e\n", phi);
    xx = 0.5*q/( cells[num].n*(*gas).m*pow( (*gas).kb*cells[num].T/(*gas).m, 1.5 ) );
    if( fabs(xx)>10.0 ){
      printf("q=%lf, q/p=%lf\n", q,xx);
    }
    //evaluation(NN, xx, lambda);
    if( (*gas).gp_nn == GP ){
      evaluation_GPR(gpr, xx, lambda);
    }
    else{
      evaluation(NN, xx, lambda);
    }
    //lambda[0]=0.0; lambda[1]=0.5; lambda[2]=0.0;
    //phi = compute_DKL_MED_MW(lambda);
    phi = compute_fisher_metric(lambda);
    cells[num].xx = phi;
    //printf("num %d, phi=%lf\n",num, phi);
    //if( phi<epsilon &&   cells[num].model == dsmc  ){
    //  cells[num].model = dsmc_to_sph;
    //}
    if( phi>epsilon &&  cells[num].model == sph){
      printf("--->  step %d phi=%lf at num=%d \n", (*box).step, phi, num);
      cells[num].model = sph_to_dsmc;
    }

    //// for testing only
    //if( num >= floor( (*box).N[1]/2 )   && (*box).step == 1){
    //  cells[num].model = sph_to_dsmc;
    //}
  }
  /*
  for(num=1; num<(*box).N[1]; num++){
    if(cells[num].model == sph_to_dsmc && cells[num-1].model == sph){
      if(num>=1){
        cells[num-1].model = sph_to_dsmc;
      }
      if(num>=2 &&  cells[num-2].model == sph){
        cells[num-2].model = sph_to_dsmc;
      }
    }
  }
  */
  // removing particles in the wrong cell
  for(i=0; i<(*gas).N; i++){
    num = floor( (*x2)[i]/(*box).delta_dim[1] );
    if( (*color)[i] == sph && (cells[num].model != sph) ){
      (*flag)[i] = -1; NremSPH++;
    }
    else if( (*color)[i] == dsmc && (cells[num].model != dsmc ) ){
      /*
      if(ind_h1 !=-1 && num < ind_h1){
        (*x2)[i] = cells[ind_h1].dim[2] + fabs((*x2)[i]-cells[ind_h1].dim[2]);
        (*U2)[i] = -(*U2)[i];
      }
      else if(ind_h2 !=-1 && num > ind_h2){
        (*x2)[i] = cells[ind_h2].dim[3] - fabs((*x2)[i]-cells[ind_h2].dim[3]);
        (*U2)[i] = -(*U2)[i];
      }
      */
      (*flag)[i] = -1; NremDSMC++;
    }
  }

  /*
  printf("index of particles being removed\n" );
  for(i=0; i<(*gas).N; i++){
    if( (*flag)[i]==-1 ){
      printf("%d ", i);
    }
  }
  printf("\n");
  */
    // count how many of each type needs to be generated
  int N_new_SPH=0, N_new_DSMC= 0,N_new_Hybrid= 0;
  for(num=0; num<(*box).N[1]; num++){
    if( cells[num].model== sph_to_dsmc ){
      Nnew_DSMC[num] = floor( cells[num].n_smooth*cells[num].volume/(*gas).Fn+ dis_x(gen));
      N_new_DSMC += Nnew_DSMC[num];
    }
    if( cells[num].model== dsmc_to_sph ){
      Nnew_SPH[num] = ceil( cells[num].n*cells[num].volume/( (*gas).mb/(*gas).m) + dis_x(gen));
      N_new_SPH += Nnew_SPH[num];
    }
  }
  /*
  // find hybrid cells & count total dummy SPH particles needed there
 for(num=1; num<(*box).N[1]-1; num++){
   if( (cells[num].model==dsmc || cells[num].model==sph_to_dsmc) &&
       (cells[num-1].model==sph || (cells[num-1].model==dsmc_to_sph ||
          cells[num+1].model==sph || cells[num+1].model==dsmc_to_sph) ) ){
            //cells[num].model==hybrid;
            Nnew_Hybrid[num] = cells[num].n*cells[num].volume/( (*gas).mb/(*gas).m );
            N_new_Hybrid += Nnew_Hybrid[num];
  }
 }
 */
  Ntot = (*gas).N + N_new_SPH + N_new_DSMC - NremSPH - NremDSMC + N_new_Hybrid;
  // resize the vectors we have
  double *x2new, *U1new, *U2new, *U3new, *x2_oldnew, *x1new, *x1_oldnew, *Tnew, *rhonew;
  int *colornew;
  //if(Ntot != (*gas).N){
    x1new = (double *) malloc( Ntot * sizeof(double) );
 	  x2new = (double *) malloc( Ntot * sizeof(double) );
 	  U1new = (double *) malloc( Ntot * sizeof(double) );
 	  U2new = (double *) malloc( Ntot * sizeof(double) );
 	  U3new = (double *) malloc( Ntot * sizeof(double) );
 	  x1_oldnew = (double *) malloc( Ntot * sizeof(double) );
 	  x2_oldnew = (double *) malloc( Ntot * sizeof(double) );
 	  colornew = (int *) malloc( Ntot * sizeof(int) );

    Tnew = (double *) malloc( Ntot * sizeof(double) );
    rhonew = (double *) malloc( Ntot * sizeof(double) );
  //}
  // store the surviving partices
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
     Tnew[j] = (*T)[i];
     rhonew[j] = (*rho)[i];
     j++;
   }
  }
  //printf("j after copy = %d, gas.N-remv=%ld\n", j, (*gas).N- NremSPH - NremDSMC);

  int N = max_int(Nnew_DSMC,(*box).N[1]);
  double *U_temp;
  double *U1_temp, *U3_temp;
  if( N>0 ){
    U_temp = (double *) malloc( N * sizeof(double) );
    U1_temp = (double *) malloc( N * sizeof(double) );
    U3_temp = (double *) malloc( N * sizeof(double) );
  }
  double std;
  double mean,stdev;
  double d0, d1, d2, d3;
  int first = 0;
  double m[3];
  double Rho;
  for(num=0; num<(*box).N[1]; num++){
    // sample dsmc particles in sph_to_dsmc cells
    if( (cells[num].model==sph_to_dsmc || cells[num].model==hybrid) && Nnew_DSMC[num]>0 ){
       printf("--- generating %d DSMC particles in  num=%d \n",Nnew_DSMC[num], num);
       xx = 0.5*cells[num].Q[1]/( cells[num].n*(*gas).m*pow( (*gas).kb*cells[num].T/(*gas).m, 1.5 ) );
       //evaluation(NN, xx, lambda);
       if( (*gas).gp_nn == GP ){
         evaluation_GPR(gpr, xx, lambda);
       }
       else{
         evaluation(NN, xx, lambda);
       }
       //evaluation_GPR(gpr, xx, lambda);
       //lambda[0] = 1.6190676471471449e-01;
       //lambda[1] = 5.8465463351464497e-01;
       //lambda[2] = -5.5633305061260560e-02;
       MetroHasting(Nnew_DSMC[num], lambda, 3, U_temp);

       //RT = R*cells[num].T;
       //q = cells[num].Q[1];
       //tau = cells[num].PIJ[3];
       //Rho = (*gas).m*cells[num].n;
       //MetroHastingChEn(Nnew_DSMC[num], Rho, RT, q, tau, U_temp);

       //MomentsofSamples(Nnew_DSMC[num], U_temp, m, 3);
       //stdev = standard_deviation(U_temp, Nnew_DSMC[num], &mean);

       for(i=0; i<Nnew_DSMC[num]; i++){
         U1_temp[i] = dis_norm(gen);
         //U_temp[i]  = dis_norm(gen);
         U3_temp[i] = dis_norm(gen);
       }
       fix_mean_variance(U1_temp, Nnew_DSMC[num]);
       fix_mean_variance(U3_temp, Nnew_DSMC[num]);

       fix_mean_variance(U_temp, Nnew_DSMC[num]);
       //fix_mean(U_temp, Nnew_DSMC[num]);
       MomentsofSamples(Nnew_DSMC[num], U_temp, m, 3);
       printf("xx %lf and m[2] %lf\n", xx, m[2]);
       printf("lambda %lf %lf %lf\n", lambda[0], lambda[1], lambda[2]);
       printf("m %lf %lf %lf\n", m[0], m[1], m[2]);
       //stdev = standard_deviation(U_temp, Nnew_DSMC[num], &mean);


      for(i=0; i<Nnew_DSMC[num]; i++){
        std = sqrt((*gas).kb*cells[num].T/(*gas).m);
        //std = sqrt((*gas).kb*273.0/(*gas).m);
        x1new[j] = (*x1)[0];
        x2new[j] = cells[num].cell_center[1] +  dis_x(gen)*(*box).delta_dim[1]*0.5;
        U1new[j] = U1_temp[i]*std+cells[num].U_space[0];//dis_norm(gen)*std+cells[num].U_space[0];
        U2new[j] = U_temp[i]*std+cells[num].U_space[1];
        U3new[j] = U3_temp[i]*std+cells[num].U_space[2];//dis_norm(gen)*std+cells[num].U_space[2];
        colornew[j] = dsmc;
        x1_oldnew[j] = (*x1_old)[0];
        x2_oldnew[j] = x2new[j];
        Tnew[j] = cells[num].T;
        rhonew[j] = cells[num].n*(*gas).m;
        if(U2new[j]!=U2new[j]){
          printf("U2new[%d]=%e\n", j, U2new[j]);
        }
        j++;
      }
      //stdev = standard_deviation(&U2new[j-Nnew_DSMC[num]], Nnew_DSMC[num], &mean);
      cells[num].model=dsmc;
    }
    // sample sph particles in dsmc_to_sph cells
    else if( cells[num].model==dsmc_to_sph ){
      for(i=0; i<Nnew_SPH[num]; i++){
        x1new[j] = (*x1)[0];
        x2new[j] = cells[num].cell_center[1] -(*box).delta_dim[1]/2.0 + (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew_SPH[num]);
        U1new[j] = cells[num].U_space[0];
        U2new[j] = cells[num].U_space[1];
        U3new[j] = cells[num].U_space[2];
        colornew[j] = sph;
        x1_oldnew[j] = (*x1_old)[0];
        x2_oldnew[j] = x2new[j];

        Tnew[j] = cells[num].T;
        rhonew[j] = cells[num].n*(*gas).m;

        if(Tnew[j]!=Tnew[j]){
          printf("Tnew[%d]=%e, at num %d \n", j, Tnew[j], num);
        }

        j++;
      }
      cells[num].model=sph;
    }
  }
  if(N>0){
    free(U_temp); free(U1_temp); free(U3_temp);
  }
  //printf("j after changed cell = %d\n", j);
  // find hybrid cells & count total dummy SPH particles needed there
  for(num=1; num<(*box).N[1]-1; num++){
   if( cells[num].model==dsmc  && (cells[num-1].model==sph || cells[num+1].model==sph) ){
            cells[num].model=hybrid;
    }
  }
  /*
  //dummy sph particles in hybrid
  for(num=0; num<(*box).N[1]; num++){
    if( cells[num].model==hybrid ){
      printf("--- generating %d dummy sph particles in hybrid cell at num=%d \n",Nnew_Hybrid[num], num);
      for(i=0; i<Nnew_Hybrid[num]; i++){
        x1new[j] = (*x1)[0];
        x2new[j] = cells[num].cell_center[1] -(*box).delta_dim[1]/2.0 + (i+0.5)*(*box).delta_dim[1]/(1.0*Nnew_Hybrid[num]);
        U1new[j] = cells[num].U_space[0];
        U2new[j] = cells[num].U_space[1];
        U3new[j] = cells[num].U_space[2];
        colornew[j] = sph;
        x1_oldnew[j] = (*x1_old)[0];
        x2_oldnew[j] = x2new[j];

        Tnew[j] = cells[num].T;
        rhonew[j] = cells[num].n*(*gas).m;
        j++;
        if(Tnew[j]!=Tnew[j]){
          printf("Tnew[%d]=%e, at num %d \n", j, Tnew[j], num);
        }
      }
    }
  }
  */
  if(j!=Ntot){
    printf("j=%d Ntot=%d\n",j,Ntot);
  }
  /////  reallocating memory
  free(*x2);free(*U1);free(*U2);free(*U3);free(*x2_old);
  free(*color); free(*x1); free(*x1_old); free(*T); free(*rho);
  (*gas).N = Ntot;
  *x1 = x1new;
  *x2 = x2new;
  *T = Tnew;
  *rho = rhonew;
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

  free(Nnew_SPH); free(Nnew_DSMC);
}
