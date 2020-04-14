#ifndef CPP_CHEADER_H_INCLUDED
#define CPP_CHEADER_H_INCLUDED
//#include <omp.h>
//#include "omp.h"

//#include <mkl.h>
//#include <mkl_cblas.h>
//#include <mkl_blas.h>
//#include <mkl_lapack.h>
//#include <mkl_lapacke.h>

//# include "ziggurat.h"

#include <math.h>
#include <random>
#include <iostream>
#include <cmath>
#include <string>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
 #include <stdio.h>
#define NSD 3

#define closed_box 1
#define couette_flow 2
#define evaporation 3
#define vacuum 4
#define dcones 5
#define relaxation 6
#define flatnose 7
#define inverted 8
#define zoomed_inverted 9
#define shock 10
#define wall 11
#define injection 12

#define sph 1
#define sph_to_dsmc 2
#define dsmc_to_sph 3
#define dsmc 4
#define hybrid 5

#define GP 0
#define net 1

/*
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif
*/
#ifdef __cplusplus
extern "C" {
#endif

void dgesv(const int *, const int *, double *, const int *, int *, double *, const int *, int *);

void dsgesv(const int *, const int *, double *, const int *, int *, const double *, const int *, double *, const int *, double *, float *, int *, int *);

/*
extern MKL_INT PARDISO
	(void *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
	double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
	MKL_INT *, double *, double *, MKL_INT *);
*/
int solving_phi(int *ia, int *ja, double *a, int n, double *b, double *x, int ) ;

/* PARDISO prototype. */
/*
//void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardisoinit (void *pt, MKL_INT *mtype, MKL_INT *iparm );

void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);
*/

uint32_t cong_seeded ( uint32_t *jcong );
double cpu_time ( );
uint32_t kiss_seeded ( uint32_t *jcong, uint32_t *jsr, uint32_t *w, uint32_t *z );
uint32_t mwc_seeded ( uint32_t *w, uint32_t *z );
float r4_exp ( uint32_t *jsr, uint32_t ke[256], float fe[256], float we[256] );
void r4_exp_setup ( uint32_t ke[256], float fe[256], float we[256] );
float r4_nor ( uint32_t *jsr, uint32_t kn[128], float fn[128], float wn[128] );
void r4_nor_setup ( uint32_t kn[128], float fn[128], float wn[128] );
float r4_uni ( uint32_t *jsr );
uint32_t shr3_seeded ( uint32_t *jsr );
void timestamp ( );

//extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
//                double* b, int* ldb, int* info );

}

#define L11 10
#define L22 20
#define L33 40
#define dimX 1
#define dimY 3

struct NerualNetwork{
  double W1[dimX*L11];
  double W2[L11*L22];
  double W3[L22*L33];
  double Wo[L33*dimY];

  double b1[dimX*L11];
  double b2[dimX*L22];
  double b3[dimX*L33];
  double bo[dimX*dimY];

  double mean_input[dimX];
  double dev_input[dimX];
  double mean_output[dimY];
  double dev_output[dimY];
};

#define Ndata 2000

struct GPR{
  double kst[Ndata*dimY];
  double X[Ndata];
  double kxX[dimX*Ndata];
  double k_lsc;
  double k_var;
  double loglike_var;

  double mean_input[dimX];
  double variance_input[dimX];
  double mean_output[dimY];
  double variance_output[dimY];
};

struct GAS
{
  int gp_nn;
  double r_cut, h, mb;
  double N_abs;
  double Lc_ratio,Lh_ratio, Lc_center, Lh_center;
  double sigc, alphac,sigh,alphah;
  double count_SP_c, count_SP_h,  count_reem_c, count_reem_h, count_evap_c, count_evap_h;
  double nvc, nvh;
  double nc, nh;
  double nv1, nv2, nv3;
  double Tc, Th, Tv1, Tv2, Tv3;
  double Lc, Lh, Lv, Lv1, Lv3;
  int LcN, LhN, LvN, LvN1, LvN3;// (double) Lc = (int) LcN*eqdist
  int Nc1,Nc2,Nc3, Nh1,Nh2,Nh3;
  int Nc, Nh;
  int Nv1, Nv2, Nv3;
  int N13;
  double phi;
  double eqdist;
  int n_inflow;
  double ast;
  double bst;
  double sigma;
  double n;
  long int N;
  double w;
  double m;
  double lambda;
  double Kn;
  double T;
  double T0;
  double delta_t;
  double frac_mft;// what fraction of mean free time
  double kb;
  double U0[NSD];
  double U[NSD];
  double u_avg;
  double mu;
  double mu_corr;
  double M[NSD*2]; //change of moments on the boxes surfaces
  std::string model;
  std::string output_name;
  double p;
  double Fn, Fn_SPH;
  double factor;
  double phi0;
  double b;
  double tau_shear[9];
  double nualpha;
  double visp;
  double crref;
  double t_post;
  double avgeraging_time_till_now;
  double times;
  double epsilon;
  int long_range;
  int n_ratio;
  double k;
};

struct BOX
{
  double drop_dia;
  double x2_prob[4];
  double ghost,thermc, thermh;
  double JR[4], JL[4], Jref[4];
  double NR[4], NL[4];
  double L1[2], L2[2];
  int N1[2], N2[2];
  double s1, s2;
  double dcones_n;
  double dy0, dx0;
  double dy1, dx1;
  double alpha[2];
  double Lx[3];
  double Ly;
  double s;
  int Ny;
  int Nx[3];
  int N_grid[3];// for MD, the grid for putting particles in
  double Len[NSD];
  double delta_dim[NSD];
  double Area[NSD*2];
  double area, dummy_nom;
  double Volume;
  int N[NSD]; // number of cells in each direction
  int num_steps;
  int num_realizations;
  int movie;
  int every;
  int after;
  int reset_every;
  int init_step;
  int num_thr;
  int problem;
  double U_wall_1;
  double U_wall_2;
  double U_wall_3;
  double U_wall_4;
  double U_wall_5;
  double U_wall_6;
  double U_wall_7;
  double T_wall_1;
  double T_wall_2;
  double T_wall_3;
  double T_wall_4;
  double T_wall_5;
  double T_wall_6;
  double T_wall_7;
  double n_wall_4;
  double n_wall_5;
  double n_wall_6;
  double n_wall_7;
  int Ndot_w7;
  int Ndot_w6;
  int Ndot_w5;
  int Ndot_w4;
  int N_tot_genenrate;
  std::string file_name;
  int reset;
  int direction[3];
  int step;
  double out;
};

struct CELLS
{
  double p_m;
  double UU, UU0, Force;
  double UUR, UU0R;
  double divU;
  int model;
  int *id_split_part;
  int n_ratio;
  int id[4];
  double n_n[4];
  double x_n[4];
  double y_n[4];
  int in;// if inside then in = 1
  int ng[2];
  double V[4];// potential on the faces
  double dx, dy;
  double dn2;
  double dn3;
  double n_smooth;
  double weight;
  double cell_center[3];
  double dim[NSD*2];
  double volume;
  int *indices_inside;// array of the indices of the particles that are inside this cell
  int num_inside;// number of particles inside this cell
  double alpha[9];
  double sum_alpha[9];
  double beta[3];
  double U_space[3]; // average velocity of cell in space only
  double T;
  double n;
  double U[3];
  ///  moments needed for FP
  double PIJ[6];
  double M4[6];
  double M5[3];
  double Q[3];
  double M3[10];
  double DM2;
  double DM4;
  double c[6];
  double gamma[3];
  double gamma_st;
//// for DSMC
  double crm;
///////   for bulk viscosity effect
  double M_n[6];
  double sum_M_n[6];
  double sum_M_np[6];
  double T_f[6];
  double M1_f[6];
  double M2_f[6];
  double M3_f[6];
  double sum_T_f[4];
  double sum_M1_f[6];
  double sum_N[4];
  double sum_M2_f[6];
  double sum_M3_f[6];
  double sum_MM_f_0[6];
  double sum_MM_f_1[6];
  double sum_MM_f_2[6];
  double sum_MM_f_3[6];
  // averaging for rho*cc ?= (1+2nbY/5)*dU/dY
  double sum_MMM[10];// 111 122 133 211 222 233
  double sum_MMMM[6];
  double sum_MM[6];
  double sum_M[3];
  double sum_weight;
  double Mp[3];
  double phi[3];
  double omega_max;
  double Lambda;
  double G_sum;
  double F1,F2;
  double xx;
};


struct Collision
{
  int num, Nc, cand_num, pair[2];
double pi, Vr_max, Vc, M_cand, Y, r, Vr,   vr_pre_coll[3], phi, theta, alpha1, alpha2, n[3];
double  v_r_prime[3], v_c[3],  vr_post_coll[3];
double  d[3], norm_post_pre;
double cos_teta,sin_teta, q;
double crm;
};

int solver_pardiso(int *ia, int *ja, double *a, int n, double *b, double *x);

void random_generator(double *xi_1, double *xi_2, double *xi_3, double *xi_1f, double *xi_2f, double *xi_3f, struct CELLS *cells, struct BOX *box,  struct GAS *gas, int *index);

void Collision_ESMC_new4(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void Collision_ESMC_new3(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);


void Collision_ESMC_new(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);

void zoomed_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag, int **n_ratio,double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3);
int reset_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color);

void min_rij(int n, double *rx, double *ry, double *rz, int *nn, double *rr);

//extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
//                double* b, int* ldb, int* info );

void calculate_potential(double *V, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void linear_fit(double n1, double r1, double n2, double r2, double n3, double r3, double *a, double *b);

void MD_8th(double *U1, double *U2, double *U3, double *x1, double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *F1, double *F2, double *F3,
                  double *F1_old, double *F2_old, double *F3_old, double *F1_m1, double *F1_m2, double *F1_m3, double *F1_m4, double *F1_m5, double *F1_m6, double *F1_m7,
                                                                  double *F2_m1, double *F2_m2, double *F2_m3, double *F2_m4, double *F2_m5, double *F2_m6, double *F2_m7,
                                                                  double *F3_m1, double *F3_m2, double *F3_m3, double *F3_m4, double *F3_m5, double *F3_m6, double *F3_m7,
                                                                  double *x1_m1, double *x2_m1, double *x3_m1, double *x1_old, double *x2_old, double *x3_old, double *phi);

void measure_pressure(double *U2,double *x2, double *x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void shock_BC(double *U2, double *x2, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void compute_evap_coeff(double *U2,double *x2,double *x2_old,struct GAS *gas, struct BOX *box, int *color,int *n_ratio);
void compute_evap_coeff_evaporation(double *U2,double *x2,double *x2_old,struct GAS *gas, struct BOX *box, int *color,int *n_ratio,  struct CELLS *cells);
void ghost_inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3);

void read_GPR(struct GPR *gpr);

void compute_SP_2D(struct GAS *gas, struct BOX *box, struct CELLS *cells, double *x1, double *x2, double *U1, double *U2, int *index);
void write_post_processing(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index, int *color, int*n_ratio,double *x2_old);

void evaluation_GPR(struct GPR *gpr, double x, double *y);
void read_NN(struct NerualNetwork *NN);

void evaluation(struct NerualNetwork *NN, double x, double *y);

void DFP_extra_streaming(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt);
void SPH_simple_hybrid(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *color);

//void Collision_DSMC_VHS_hybrid(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC);
void Collision_DSMC_VHS_hybrid(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC, int *color);

void cell_update_hybrid(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double * T, double * rho, int *color);
void add_remove_hybrid(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho, struct NerualNetwork *NN, struct GPR *gpr);

void fix_mean(double *U_temp, int N);

int solver_pardiso();

void SPH_simple(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
double compute_DKL_MED_MW(double *lambda);

void MomentsofSamples(int nSamples, double *x, double *m, int Nm);

void reset_inverted_ghost(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color);
//void evap_BC(double *U2, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void evap_BC(double *U1,double *U2, double *U3, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells);

int max_int(int *a,int n);

void BC_wall(double *U1,double *U2,double *U3, double *x2, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x2_old);


void fix_mean_variance(double *U_temp, int N, double U, double T);

void solver( );
void addrem_BC_hybrid(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho);
void BC_specular_periodic_reflection(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void Collision_ESMC_dcones(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void MD_potential(double *dis, double *phi, struct GAS *gas);
double sum(double *data, int n);
void thermostat(double *U1,double *U2,double *U3, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void BC_inverted(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, double *x1_old, double *x2_old, double *x3_old);
void writting_position_to_file2D(double *x1, double *x2, int N, double *t, int len_t, int step);
void post_processing(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index, int *color, int*n_ratio,double *x2_old);
void  generate_new_particles(double *x1n,double *x2n,double *U1n,double *U2n,double *U3n,double *x1_oldn,double *x2_oldn, struct GAS *gas, struct BOX *box, double *dtn,  struct CELLS *cells);
double volume_conical_frustum(double R1, double R2, double h);
int new_particles(struct GAS *gas, struct BOX *box,  struct CELLS *cells);
//void poisson_solver(double *U2, double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void poisson_solver(double *U2, double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
int dcone_index(double x, double y, struct CELLS *cells,  struct BOX *box, int *j, int *k);
void write_relaxation( struct CELLS *cells,  struct BOX *box, int n);
double geometric_sum(int n, double s);
void store_old_positions(double *x1,double *x2, double *x3, double *x1_old, double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box);
double generateGaussianNoise(double mu, double sigma);

void inverted_addremove_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag, int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color);

void Velocity_FP_opt(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3);
void solver(double *a, double *b, int n, int nrhs );
void read_pos_vel(double *x1, double *x2, double *U1, double *U2, double *U3, int N);
void write_pos_vel(double *x1, double *x2, double *U1, double *U2, double *U3, int N);
void derivative_T_dy( struct BOX *box, struct CELLS *cells);
void write_averaging(struct CELLS *cells,  struct BOX *box, struct aveg *avg);
void derivative_U_dy( struct BOX *box, struct CELLS *cells);
void averaging (double *U1,double *U2,double *U3, struct CELLS *cells,  struct BOX *box, int step, struct GAS *gas, int *index, struct aveg *avg);
void c_gamma_calculation( double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cell, struct BOX *box);
void Velocity_FP_new(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index);
void calculate_pressure_tensor(struct GAS *gas, struct BOX *box, struct CELLS *cells);
void cell_poly_velocity(double *U1,double *U2,double *U3, struct CELLS *cells,  struct BOX *box, int step, struct GAS *gas, int *index);
double min(double a, double b);
void stream_axisym(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt);
void update_face_info(double *U1,double *U2,double *U3, double *x1,double *x2, double *x1_old,double *x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
double line_interesect_y(double x, double x_old, double y_old, double x_new, double y_new);
double line_interesect_x(double y, double x_old, double y_old, double x_new, double y_new);
int flatnose_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x1_old,double *x2_old, int* flag, struct GAS *gas, struct BOX *box, struct CELLS *cells,
double *U1n,double *U2n,double *U3n,double *x1n,double *x2n,double *x1_oldn,double *x2_oldn, int* flagn, double *dtn);
double sign(double x);
int flatnose_index(double x, double y, struct CELLS *cells,  struct BOX *box, int *j, int *k);
void dcones_vlasov_integral(double *U1, double *U2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void Collision_ESMC_new2(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void Collision_ESMC(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *ESMC, int *index);

void MD(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);
void MD_force(double *dis, double *F1,double *F2,double *F3, struct GAS *gas);

void MD_new(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index,
         double *F1, double *F2, double *F3);

void Position_FP_ideal_gas(double *U1,double *U2,double *U3, double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt);


void BC_for_CBA(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *flag);

void stream_particles(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);


//void apply_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *flag,double *xi_1, double *xi_2, double *xi_3, double *xi_1f, double *xi_2f, double *xi_3f, double *Mp1, double *Mp2, double *Mp3, double *T);
void apply_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, double *x1_old,double *x2_old, double *x3_old, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, int *flag);

void velocity_update(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, struct Collision *CBA, double *xi_1f, double *xi_2f, double *xi_3f, double *Mp1, double *Mp2, double *Mp3, double *F1, double *F2, double *F3, double *rho, double *T, int* color);

void cell_update(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double ** T, double * rho, int *color);

void BC_evaporation2(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, double *x1_old, double *x2_old, double *x3_old);
void scale_mu(double T_want, double T, double *mu_want, double mu, struct GAS *gas);
double increase_collision_rate(double n, struct GAS *gas);

void Collision_DSMC_HS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC);

void Thompson_3d(double *a,double *b,double *c,double *x, double *RHS,int n);

void Velocity_FP_opt_mod(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *xi2_1, double *xi2_2, double *xi2_3);

void Position_FP_new3(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double dt);

void vlasov_term(double *U2, double *V, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);

void Velocity_FP_linear(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *Mp1, double *Mp2, double *Mp3);

void Position_FP_new(double *U1,double *U2,double *U3, double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index,  double *xi_1, double *xi_2, double *xi_3);
void shear_stress_in_total( double *M1,double *M2,double *M3, struct BOX *box, struct CELLS *cells, int new_step, struct GAS *gas);
void average_velocity(double *U1,double *U2,double *U3, struct GAS *gas);
void cell_temperature_mu(double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int step, int after);
void mean_mu_cell(struct CELLS *cells,  struct BOX *box, struct GAS *gas, double *mean_mu_s, int new_step);

void curve_fit_U(struct CELLS *cells,  struct BOX *box, struct GAS *gas, double *slope, double *a);

void gas_update(double *U1,double *U2,double *U3, struct GAS *gas, int new_step);

int dcones_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x1_old,double *x2_old, int* flag, struct GAS *gas, struct BOX *box, struct CELLS *cells,
double *U1n,double *U2n,double *U3n,double *x1n,double *x2n, double *x1n_old,double *x2n_old, int* flagn, double *dtn);

void update_face_info_BC(double U1,double U2,double U3, double x1,double x2, double x1_old,double x2_old, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void cells_mean_velocity(double *U1,double *U2,double *U3,struct BOX *box, struct CELLS *cells, int new_step);

//void hybrid_surface_flux(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho,  struct GPR *gpr);

void BC_thermal_wall(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag,double *x1_old, double *x2_old);

  void shear_stress_in_each_cell( double *M1,double *M2,double *M3, struct BOX *box, struct CELLS *cells, int new_step, struct GAS *gas);
void Gauss_Jordan(double **A,double *RHS, int n, double *x);
void alpha_calculation( double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells, struct BOX *box);

void curve_fit(double *x, double *y, int n, int *indices, double *a, double *b);

void add_ghost_particles(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho,  struct GPR *gpr, struct NerualNetwork *NN);

void hybrid_surface_flux(double **U1,double **U2,double **U3,double **x1,double **x1_old,double **x2,double **x2_old,struct GAS *gas, struct BOX *box, struct CELLS *cells, int **index, int **flag,  int **n_ratio, double **xi_1, double **xi_2, double **xi_3, double **xi_1f, double **xi_2f, double **xi_3f, double **Mp1, double **Mp2, double **Mp3, int **color,  double ** T, double ** rho,  struct GPR *gpr, struct NerualNetwork *NN);





void CSR_Create(int M,int N, int nz, int *I, int *J, double *val,int *AI_CSR);
void A_in_Ax1(double *A,int *I, int *J, int nonzero, double N0, double N1, double N2, struct BOX *box);
void A_in_Ax2(double *A,int *I, int *J, int nonzero, double N0, double N1, double N2, struct BOX *box,  struct GAS *gas, struct CELLS *cells);
void A_in_Ax3(double *A,int *I, int *J, int nonzero, double N0, double N1, double N2, struct BOX *box,  struct GAS *gas, struct CELLS *cells);
void A_in_Ax4(double *A,int *I, int *J, int nonzero, double N0, double N1, double N2, struct BOX *box,  struct GAS *gas, struct CELLS *cells);
void rhs_in_Ax_eq_rhs1(double *phi, double *rhs, struct GAS *gas,  struct CELLS *cells, struct BOX *box);
void rhs_in_Ax_eq_rhs2(double *rhs, struct GAS *gas,  struct CELLS *cells, struct BOX *box);
void rhs_in_Ax_eq_rhs3(double *rhs, struct GAS *gas,  struct CELLS *cells, struct BOX *box);
void rhs_in_Ax_eq_rhs4(double *rhs, struct GAS *gas,  struct CELLS *cells, struct BOX *box);

void rotate(double *x0,double *y0, double theta, double *x, double *y);

void post_vtk(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index);
void initialization(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void BC_specular_reflection(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box);
void find_error_psi(double *x1,double *x2,double *x3, struct GAS *gas, double *psi);
void find_inside_cells(double *x1,double *x2,double *x3,double gas_N, struct BOX *box, struct CELLS *cells, int *index);
void avg_psi_in_cell(double *psi, struct BOX *box, struct CELLS *cells);
void particle_to_cell_phi(double *x1,double *x2,double *x3, struct BOX *box, struct CELLS *cells, double *psi, double *phi);

void dphi_error(struct BOX *box, struct CELLS *cells, double *phi, double *dphi_dx, double *dphi_dy, double *dphi_dz);
void dphi_cell_to_particle(struct BOX *box, struct CELLS *cells, double *x1,double *x2,double *x3, double *dphi_dx, double *dphi_dy, double *dphi_dz , double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p, double *psi, struct GAS *gas);

void BC_evaporation(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, double *x3_old);

void BC_thermal_wall_3D(double *U1,double *U2,double *U3,double *x1,double *x2, double *x3, struct GAS *gas, struct BOX *box, int *index, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, double *x3_old);
void intersect_3D_x(double x, double *y, double *z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new);
void intersect_3D_y(double *x, double y, double *z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new);
void intersect_3D_z(double *x, double *y, double z, double x_old, double y_old, double z_old, double x_new, double y_new, double z_new);

void BC_periodic(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void MD(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index,
         double *F1, double *F2, double *F3, double *F1_old, double *F2_old, double *F3_old);

void find_error_dphi_p_directly3(double *x1,double *x2,double *x3, struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p, double *psi);
void Collision_ESMC_1D(double *U1,double *U2,double *U3,double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *ESMC, int *index);

//void Velocity_FP(double *M1,double *M2,double *M3, struct GAS *gas);
void Velocity_FP(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index);
void Velocity_FP_with_dia(double *M1,double *M2,double *M3, struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p);
void dph_p_correction(struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p);


void find_error_dphi_p_directly(double *x1,double *x2,double *x3, struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p);

void Position_FP_new2(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);

void vlasov_term(double *U2, double *V, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *x1, double *x2, double *x3);

void Position_FP(double *U1,double *U2,double *U3, double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells, double *phi);
void Position_FP_idealgas(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void Position_FP2(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells);
void Position_FP3(double *U1,double *U2,double *U3, double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void Velocity_FP_with_dia_directly(double *M1,double *M2,double *M3, struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p);

void Collision_CBA_HS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells , struct Collision *CBA);
void Collision_CBA_VHS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *CBA);
void Collision_DSMC_VHS(double *U1,double *U2,double *U3,double *x1,double *x2, struct GAS *gas, struct BOX *box, struct CELLS *cells, struct Collision *DSMC);

void shuffle(double *x1,double *x2,double *x3, struct GAS *gas);

void find_error_psi_string_model(double *x1,double *x2,double *x3, struct GAS *gas, double *dx, double *dy, double *dz);

void making_inputs(struct GAS *gas,struct BOX *Box);
void make_cells(struct BOX *Box, struct CELLS *cells, struct GAS *gas);
void writting_position_to_file(double *x1, double *x2, double*x3, int N, double *t, int len_t, int step);
void writting_data_to_file(int num_steps, double delta_t, double *p, double* p_kin, double *T, struct GAS *gas, double* E_kin, int realization, double* p_mod);
void writting_factor_to_file(double T, double n, double factor);
void generating_Python_paraview_script(int num_step, struct BOX *Box);
void readSettingsFile(struct GAS *gas,struct BOX *box);

void find_error_dphi_p_directly2(double *x1,double *x2,double *x3,double *U1,double *U2,double *U3, struct GAS *gas, double *dphi_dx_p, double *dphi_dy_p, double *dphi_dz_p);

void vlasov_taking_the_integral(double *U2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index, double *x1, double *x2, double *x3);

void update_index_number_inside(double *x1, double *x2, double *U1, double *U2, double *U3, double *Mp1, double *Mp2, double *Mp3, int *index, struct CELLS *cells,  struct BOX *box, struct GAS *gas);

void Velocity_FP_parallel(double *M1,double *M2,double *M3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *xi_1, double *xi_2, double *xi_3, double *Mp1, double *Mp2, double *Mp3);

void Velocity_FP_string_model(double *M1,double *M2,double *M3, struct GAS *gas, double *dx, double *dy, double *dz);
  void find_error_psi_elastic_collision(double *x1,double *x2,double *x3,double *U1,double *U2,double *U3, struct GAS *gas, double *dv_x, double *dv_y, double *dv_z);
  void Velocity_FP_elastic_collision(double *M1,double *M2,double *M3, struct GAS *gas, double *dv_x, double *dv_y, double *dv_z);

  void SPH_velocity_temperature_update(double *x2, double *U2, double *rho, double *T, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index);

void  Velocity_FP_sphere(double *M1,double *M2,double *M3, struct GAS *gas, double *dv_x, double *dv_y, double *dv_z);
void find_error_psi_sphere_number(double *x1,double *x2,double *x3,double *U1,double *U2,double *U3, struct GAS *gas, double *dv_x, double *dv_y, double *dv_z);

void find_error_psi_elastic_collision2(double *x1,double *x2,double *x3,double *U1,double *U2,double *U3, struct GAS *gas, double *dv_x, double *dv_y, double *dv_z);
// void Position_FP_new2(double *U1,double *U2,double *U3, double *x1,double *x2,double *x3, struct GAS *gas, struct BOX *box, struct CELLS *cells);

void cell_update_post(double *x1,double *x2, double *U1,double *U2,double *U3, struct GAS *gas, struct CELLS *cells,  struct BOX *box, int *index, double *Mp1,double *Mp2,double *Mp3);

double standard_deviation(double data[], int n, double *mean);
double mean(double *data, int n);
double max(double a, double b);
void CG_Solver(
	 int M,
	 int N,
	 int nz,
	 int *AI_CSR,
	 int *RJ,
	 double *Rval,
	 int *AJ_CSC,
	 int *CI,
	 double *Cval,
	 double *x0,
	 double *x_m,
	 double *x,
	 double *b,
	 int guessed_m,
	 char SG,
	 int M_kind,
	 double tol,
	 double *Err);
void Residual_Ax_b(int M,int N, int nz, int *AI_CSR, int *RJ, double *Rval, int *AJ_CSC, int *CI, double *Cval, double *x, double *b,double *r, char SG);
void Ax(int M,int N, int nz, int *AI_CSR, int *RJ, double *Rval, int *AJ_CSC, int *CI, double *Cval, double *x, double *y, char SG);
double scalar_prod(int n,double *x,double *y);
double norm_2(int n,double *x);
double dot_prod(double *val,double *x,int i1,int i2, int *J);
double *Left_Preconditioned(int M, int N, int nz, int *AI_CSR, int *RJ, double *Rval, double *b, int M_kind);
void CSC_Create(int M,int N, int nz, int *I, int *J, double *val,int *AJ_CSC, int *CI, int *CJ, double *Cval);
void get_krylov
	(double *v,
	 double *h,
	 int M,
	 int N,
	 int nz,
	 int *AI_CSR,
	 int *RJ,
	 double *Rval,
	 int *AJ_CSC,
	 int *CI,
	 double *Cval,
	 int j,
	 char SG,
	 int *happy_break,
	 int M_kind);
void GMRES
	(int M,
	 int N,
	 int nz,
	 int *AI_CSR,
	 int *RJ,
	 double *Rval,
	 int *AJ_CSC,
	 int *CI,
	 double *Cval,
	 double *x0,
	 double *x_m,
	 double *b,
	 int guessed_m,
	 char SG,
	 int M_kind,
	 double tol,
	 int *act_m);

void Restarted_GMRES
	(int M,
	 int N,
	 int nz,
	 int *AI_CSR,
	 int *RJ,
	 double *Rval,
	 int *AJ_CSC,
	 int *CI,
	 double *Cval,
	 double *x0,
	 double *x,
	 double *b,
	 int guessed_m,
	 char SG,
	 int M_kind,
	 double tol,
	 int *act_m,
	 int *iter);
#endif /* CHEADER_H_INCLUDED */
