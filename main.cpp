#include "cpp_headers.h"

int main(int argc, const char **argv)
{
  double start,end;
  struct GAS gas;
  struct BOX box;
  struct CELLS *cells;
  struct Collision CBA;
  struct NerualNetwork NN;
  struct GPR gpr;

  // check if the input script is given as an argument or not
  FILE *fp;
  int i,j,k, ii, iii;
  if (argc > 1){
    box.file_name = argv[1];
  } else {
    std::cout << " problem reading comand line \n";
    return 1;
  };
  ////					 reading input script
  readSettingsFile(&gas, &box);
  ////          creating cells
  cells = (struct CELLS *) malloc( box.N[0]*box.N[1]*box.N[2] * sizeof(struct CELLS) );
  make_cells(&box, cells, &gas);
  for (int j=0; j < box.N[0]; j++)
    for (int k=0; k < box.N[1]; k++)
        for (int l=0; l < box.N[2]; l++){
	  int num = l*box.N[0]*box.N[1] + k*box.N[0] +j;
	  cells[num].indices_inside = (int *) malloc( 1 * sizeof(int) );
	  for(int i=0; i<1; i++){
	    cells[num].indices_inside[i] = -1;}
	}
////////			allocating memory for variables
// U1, U2, U3 are the velocities
// and x1, x2 are the position of each particle
  double *U1, *U2, *U3, *x1, *x2,*x3, *x3_old;
// index stores the cell id for each particle
  int *index, *flag, done;
  int *color;
  U1 = (double *) malloc( gas.N * sizeof(double) );
  U2 = (double *) malloc( gas.N * sizeof(double) );
  U3 = (double *) malloc( gas.N * sizeof(double) );
  x1 = (double *) malloc( gas.N * sizeof(double) );
  x2 = (double *) malloc( gas.N * sizeof(double) );
  index = (int *) malloc( gas.N * sizeof(int) );
  flag = (int *) malloc( gas.N * sizeof(int) );
  color = (int *) malloc( gas.N * sizeof(int) );
  x3 = (double *) malloc( gas.N * sizeof(double) );
	x3_old = (double *) malloc( gas.N * sizeof(double) );
  double * x2_old = (double *) malloc( gas.N * sizeof(double) );
  double * x1_old = (double *) malloc( gas.N * sizeof(double) );
   double *T, *rho;
   double *F1, *F2, *F3, *F1_old, *F2_old, *F3_old;
double *xi_1f, *xi_2f, *xi_3f, *xi_1, *xi_2, *xi_3;
double *V, *V0, *Ek, *Ek0;
double *Mp1, *Mp2, *Mp3;
int *n_ratio;
/////////////////////////////////////////////////////////////////////

// initialize variables
 initialization(U1, U2, U3, x1, x2, x3, &gas, &box, cells);
 for(i=0; i<gas.N; i++)
       x1_old[i] = x1[i];
 for(i=0; i<gas.N; i++)
       x2_old[i] = x2[i];
 printf("nkT of the system = %e\n",gas.n*gas.kb*gas.T);
 for( i=0; i<gas.N; i++){
 flag[i] = 0;
 color[i] = 0;
 }

  double std = sqrt(gas.kb*gas.T/gas.m);
  double mean;
  std::cout<<"std of U suppoesed to be = "<<std<<"\n";
  std::cout<<"std of U1 = "<<standard_deviation(U1,gas.N,&mean)<<"\n";

  printf("T0 = %e\n", gas.m*pow(standard_deviation(U1, gas.N, &mean),2)/gas.kb);
  printf("T0 = %e\n", gas.m*pow(standard_deviation(U2, gas.N, &mean),2)/gas.kb);
  printf("T0 = %e\n", gas.m*pow(standard_deviation(U3, gas.N, &mean),2)/gas.kb);
  printf("gas.T0 = %e\n", gas.T0);
  printf("gas.T = %e\n", gas.T);
  printf("gas.crref = %e\n", gas.crref);

  printf("\n\n -PostProc:\nevery %d\nafter %d\n%d steps\n\n",box.every, box.after, box.num_steps);

  double thermal;
  thermal =  sqrt(gas.kb*gas.T/gas.m);
  double t, total_time, st0;
  int step;
  int N_cell;
  int num;

  step = 0;
  box.step = step;
  t = 0.0;
  gas.times = 1.0;
  gas.N_abs = 0.0;

  // update cell info
  cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);
	N_cell = box.N[0]*box.N[1]*box.N[2];
  // set <M> and <MM> to zero
  double pi = acos(-1);
  for( num=0; num<N_cell; num++){
  thermal =  sqrt(gas.kb*cells[num].T/gas.m);
 	cells[num].crm=10.0*( (thermal)) ;//*pow(gas.crref/thermal,2.0*gas.visp-1.0) );
  cells[num].omega_max = 4.0*pi*gas.sigma*gas.sigma*cells[num].crm*cells[num].n*gas.delta_t;
	cells[num].U_space[0]=0.0;
	cells[num].U_space[1]=0.0;
	cells[num].U_space[2]=0.0;
	for(ii=0; ii<9; ii++)
	  cells[num].alpha[ii] = 0.0;
	for(ii=0; ii<3; ii++)
	  cells[num].beta[ii] = 0.0;
 	for(ii=0; ii<6; ii++)
	  cells[num].sum_MM[ii] = 0.0;
 	for(ii=0; ii<3; ii++)
	  cells[num].sum_M[ii] = 0.0;
        cells[num].sum_weight = 0.0;
  	if(gas.model == "ESMC")
        	cells[num].omega_max =cells[num].omega_max*increase_collision_rate(gas.n, &gas);
 }
  total_time = box.num_steps*gas.delta_t;

std::cout << "====== Start of Simulation =============================================================" << std::endl;
step = 1;
box.step = step;
while(t<total_time){
	box.step = step;
	if (step % box.every == 0)
		printf("%%%%%%%%%%  Step= %d, Np=%ld\n",step, gas.N);

  cell_update(x1, x2, U1, U2, U3, &gas, cells, &box, index, &T, rho, color);
  post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);
  write_post_processing(step, cells,  &box, &done, &gas, U1, U2, U3, x1, x2, flag, index, color,  n_ratio, x2_old);

	velocity_update(U1, U2, U3, x1, x2, x3, &gas, &box, cells, index, &CBA, xi_1f, xi_2f, xi_3f, Mp1, Mp2, Mp3, F1, F2, F3, rho, T, color);
  store_old_positions(x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box);
	stream_particles(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index);
	apply_BC(U1, U2, U3, x1, x2, x3, x1_old, x2_old, x3_old, &gas, &box, cells, index, flag);

	t = t + gas.delta_t;
  step++;
}

// release the memory
free(x1);
free(x2);
free(U1);
free(U2);
free(U3);
free(flag);
free(index);
free(color);
for(int num=0; num < box.N[0]*box.N[1]*box.N[2]; num++)
	free(cells[num].indices_inside);
free(cells);
return 0;
}
