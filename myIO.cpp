
#include "cpp_headers.h"
#include <iostream>
#include <cstring>
#include <string>



void read_GPR(struct GPR *gpr){
  int i,j;
  char comma;
  int Ni, Nj;
  double *x;
  FILE *fp;
  // read W1
  //fp = fopen("models/0.1/kst_2000.dat", "r");
  fp = fopen("models/2.0/kst_2000.dat", "r");

  Ni = Ndata; Nj = dimY;
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*gpr).kst[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("kst has been read!\n" );



  // read b1
  //fp = fopen("models/0.1/X_2000.dat", "r");
  fp = fopen("models/2.0/X_2000.dat", "r");

  Ni = Ndata; Nj = 1;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*gpr).X[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("X has been read!\n" );


  double dummy = 0.0;
  //fp = fopen("models/0.1/k_lsc_2000.dat", "r");
  fp = fopen("models/2.0/k_lsc_2000.dat", "r");
  fscanf(fp, "%lf", &dummy);
  (*gpr).k_lsc = dummy;
  fclose(fp);

  //fp = fopen("models/0.1/k_var_2000.dat", "r");
  fp = fopen("models/2.0/k_var_2000.dat", "r");
  fscanf(fp, "%lf", &dummy);
  (*gpr).k_var = dummy;
  fclose(fp);

  //fp = fopen("models/0.1/loglike_var_2000.dat", "r");
  fp = fopen("models/2.0/loglike_var_2000.dat", "r");
  fscanf(fp, "%lf", &dummy);
  (*gpr).loglike_var = dummy;
  fclose(fp);

  // read mean input
  //fp = fopen("models/0.1/GPR_mean_input_0.1_new.dat", "r");
  fp = fopen("models/2.0/GPR_mean_input_2.0_new.dat", "r");
  Ni = dimX;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*gpr).mean_input[i] = x[i];
  }
  printf("GPR mean input has been read!\n" );

  // read mean input
  //fp = fopen("models/0.1/GPR_mean_output_0.1_new.dat", "r");
  fp = fopen("models/2.0/GPR_mean_output_2.0_new.dat", "r");
  Ni = dimY;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*gpr).mean_output[i] = x[i];
  }
  printf("GPR mean output has been read!\n" );



  // read dev input
  //fp = fopen("models/0.1/GPR_variance_input_0.1_new.dat", "r");
  fp = fopen("models/2.0/GPR_variance_input_2.0_new.dat", "r");
  Ni = dimX;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*gpr).variance_input[i] = x[i];
  }
  printf("GPR var input has been read!\n" );



  // read dev input
  //fp = fopen("models/0.1/GPR_variance_output_0.1_new.dat", "r");
  fp = fopen("models/2.0/GPR_variance_output_2.0_new.dat", "r");
  Ni = dimY;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*gpr).variance_output[i] = x[i];
  }
  printf("GPR var output has been read!\n" );


  free(x);
}



void read_NN(struct NerualNetwork *NN){
  int i,j;
  char comma;
  int Ni, Nj;
  double *x;
  FILE *fp;
  // read W1
  //fp = fopen("models/functions_3_0.1_new/W1.txt", "r");
  //fp = fopen("models/functions_3/W1.txt", "r");
  fp = fopen("models/functions_3_2.0_new/W1.txt", "r");

  Ni = dimX; Nj = L11;
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).W1[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("W1 has been read!\n" );




  // read W2
  //fp = fopen("models/functions_3_0.1_new/W2.txt", "r");
  //fp = fopen("models/functions_3/W2.txt", "r");
  fp = fopen("models/functions_3_2.0_new/W2.txt", "r");
  Ni = L11; Nj = L22;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).W2[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("W2 has been read!\n" );


  // read W3
  //fp = fopen("models/functions_3_0.1_new/W3.txt", "r");
  //fp = fopen("models/functions_3/W3.txt", "r");
  fp = fopen("models/functions_3_2.0_new/W3.txt", "r");


  Ni = L22; Nj = L33;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).W3[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("W3 has been read!\n" );


  // read Wo
  //fp = fopen("models/functions_3_0.1_new/Wo.txt", "r");
  //fp = fopen("models/functions_3/Wo.txt", "r");
  fp = fopen("models/functions_3_2.0_new/Wo.txt", "r");

  Ni = L33; Nj = dimY;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).Wo[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("Wo has been read!\n" );


  // read b1
  //fp = fopen("models/functions_3_0.1_new/b1.txt", "r");
  //fp = fopen("models/functions_3/b1.txt", "r");
  fp = fopen("models/functions_3_2.0_new/b1.txt", "r");

  Ni = dimX; Nj = L11;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).b1[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("b1 has been read!\n" );



  // read b2
  //fp = fopen("models/functions_3_0.1_new/b2.txt", "r");
  //fp = fopen("models/functions_3/b2.txt", "r");
  fp = fopen("models/functions_3_2.0_new/b2.txt", "r");

  Ni = dimX; Nj = L22;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).b2[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("b2 has been read!\n" );


  // read b2
  //fp = fopen("models/functions_3_0.1_new/b3.txt", "r");
  //fp = fopen("models/functions_3/b3.txt", "r");
  fp = fopen("models/functions_3_2.0_new/b3.txt", "r");


  Ni = dimX; Nj = L33;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).b3[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("b3 has been read!\n" );


  // read b2
  //fp = fopen("models/functions_3_0.1_new/bo.txt", "r");
  //fp = fopen("models/functions_3/bo.txt", "r");
  fp = fopen("models/functions_3_2.0_new/bo.txt", "r");

  Ni = dimX; Nj = dimY;
  free(x);
  x = (double *) malloc( Ni*Nj * sizeof(double) );
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      if(j<Nj-1)
        fscanf(fp, "%lf%c", &x[i*Nj+j],&comma);
      else
        fscanf(fp, "%lf", &x[i*Nj+j]);
    }
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
    for(j=0; j<Nj; j++){
      (*NN).bo[i*Nj+j] = x[i*Nj+j];
    }
  }
  printf("bo has been read!\n" );


  // read mean input
  //fp = fopen("models/functions_3_0.1_new/mean_input.dat", "r");
  //fp = fopen("models/functions_3/mean_input.dat", "r");
  fp = fopen("models/functions_3_2.0_new/mean_input.dat", "r");


  Ni = dimX;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*NN).mean_input[i] = x[i];
  }
  printf("mean input has been read!\n" );

  // read mean input
  //fp = fopen("models/functions_3_0.1_new/mean_output.dat", "r");
  //fp = fopen("models/functions_3/mean_output.dat", "r");
  fp = fopen("models/functions_3_2.0_new/mean_output.dat", "r");


  Ni = dimY;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*NN).mean_output[i] = x[i];
  }
  printf("mean output has been read!\n" );



  // read dev input
  //fp = fopen("models/functions_3_0.1_new/deviation_input.dat", "r");
  //fp = fopen("models/functions_3/deviation_input.dat", "r");
  fp = fopen("models/functions_3_2.0_new/deviation_input.dat", "r");


  Ni = dimX;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*NN).dev_input[i] = x[i];
  }
  printf("dev input has been read!\n" );



  // read dev input
  //fp = fopen("models/functions_3_0.1_new/deviation_output.dat", "r");
  //fp = fopen("models/functions_3/deviation_output.dat", "r");
  fp = fopen("models/functions_3_2.0_new/deviation_output.dat", "r");


  Ni = dimY;
  free(x);
  x = (double *) malloc( Ni * sizeof(double) );
  for(i=0;i<Ni;i++){
      fscanf(fp, "%lf", &x[i]);
  }
  fclose(fp);
  for(i=0;i<Ni;i++){
      (*NN).dev_output[i] = x[i];
  }
  printf("dev output has been read!\n" );
}
/*
void read_pos_vel(double *x1, double *x2, double *U1, double *U2, double *U3, int N){
      int i;
      FILE *fp;
      fp = fopen("pos_vel.txt", "r");
      for( i = 0; i < N ;i++){

	fscanf(fp, "%f", &x1[i]);
	fscanf(fp, "%lf", &x2[i]);
	fscanf(fp, "%lf", &U1[i]);
	fscanf(fp, "%lf", &U2[i]);
	fscanf(fp, "%lf", &U3[i]);
      }
      fclose(fp);
}
*/
void write_pos_vel(double *x1, double *x2, double *U1, double *U2, double *U3, int N){
      int i;
      FILE *fp;
      fp = fopen("pos_vel.txt", "w");
      for( i = 0; i < N ;i++){
	fprintf(fp, "%20e%20e%20e%20e%20e\n",
			x1[i],
			x2[i],
			U1[i],
			U2[i],
			U3[i]
 	      );
      }
      fclose(fp);
}

void write_relaxation( struct CELLS *cells,  struct BOX *box, int n){
      int i;
      FILE *fp;
      fp = fopen("relaxation.txt", "a");

       fprintf(fp, "%20d%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e\n",
		    n,
		    cells[0].weight,
		    cells[0].U_space[0],
		    cells[0].U_space[1],
		    cells[0].U_space[2],
		    cells[0].PIJ[0],
		    cells[0].PIJ[1],
		    cells[0].PIJ[2],
		    cells[0].PIJ[3],
		    cells[0].PIJ[4],
		    cells[0].PIJ[5]);

      fclose(fp);
}

void post_vtk(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index){
const char * a = (*gas).output_name.c_str();
const char * b;
if( (*gas).long_range == 1)
	b = "direct";
else if( (*gas).long_range == 2 )
	b = "denexp";
else if( (*gas).long_range == 3 )
	b = "screened-Poisson";
else
	b = "noatt";
  FILE *fp;
  char name_output[50];
  snprintf(name_output, 50, "%s_%s_%d.vtk",a,b,step);
  fp = fopen(name_output, "w");
  fprintf(fp, "# vtk DataFile Version 3.1\nthis is an example file created for VisuSimple\nASCII\nDATASET UNSTRUCTURED_GRID\n\n");
  int Nx, Ny, Ncells, Npoints;
  int j, k, vt, num, id;
  if((*box).problem == flatnose){
	Nx = (*box).N1[0]+(*box).N1[1];
	Ny = (*box).N2[0]+(*box).N2[1];
  }
  Ncells = 0;
  Npoints = (Ny+1)*(Nx+1)-(*box).N1[1]*(*box).N2[0];
  fprintf(fp, "POINTS %d DOUBLE\n",Npoints);
  id = 0;
  for(k=0; k<Ny; k++){
	for(j=0; j<Nx; j++){
		num = k*Nx+j;
		if(cells[num].in==1){
		  Ncells ++;
  		  if(j==0 && k==0){
			cells[num].id[0] = id;
			cells[num].id[1] = id+1;
			cells[num].id[2] = id+2;
			cells[num].id[3] = id+3;
			id = id+4;
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]-cells[num].dx/2.0, cells[num].dim[0]);
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[1]);
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]-cells[num].dx/2.0, cells[num].dim[2]);
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[3]);
		  }
		  else if(j==0 && k!=0){
			cells[num].id[0] = cells[(k-1)*Nx].id[2];
			cells[num].id[1] = cells[(k-1)*Nx].id[3];
			cells[num].id[2] = id;
			cells[num].id[3] = id+1;
			id = id+2;
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]-cells[num].dx/2.0, cells[num].dim[2]);
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[3]);
		  }
		  else if(j!=0 && k==0){
			cells[num].id[0] = cells[k*Nx+j-1].id[1];
			cells[num].id[1] = id;
			cells[num].id[2] = cells[k*Nx+j-1].id[3];
			cells[num].id[3] = id+1;
			id = id+2;
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[1]);
			fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[3]);
		  }
		  else{
			if(k==(*box).N2[0] && j>=(*box).N1[0]){
				cells[num].id[0] = cells[k*Nx+j-1].id[1];
				cells[num].id[1] = id;
				cells[num].id[2] = cells[k*Nx+j-1].id[3];
				cells[num].id[3] = id+1;
				id = id+2;
				fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[1]);
				fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[3]);
			}
			else{
				cells[num].id[0] = cells[k*Nx+j-1].id[1];
				cells[num].id[1] = cells[(k-1)*Nx+j].id[3];
				cells[num].id[2] = cells[k*Nx+j-1].id[3];
				cells[num].id[3] = id;
				id = id+1;
				fprintf(fp, "%lf  %lf  0.0\n",cells[num].cell_center[0]+cells[num].dx/2.0, cells[num].dim[3]);
			}
		  }
		}
	}
  }
  fprintf(fp, "\n\nCELLS %d %d", Ncells, Ncells*5);
  for(num=0; num<Nx*Ny; num++){
	if(cells[num].in==1){
		fprintf(fp, "\n4 ");
		fprintf(fp, "%d ", cells[num].id[0]);
		fprintf(fp, "%d ", cells[num].id[1]);
		fprintf(fp, "%d ", cells[num].id[3]);
		fprintf(fp, "%d ", cells[num].id[2]);
	}
  }
  fprintf(fp, "\n\nCELL_TYPES %d\n", Ncells);
  for(num=0; num<Nx*Ny; num++){
	if(cells[num].in==1)
		fprintf(fp, "9 ");
  }
  fprintf(fp, "\n\nCELL_DATA %d\nSCALARS Cell_Temperature DOUBLE\nLOOKUP_TABLE default", Ncells);
  for(num=0; num<Nx*Ny; num++){
	if(cells[num].in==1){
		fprintf(fp, "\n%lf ", cells[num].T);
	}
  }
  fclose(fp);
}


void write_post_processing(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index, int *color, int*n_ratio,double *x2_old){
  int i,j,jj,k,iii,num;
  FILE *fp;


  if ( ( (step > (*gas).times*(*box).after ) && step % (*box).every == 0) || step==(*box).init_step ){

 const char * a = (*gas).output_name.c_str();
 const char * b;
 if( (*gas).long_range == 1)
 	b = "direct";
 else if( (*gas).long_range == 2 )
 	b = "denexp";
 else if( (*gas).long_range == 3 )
 	b = "screened-Poisson";
 else
 	b = "noatt";
 //FILE *fp;
   char name_output[100];
   snprintf(name_output, 100, "PostProc/%s_%s_%d.dat",a,b,step);
   fp = fopen(name_output, "w");

     //  fp = fopen(a, "w");
       fprintf(fp, "%10s%10s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s%20s\n",
       "#variables=",
       "1-x1",
       "x2",
       "dx",
       "dy",
       "dz",
       "tavg",
       "dt",
       "U_space[0]",
       "U_space[1]",
       "10-U_space[2]",
       "T",
       "n",
       ///////  cell
       "sum_weight",
       "sum_M[0]",
       "sum_M[1]",
       "sum_M[2]",
       "sum_MM[0]",
       "sum_MM[1]",
       "sum_MM[2]",
       "20-sum_MM[3]",
       "sum_MM[4]",
       "sum_MM[5]",
       "sum_MMM[0]",
       "sum_MMM[1]",
       "sum_MMM[2]",
       "sum_MMM[3]",
       "sum_MMM[4]",
       "sum_MMM[5]",
       "sum_MMM[6]",
       "30-sum_MMM[7]",
       "sum_MMM[8]",
       "sum_MMM[9]",
       /// face
       "sum_T_f[0]",
       "sum_T_f[1]",
       "sum_T_f[2]",
       "sum_T_f[3]",
       "sum_M_n[0]",
       "sum_M_n[1]",
       "sum_M_n[2]",
       "40-sum_M_n[3]",
       "sum_M1_f[0]",
       "sum_M1_f[1]",
       "sum_M1_f[2]",
       "sum_M1_f[3]",
       "sum_M1_f[4]",
       "sum_M1_f[5]",
       "sum_M2_f[0]",
       "sum_M2_f[1]",
       "sum_M2_f[2]",
       "50-sum_M2_f[3]",
       "sum_M2_f[4]",
       "sum_M2_f[5]",
       "sum_M3_f[0]",
       "sum_M3_f[1]",
       "sum_M3_f[2]",
       "sum_M3_f[3]",
       "sum_N[2]",
       "sum_N[3]",
       "sum_M_np[0]",
       "60-sum_M_np[1]",
       "sum_M_np[2]",
       "sum_M_np[3]",
       "sum_M_np[4]",
       "sum_M_np[5]",
       "sum_MM_f_1[0]",
       "sum_MM_f_1[1]",
       "sum_MM_f_1[2]",
       "sum_MM_f_1[3]",
       "sum_MM_f_1[4]",
       "70-sum_MM_f_1[5]",
       "sum_MM_f_2[0]",
       "sum_MM_f_2[1]",
       "sum_MM_f_2[2]",
       "sum_MM_f_2[3]",
       "sum_MM_f_2[4]",
       "sum_MM_f_2[5]",
       "sum_MM_f_3[0]",
       "sum_MM_f_3[1]",
       "sum_MM_f_3[2]",
       "80-sum_MM_f_3[3]",
       "sum_MM_f_3[4]",
       "sum_MM_f_3[5]",
       "volume",
       "n_smooth",
       "N",
       "ptot",
       "Fn",
       "sum_MMMM[0]",
       "sum_MMMM[1]",
       "90-sum_MMMM[2]",
       "sum_MMMM[3]",
       "sum_MMMM[4]",
       "sum_MMMM[5]",
       "model",
       "xx",
       "PIJ[3]",
       "Q[1]",
        "F2",
        "divU",
        "100-UU",
        "UU0",
        "UUR",
        "UU0R",
        "alpha1",
        "alpha2",
        "alpha3",
        "alpha4",
        "alpha5",
        "alpha6",
        "110-alpha7",
        "alpha8",
        "alpha9",
        "beta1",
        "beta2",
        "beta3",
        "gamma_st",
        "phi[0]",
        "phi[1]",
        "phi[2]",
        "p_m",
        "p_rep"
          );
 int Nx,Ny;
 if((*box).problem == dcones){
 	Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
 	Ny = (*box).Ny;
 }
 else if((*box).problem == flatnose){
 	Nx = (*box).N1[0]+(*box).N1[1];
 	Ny = (*box).N2[0]+(*box).N2[1];
 }
 else{
 	Nx = (*box).N[0];
 	Ny = (*box).N[1];
 }




for( i = 0; i < Nx ;i++){
 	for( j = 0; j < Ny ;j++){
 	  num =  j*Nx + i;
 	  if(cells[num].in ==1){
 	      fprintf(fp, "%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20d%20e%20e%20e%20e%20e%20e%20e%20e%20d%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e%20e\n",
 		  cells[num].cell_center[0],
 		  cells[num].cell_center[1],
                   cells[num].dx,
                   cells[num].dy,
                   cells[num].dim[5]-cells[num].dim[4],
                   (*gas).t_post,
                   (*gas).delta_t,
 		      cells[num].U_space[0],
 		      cells[num].U_space[1],
 		      cells[num].U_space[2],
 		      cells[num].T,
 		      cells[num].n,
                   ///////  cell
                   cells[num].sum_weight,
                   cells[num].sum_M[0],
                   cells[num].sum_M[1],
                   cells[num].sum_M[2],
                   cells[num].sum_MM[0],
                   cells[num].sum_MM[1],
                   cells[num].sum_MM[2],
                   cells[num].sum_MM[3],
                   cells[num].sum_MM[4],
                   cells[num].sum_MM[5],
                   cells[num].sum_MMM[0],
                   cells[num].sum_MMM[1],
                   cells[num].sum_MMM[2],
                   cells[num].sum_MMM[3],
                   cells[num].sum_MMM[4],
                   cells[num].sum_MMM[5],
                   cells[num].sum_MMM[6],
                   cells[num].sum_MMM[7],
                   cells[num].sum_MMM[8],
                   cells[num].sum_MMM[9],
                   /// face
                   cells[num].sum_T_f[0],//((*gas).m/(3.0*(*gas).kb))*cells[num].T_f[0]/cells[num].M_n[0],
                   cells[num].sum_T_f[1],//((*gas).m/(3.0*(*gas).kb))*cells[num].T_f[1]/cells[num].M_n[1],
                   cells[num].sum_T_f[2],//((*gas).m/(3.0*(*gas).kb))*cells[num].T_f[2]/cells[num].M_n[2],
                   cells[num].sum_T_f[3],//((*gas).m/(3.0*(*gas).kb))*cells[num].T_f[3]/cells[num].M_n[3],
                   cells[num].sum_M_n[0],
                   cells[num].sum_M_n[1],
                   cells[num].sum_M_n[2],
                   cells[num].sum_M_n[3],
                   cells[num].sum_M1_f[0],
                   cells[num].sum_M1_f[1],
                   cells[num].sum_M1_f[2],
                   cells[num].sum_M1_f[3],
                   cells[num].sum_M1_f[4],
                   cells[num].sum_M1_f[5],
                   cells[num].sum_M2_f[0],
                   cells[num].sum_M2_f[1],
                   cells[num].sum_M2_f[2],
                   cells[num].sum_M2_f[3],
                   cells[num].sum_M2_f[4],
                   cells[num].sum_M2_f[5],
                   cells[num].sum_M3_f[0],
                   cells[num].sum_M3_f[1],
                   cells[num].sum_M3_f[2],
                   cells[num].sum_M3_f[3],
                   cells[num].sum_N[2],
                   cells[num].sum_N[3],
                   cells[num].sum_M_np[0],
                   cells[num].sum_M_np[1],
                   cells[num].sum_M_np[2],
                   cells[num].sum_M_np[3],
                   cells[num].sum_M_np[4],
                   cells[num].sum_M_np[5],
                   cells[num].sum_MM_f_1[0],
                   cells[num].sum_MM_f_1[1],
                   cells[num].sum_MM_f_1[2],
                   cells[num].sum_MM_f_1[3],
                   cells[num].sum_MM_f_1[4],
                   cells[num].sum_MM_f_1[5],
                   cells[num].sum_MM_f_2[0],
                   cells[num].sum_MM_f_2[1],
                   cells[num].sum_MM_f_2[2],
                   cells[num].sum_MM_f_2[3],
                   cells[num].sum_MM_f_2[4],
                   cells[num].sum_MM_f_2[5],
                   cells[num].sum_MM_f_3[0],
                   cells[num].sum_MM_f_3[1],
                   cells[num].sum_MM_f_3[2],
                   cells[num].sum_MM_f_3[3],
                   cells[num].sum_MM_f_3[4],
                   cells[num].sum_MM_f_3[5],
                   cells[num].volume,
 		  cells[num].n_smooth,
 		  cells[num].num_inside,
 		  cells[num].n*(*gas).kb*cells[num].T*( 1.0+cells[num].n*(*gas).b*increase_collision_rate(cells[num].n,gas) ) - 2.0*(*gas).b*cells[num].n*cells[num].n*(*gas).phi,\
 		  cells[num].n_ratio*(*gas).Fn,
       cells[num].sum_MMMM[0],
       cells[num].sum_MMMM[1],
       cells[num].sum_MMMM[2],
       cells[num].sum_MMMM[3],
       cells[num].sum_MMMM[4],
       cells[num].sum_MMMM[5],
       cells[num].model,
       cells[num].xx,
       cells[num].PIJ[3],
       cells[num].Q[1],
       cells[num].F2,
       cells[num].divU,
       cells[num].UU,
       cells[num].UU0,
       cells[num].UUR,
       cells[num].UU0R,
       cells[num].alpha[0],
       cells[num].alpha[1],
       cells[num].alpha[2],
       cells[num].alpha[3],
       cells[num].alpha[4],
       cells[num].alpha[5],
       cells[num].alpha[6],
       cells[num].alpha[7],
       cells[num].alpha[8],
       cells[num].beta[0],
       cells[num].beta[1],
       cells[num].beta[2],
       cells[num].gamma_st,
       cells[num].phi[0],
       cells[num].phi[1],
       cells[num].phi[2],
       cells[num].p_m,
       cells[num].n*(*gas).kb*cells[num].T*(1.0+cells[num].n*(*gas).b*increase_collision_rate(cells[num].n,gas))
  	      );
 	  }
         }
         fprintf(fp,"\n");
  }

  fclose(fp);



     double totalc = (*gas).count_SP_c + (*gas).count_reem_c + (*gas).count_evap_c;
     double totalh = (*gas).count_SP_h + (*gas).count_reem_h + (*gas).count_evap_h;
     //if(totalc>1e-15 && totalh > 1e-15 && totalc-(*gas).count_evap_c > 1e-15 && totalh-(*gas).count_evap_h > 1e-15){

     if(fabs((*gas).N_abs)>1e-15 && fabs((*gas).t_post)>1e-15 ){
       char name_output2[100];
       snprintf(name_output2, 100, "PostProc/%s.dat",a);
       fp = fopen(name_output2, "a");
       fprintf(fp,"%d  ",step);
       fprintf(fp,"%20e%20e%20e", (*gas).N_abs, (*gas).t_post, (*gas).N_abs/ (*gas).t_post);
       fprintf(fp,"\n");
       fclose(fp);

       printf("also in %s\n",name_output2);
     }

     /*
     if(fabs((*box).JR[0])>1e-15){
 	     fprintf(fp, "%20e", (*box).JR[1]/(*box).JR[0]);
     	  for(i=0; i<4; i++){
 		         fprintf(fp, "%20e%20e%20e%20e", (*box).JR[i], (*box).NR[i], (*box).JL[i],(*box).NL[i]);
     	  }
     	  fprintf(fp,"\n");
     	  fclose(fp);
     }
     */

     if (step % (*box).every == 0 || step == (*box).init_step)
 	     printf("Output written in PostProc/%s_%s_%d.dat\n",a,b,step);
     }


}

void post_processing(int step, struct CELLS *cells,  struct BOX *box, int *done, struct GAS *gas, double *U1, double *U2, double *U3, double *x1, double *x2, int *flag, int *index, int *color, int*n_ratio,double *x2_old){
  int i,j,jj,k,iii,num;
  FILE *fp;
//  double times = 4.0;
  double dummy=0.0,Y,mu;
double start, end;
//start = omp_get_wtime();




	if ((*box).step % (*box).every == 0){
              double T_avg = 0.0;
              for(num=0; num<(*box).N[0]*(*box).N[1]*(*box).N[2]; num++)
                  T_avg += cells[num].T;
              T_avg = T_avg/((*box).N[0]*(*box).N[1]*(*box).N[2]);
             if((*box).N[1]>1){
                   printf("T   = %e,   %e,   %e.... %e ....  %e   %e\n",cells[0].T, cells[1].T,cells[2].T, cells[(*box).N[1]/2].T, cells[(*box).N[1]-2].T, cells[(*box).N[1]-1].T);
                   printf("N   = %d,   %d,  %d .... %d ....  %d   %d\n",cells[0].num_inside, cells[1].num_inside,cells[2].num_inside, cells[(*box).N[1]/2].num_inside, cells[(*box).N[1]-2].num_inside, cells[(*box).N[1]-1].num_inside);
		   printf("n*sigma^3 = %e,   %e,   %e.... %e ....  %e   %e\n",cells[0].n*pow((*gas).sigma,3), cells[1].n*pow((*gas).sigma,3),cells[2].n*pow((*gas).sigma,3), cells[(*box).N[1]/2].n*pow((*gas).sigma,3), cells[(*box).N[1]-2].n*pow((*gas).sigma,3), cells[(*box).N[1]-1].n*pow((*gas).sigma,3));
                   printf("p   = %e,   %e .... %e ....  %e   %e\n",(*gas).m*cells[0].n*cells[0].DM2/3.0, (*gas).m*cells[1].n*cells[1].DM2/3.0, (*gas).m*cells[(*box).N[1]/2].n*cells[(*box).N[1]/2].DM2/3.0, (*gas).m*cells[(*box).N[1]-2].n*cells[(*box).N[1]-1].DM2/3.0, (*gas).m*cells[(*box).N[1]-2].n*cells[(*box).N[1]-1].DM2/3.0);
                   printf("T_avg = %e\n", T_avg);
                   	printf("(*gas).Fn=%e\n",(*gas).Fn);
             }
             else{
                   printf("T   = %e\n",cells[0].T);
                   printf("N   = %d\n",cells[0].num_inside);
                   printf("p   = %e\n",(*gas).m*cells[0].n*cells[0].DM2/3.0);
                   //printf("potential energy = %e\n", sum(phi,gas.N));
                   printf("kinetic energy = %e\n", 0.5*(*gas).m*cells[0].num_inside*cells[0].DM2);
                   if((*box).problem == dcones)
			printf("N7=%d, N6=%d, N5=%d, N4=%d\n",(*box).Ndot_w7,(*box).Ndot_w6,(*box).Ndot_w5,(*box).Ndot_w4);
             }
         }





  if ( step > (*gas).times*(*box).after){
      (*gas).t_post += (*gas).delta_t;

      cell_poly_velocity(U1, U2, U3, cells,  box, step-(*gas).times*(*box).after, gas, index);
      //measure_pressure(U2, x2, x2_old, gas, box, cells);
      //compute_evap_coeff(U2, x2, x2_old, gas, box, color, n_ratio);
      for(i=0; i<(*box).N[0]*(*box).N[1]*(*box).N[2]; i++){


         dummy = cells[i].n*(*gas).m*cells[i].PIJ[1];

         Y = increase_collision_rate(cells[i].n, gas);
         mu = (*gas).mu*pow(cells[i].T/(*gas).T0,(*gas).visp);

         //////////////////////////

        }
     }




//##########################   Incrementing time        ##########################

    if (step % (*box).every == 0 && (*box).problem != dcones && (*box).problem != flatnose){
       if(*done == 0 && (*box).reset !=1){
         *done = 1;
         for(i=0; i<(*gas).N; i++){
             if(flag[i] == 0){
	    	*done=0;
		break;
	     }
         }
         if(*done == 1){
		(*box).after = step;
		//write_pos_vel(x1, x2, x3, U1, U2, U3, (*gas).N);
		printf("______>     positions and velocities are written in vel_pos.txt\n");
         }
       }
       if(*done==1)
	 printf("both walls are already hit hit by all particles\n");
    }




           /////////////////////////////////////////////////////////////////////////////////
           //##########################       Sample pressure     ##########################
           /////////////////////////////////////////////////////////////////////////////////
           /*
           if ( (*box).step > (*gas).times*(*box).after && (*box).problem != dcones && (*box).problem != flatnose){
           	for( i=0; i<6; i++){
	       	      (*box).dummy_nom += (*gas).M[i];
             	}
             	(*gas).p = (*box).dummy_nom/( (*box).area*(*gas).avgeraging_time_till_now );
             	if ( (*box).step % (*box).every == 0){
			//for( i=0; i<6; i++)
			//	printf("(*gas).M[%d]=%lf",i,(*gas).M[i]);
                	printf("p = %e\n", (*gas).p);
             	}
           }
           */


         //writting_position_to_file(x1, x2, x3, gas.N, &t, 1, step);
         //generating_Python_paraview_script(step+1, &box);

//	end = omp_get_wtime();
//	if ( (*box).step % (*box).every == 0)
//		printf("post-proc took %e s\n", end - start);
}
void making_inputs(struct GAS *gas,struct BOX *Box)
{
  //// some constants
  double pi = acos(-1.0);

  ////	BOX Data
  double scale = 1.0;
  // BOX Dimensions
  (*Box).Len[0] = scale*100;
  (*Box).Len[1] = scale*100;
  (*Box).Len[2] = scale*100;

  (*Box).delta_dim[0] = (*Box).Len[0]/1.0;
  (*Box).delta_dim[1] = (*Box).Len[1]/1.0;
  (*Box).delta_dim[2] = (*Box).Len[2]/1.0;

  (*Box).N[0] = (int )( (*Box).Len[0]/ (*Box).delta_dim[0]);
  (*Box).N[1] = (int )( (*Box).Len[1]/ (*Box).delta_dim[1]);
  (*Box).N[2] = (int )( (*Box).Len[2]/ (*Box).delta_dim[2]);

  (*Box).Area[0] = (*Box).delta_dim[1]*(*Box).delta_dim[2];
  (*Box).Area[1] = (*Box).delta_dim[1]*(*Box).delta_dim[2];
  (*Box).Area[2] = (*Box).delta_dim[0]*(*Box).delta_dim[2];
  (*Box).Area[3] = (*Box).delta_dim[0]*(*Box).delta_dim[2];
  (*Box).Area[4] = (*Box).delta_dim[0]*(*Box).delta_dim[1];
  (*Box).Area[5] = (*Box).delta_dim[0]*(*Box).delta_dim[1];

  (*Box).Volume = (*Box).Len[0]*(*Box).Len[1]*(*Box).Len[2];


  ///Gas Data
  (*gas).kb = 1.38064852*1e-23;
  (*gas).sigma = 1;
  (*gas).T = 1.0/(*gas).kb;
  (*gas).n = 1e-2;
  (*gas).w = 1.0/20;
  (*gas).m = 1.0;
  (*gas).N = (int) ( (*gas).n*(*Box).Volume )+1;
  (*gas).lambda = 1.0/(sqrt(2)*pi*(*gas).sigma*(*gas).sigma*(*gas).n);
  (*gas).Kn = (*gas).lambda/(*Box).Len[1];
  (*gas).delta_t =  0.005*(*gas).lambda/sqrt((*gas).kb*(*gas).T/(*gas).m);

  (*gas).U0[0] = 1.0;
  (*gas).U0[1] = 1.0;
  (*gas).U0[2] = 1.0;

  printf("N = %ld\n", (*gas).N);
  printf("Lambda = %lf\n", (*gas).lambda);
};
void make_cells(struct BOX *Box, struct CELLS *cells, struct GAS *gas){
  int j,k, l,kk,jj;
  int num = 0;

int ii;
(*Box).out = 0.0;
int Nx = (*Box).Nx[0]+(*Box).Nx[1]+(*Box).Nx[2];
int Ny;
double dx0, dy0;
double dx1, dy1;
double Ly;
dy0 = (*Box).Ly/geometric_sum((*Box).Ny-1,(*Box).s);
(*Box).dy0 = dy0;
if((*Box).problem == dcones){
	for (k=0; k < (*Box).Ny; k++){
		for(j=0; j<Nx; j++){
			num = k*Nx + j;
			cells[num].in = 1;
			if(j<(*Box).Nx[0]){
// first dx and dy
				cells[num].dx = (*Box).Lx[0]/(1.0*(*Box).Nx[0]);
				cells[num].dy = dy0*pow((*Box).s,k);
// dim here is y of 4 corners of the cells
		        	cells[ num ].dim[0] = dy0*geometric_sum(k-1,  (*Box).s);
				cells[ num ].dim[1] = dy0*geometric_sum(k-1,  (*Box).s);
		        	cells[ num ].dim[2] = dy0*geometric_sum(k,(*Box).s);
				cells[ num ].dim[3] = dy0*geometric_sum(k,(*Box).s);
// cell centers needed for postprocessing
          			cells[num].cell_center[0] = (j+0.5)*cells[num].dx;
          			cells[num].cell_center[1] = 0.5*(cells[num].dim[0]+cells[num].dim[2]);
			}
			else if(j >= (*Box).Nx[0] && j < (*Box).Nx[0]+(*Box).Nx[1]){
// first dx and dy
				cells[num].dx = (*Box).Lx[1]/(1.0*(*Box).Nx[1]);
				cells[num].dy = dy0*pow((*Box).s,k);
// dim here is y of 4 corners of the cells
		        	cells[ num ].dim[0] = (j-(*Box).Nx[0])  *cells[num].dx*tan((*Box).alpha[0])+dy0*geometric_sum(k-1,  (*Box).s);
				cells[ num ].dim[1] = (j-(*Box).Nx[0]+1)*cells[num].dx*tan((*Box).alpha[0])+dy0*geometric_sum(k-1,  (*Box).s);
		        	cells[ num ].dim[2] = (j-(*Box).Nx[0])  *cells[num].dx*tan((*Box).alpha[0])+dy0*geometric_sum(k,(*Box).s);
				cells[ num ].dim[3] = (j-(*Box).Nx[0]+1)*cells[num].dx*tan((*Box).alpha[0])+dy0*geometric_sum(k,(*Box).s);
// cell centers needed for postprocessing
          			cells[num].cell_center[0] = (*Box).Lx[0] + (j-(*Box).Nx[0] +0.5)*cells[num].dx;
          			cells[num].cell_center[1] = 0.5*(cells[num].dim[0]+cells[num].dim[3]);
			}
			else{
// first dx and dy
				cells[num].dx = (*Box).Lx[2]/(1.0*(*Box).Nx[2]);
				cells[num].dy = dy0*pow((*Box).s,k);
// dim here is y of 4 corners of the cells
				Ly = (*Box).Lx[1]*tan((*Box).alpha[0]);
		        	cells[ num ].dim[0] = Ly+(j-(*Box).Nx[1]-(*Box).Nx[0])  *cells[num].dx*tan((*Box).alpha[1])+dy0*geometric_sum(k-1,  (*Box).s);
				cells[ num ].dim[1] = Ly+(j-(*Box).Nx[1]-(*Box).Nx[0]+1)*cells[num].dx*tan((*Box).alpha[1])+dy0*geometric_sum(k-1,  (*Box).s);
		        	cells[ num ].dim[2] = Ly+(j-(*Box).Nx[1]-(*Box).Nx[0])  *cells[num].dx*tan((*Box).alpha[1])+dy0*geometric_sum(k,(*Box).s);
				cells[ num ].dim[3] = Ly+(j-(*Box).Nx[1]-(*Box).Nx[0]+1)*cells[num].dx*tan((*Box).alpha[1])+dy0*geometric_sum(k,(*Box).s);
// cell centers needed for postprocessing
          			cells[num].cell_center[0] = (*Box).Lx[0] + (*Box).Lx[1] + (j-(*Box).Nx[0]-(*Box).Nx[1] +0.5)*cells[num].dx;
          			cells[num].cell_center[1] = 0.5*(cells[num].dim[0]+cells[num].dim[3]);
			}
////// Initializing every other elemet of cells to zero
	   cells[ num ].dim[4] = 0.0;
	   cells[ num ].dim[5] = 0.0;
	   cells[ num ].ng[0] =0;
	   cells[ num ].ng[1] =0;
	   cells[ num ].volume = volume_conical_frustum(cells[num].dim[2],cells[num].dim[3],cells[num].dx) - volume_conical_frustum(cells[num].dim[0],cells[num].dim[1],cells[num].dx);
	   cells[ num ].indices_inside = NULL;
           cells[ num ].num_inside = 0;
	   cells[ num ].weight = 0.0;
	   cells[ num ].F2 = 0.0;
	   cells[ num ].dn2 = 0.0;
	   cells[ num ].dn3 = 0.0;
	   cells[ num ].n_smooth = 0.0;
     cells[num].UU = 0.0;
     cells[num].UU0 = 0.0;
     cells[num].UUR = 0.0;
     cells[num].UU0R = 0.0;
     cells[num].p_m = 0.0;
           for(ii=0; ii<9; ii++)
		cells[ num ].alpha[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].beta[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].U_space[ii] = 0.0;
  	   cells[ num ].T = 0.0;
  	   cells[ num ].n = 0.0;
         for(ii=0; ii<10; ii++)
            cells[ num ].sum_MMM[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].U[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].PIJ[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].M4[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].Q[ii] = 0.0;

  	   for(ii=0; ii<10; ii++)
		cells[ num ].M3[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].c[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].gamma[ii] = 0.0;

        cells[ num ].G_sum = 0.0;
  	  cells[ num ].DM2 = 0.0;
  	  cells[ num ].DM4 = 0.0;
  	  cells[ num ].crm = 0.0;
	  cells[ num ].omega_max = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M_n[ii] = -1.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M_n[ii] = 0.0;
    for(ii=0; ii<6; ii++)
  		cells[ num ].sum_M_np[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M1_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M1_f[ii] = 0.0;
  for(ii=0; ii<4; ii++)
    cells[num].sum_N[ii] = 0.0;
    for(ii=0; ii<4; ii++)
  		cells[ num ].sum_T_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M2_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
            cells[ num ].sum_M3_f[ii] = 0.0;
      for(ii=0; ii<6; ii++)
		cells[ num ].M2_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M3_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].T_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].sum_MM[ii] = 0.0;
      for(ii=0; ii<3; ii++)
		cells[ num ].sum_M[ii] = 0.0;

   cells[num].sum_weight = 0.0;
   cells[num].xx = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_0[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_1[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_2[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_3[ii] = 0.0;


   for(ii=0; ii<3; ii++)
       cells[num].phi[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].Mp[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].M5[ii]=0.0;
   cells[num].gamma_st = 0.0;
		}
	}
	printf("done!\n");
}
else if((*Box).problem == flatnose){
	dx0 = (*Box).L1[0]/geometric_sum((*Box).N1[0]-1,(*Box).s1);
	dx1 = (*Box).L1[1]/geometric_sum((*Box).N1[1]-1,2.0-(*Box).s1);
	dy0 = (*Box).L2[0]/geometric_sum((*Box).N2[0]-1,2.0-(*Box).s2);
	dy1 = (*Box).L2[1]/geometric_sum((*Box).N2[1]-1,(*Box).s2);
	(*Box).dx0 = dx0;
	(*Box).dx1 = dx1;
	(*Box).dy0 = dy0;
	(*Box).dy1 = dy1;
	Nx = (*Box).N1[0]+(*Box).N1[1];
	Ny = (*Box).N2[0]+(*Box).N2[1];
	double dummy = (*Box).L1[0];
	for (k=0; k < Ny; k++){
		for(j=0; j<Nx; j++){
			num = k*Nx + j;
			if(k<(*Box).N2[0] && j>=(*Box).N1[0])
				cells[num].in = 0;
			else
				cells[num].in = 1;
			if( j<(*Box).N1[0] ){
				cells[num].dx = dx0*pow((*Box).s1, j);
				cells[num].cell_center[0] = dx0*geometric_sum(j-1,(*Box).s1)+cells[num].dx/2.0;
				//cells[num].dx = (*Box).L1[0]/(1.0*(*Box).N1[0]);
				//cells[num].cell_center[0] = (j+0.5)*cells[num].dx;
			}
			else{
				jj = j-(*Box).N1[0];
				cells[num].dx = dx1*pow(2.0-(*Box).s1, jj);
				cells[num].cell_center[0] = (*Box).L1[0] +  dx1*geometric_sum(jj-1,2.0-(*Box).s1)+cells[num].dx/2.0;
/*
				cells[num].dx = (*Box).L1[1]/(1.0*(*Box).N1[1]);
				cells[num].cell_center[0] = (*Box).L1[0] + (j-(*Box).N1[0] +0.5)*cells[num].dx;
*/
			}
			if( k<(*Box).N2[0] ){
				cells[num].dy = dy0*pow(2.0-(*Box).s2,k);
		        	cells[ num ].dim[0] =  dy0*geometric_sum(k-1,  2.0-(*Box).s2);
				cells[ num ].dim[1] =  dy0*geometric_sum(k-1,  2.0-(*Box).s2);
		        	cells[ num ].dim[2] =  dy0*geometric_sum(k,    2.0-(*Box).s2);
				cells[ num ].dim[3] =  dy0*geometric_sum(k,    2.0-(*Box).s2);
/*
				cells[num].dy = (*Box).L2[0]/(1.0*(*Box).N2[0]);
				cells[ num ].dim[0] = k*cells[num].dy;
				cells[ num ].dim[1] = k*cells[num].dy;
				cells[ num ].dim[2] = (k+1)*cells[num].dy;
				cells[ num ].dim[3] = (k+1)*cells[num].dy;
*/
				cells[num].cell_center[1] = 0.5*(cells[num].dim[0]+cells[num].dim[2]);
			}
			else{
				kk = k-(*Box).N2[0];
				cells[num].dy = dy1*pow((*Box).s2,kk);
		        	cells[ num ].dim[0] = (*Box).L2[0] + dy1*geometric_sum(kk-1,  (*Box).s2);
				cells[ num ].dim[1] = (*Box).L2[0] + dy1*geometric_sum(kk-1,  (*Box).s2);
		        	cells[ num ].dim[2] = (*Box).L2[0] + dy1*geometric_sum(kk,(*Box).s2);
				cells[ num ].dim[3] = (*Box).L2[0] + dy1*geometric_sum(kk,(*Box).s2);

//				cells[num].dy = (*Box).L2[1]/(1.0*(*Box).N2[1]);
//				cells[ num ].dim[0] = (*Box).L2[0]+(k-(*Box).N2[0])*cells[num].dy;
//				cells[ num ].dim[1] = (*Box).L2[0]+(k-(*Box).N2[0])*cells[num].dy;
//				cells[ num ].dim[2] = (*Box).L2[0]+(k-(*Box).N2[0]+1)*cells[num].dy;
//				cells[ num ].dim[3] = (*Box).L2[0]+(k-(*Box).N2[0]+1)*cells[num].dy;
				cells[num].cell_center[1] = 0.5*(cells[num].dim[0]+cells[num].dim[2]);
			}

////// Initializing every other elemet of cells to zero
	   if(dummy > cells[num].dx)
		dummy = cells[num].dx;
	   if(dummy > cells[num].dy)
		dummy = cells[num].dy;
	   cells[ num ].dim[4] = 0.0;
	   cells[ num ].dim[5] = 0.0;
	   cells[ num ].ng[0] =0;
	   cells[ num ].ng[1] =0;
	   cells[ num ].volume = acos(-1.0)*(cells[num].dim[2]*cells[num].dim[2]-cells[num].dim[0]*cells[num].dim[0])*cells[num].dx;
	   cells[ num ].indices_inside = NULL;
           cells[ num ].num_inside = 0;
	   cells[ num ].weight = 0.0;
	   cells[ num ].F2 = 0.0;
	   cells[ num ].dn2 = 0.0;
	   cells[ num ].dn3 = 0.0;
	   cells[ num ].n_smooth = 0.0;
     cells[num].UU = 0.0;
     cells[num].UU0 = 0.0;
     cells[num].UUR = 0.0;
     cells[num].UU0R = 0.0;
     cells[num].p_m = 0.0;
           for(ii=0; ii<9; ii++)
		cells[ num ].alpha[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].beta[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].U_space[ii] = 0.0;
  	   cells[ num ].T = 0.0;
  	   cells[ num ].n = 0.0;
       cells[ num ].divU = 0.0;
         for(ii=0; ii<10; ii++)
            cells[ num ].sum_MMM[ii] = 0.0;
            for(ii=0; ii<6; ii++)
                cells[ num ].sum_MMMM[ii] = 0.0;
  	   for(ii=0; ii<3; ii++)
		cells[ num ].U[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].PIJ[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].M4[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].Q[ii] = 0.0;

  	   for(ii=0; ii<10; ii++)
		cells[ num ].M3[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].c[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].gamma[ii] = 0.0;

        cells[ num ].G_sum = 0.0;
  	  cells[ num ].DM2 = 0.0;
  	  cells[ num ].DM4 = 0.0;
  	  cells[ num ].crm = 0.0;
	  cells[ num ].omega_max = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M_n[ii] = -1.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M_n[ii] = 0.0;
    for(ii=0; ii<6; ii++)
  		cells[ num ].sum_M_np[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M1_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M1_f[ii] = 0.0;
    for(ii=0; ii<4; ii++)
      cells[num].sum_N[ii] = 0.0;
    for(ii=0; ii<4; ii++)
      cells[ num ].sum_T_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M2_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
            cells[ num ].sum_M3_f[ii] = 0.0;
      for(ii=0; ii<6; ii++)
		cells[ num ].M2_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M3_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].T_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].sum_MM[ii] = 0.0;
      for(ii=0; ii<3; ii++)
		cells[ num ].sum_M[ii] = 0.0;
   if(cells[num].in == 0)
	   cells[num].sum_weight = -1.0;
   else
	   cells[num].sum_weight = 0.0;
   cells[num].xx = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_0[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_1[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_2[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_3[ii] = 0.0;


   for(ii=0; ii<3; ii++)
       cells[num].phi[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].Mp[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].M5[ii]=0.0;
   cells[num].gamma_st = 0.0;
		}
	}
        (*gas).delta_t =  (*gas).frac_mft*min(  (*gas).lambda, dummy)/sqrt( (*gas).kb*(*Box).T_wall_5/(*gas).m );
	printf("dt = %e\n",(*gas).delta_t);
	printf("done!\n");
}
else{
  for (l=0; l < (*Box).N[2]; l++)
    for (k=0; k < (*Box).N[1]; k++)
      for (j=0; j < (*Box).N[0]; j++){


	    num = l*(*Box).N[0]*(*Box).N[1] + k*(*Box).N[0] + j;
      cells[num].id[0] = k*((*Box).N[0]+1) + j;
      cells[num].id[1] = k*((*Box).N[0]+1) + j + 1;
      cells[num].id[2] = (k+1)*((*Box).N[0]+1) + j;
      cells[num].id[3] = (k+1)*((*Box).N[0]+1) + j + 1;



	    cells[num].in = 1;
            cells[ num ].dim[0] = (double)(j)  *(*Box).delta_dim[0];
	    cells[ num ].dim[1] = (double)(j+1)*(*Box).delta_dim[0];
            cells[ num ].dim[2] = (double)(k)  *(*Box).delta_dim[1];
	    cells[ num ].dim[3] = (double)(k+1)*(*Box).delta_dim[1];
	    cells[ num ].dim[4] = (double)(l)  *(*Box).delta_dim[2];
	    cells[ num ].dim[5] = (double)(l+1)*(*Box).delta_dim[2];
	    cells[num].dx = cells[ num ].dim[1] - cells[ num ].dim[0];
	    cells[num].dy = cells[ num ].dim[3] - cells[ num ].dim[2];
	    cells[ num ].volume = (cells[num].dim[1]-cells[num].dim[0])*(cells[num].dim[3]-cells[num].dim[2])*(cells[num].dim[5]-cells[num].dim[4]);
          cells[num].cell_center[0] = 0.5*(cells[num].dim[0]+cells[num].dim[1]);
          cells[num].cell_center[1] = 0.5*(cells[num].dim[2]+cells[num].dim[3]);
          cells[num].cell_center[2] = 0.5*(cells[num].dim[4]+cells[num].dim[5]);
          cells[num].n_n[0] = 0.0;
          cells[num].n_n[1] = 0.0;
          cells[num].n_n[2] = 0.0;
          cells[num].n_n[3] = 0.0;
          cells[num].x_n[0] = cells[num].cell_center[0]-cells[num].dx/2.0;
          cells[num].x_n[1] = cells[num].cell_center[0]+cells[num].dx/2.0;
          cells[num].x_n[2] = cells[num].cell_center[0]-cells[num].dx/2.0;
          cells[num].x_n[3] = cells[num].cell_center[0]+cells[num].dx/2.0;
          cells[num].y_n[0] = cells[num].cell_center[1]-cells[num].dy/2.0;
          cells[num].y_n[1] = cells[num].cell_center[1]-cells[num].dy/2.0;
          cells[num].y_n[2] = cells[num].cell_center[1]+cells[num].dy/2.0;
          cells[num].y_n[3] = cells[num].cell_center[1]+cells[num].dy/2.0;
	   double a = 0.0;

	int sg1, sg3;
	if((*gas).LvN1>0)
		sg1 = 1;
	else
		sg1 = 0;

	if((*gas).LvN3>0)
		sg3 = 1;
	else
		sg3 = 0;

	   if( (*Box).problem == inverted && (*gas).model != "MD"){
		if( cells[num].cell_center[1] < (*gas).Lc_center-(*gas).Lc_ratio/2.0  && sg1 == 1){
			cells[ num ].n_ratio = 1;
		}
		else if( cells[num].cell_center[1] < (*gas).Lc_center+(*gas).Lc_ratio/2.0 ){
			cells[ num ].n_ratio = (*gas).n_ratio;
		}
		else if( cells[num].cell_center[1] < (*gas).Lh_center-(*gas).Lh_ratio/2.0 ){
			cells[ num ].n_ratio = 1;
		}
		else if( cells[num].cell_center[1] < (*gas).Lh_center+(*gas).Lh_ratio/2.0 ){
			cells[ num ].n_ratio = (*gas).n_ratio;
		}
		else if( sg1 == 1){
			cells[ num ].n_ratio = 1;
		}
		else if( sg1 == 0){
			cells[ num ].n_ratio = (*gas).n_ratio;
		}
	   }
	   else if( (*Box).problem == evaporation && (*gas).model != "MD"){
		if( cells[num].cell_center[1] < (*gas).Lc_center-(*gas).Lc_ratio/2.0  && sg1 == 1){
			cells[ num ].n_ratio = 1;
		}
		else if( cells[num].cell_center[1] < (*gas).Lc_center+(*gas).Lc_ratio/2.0 ){
			cells[ num ].n_ratio = (*gas).n_ratio;
		}
		else{
			cells[ num ].n_ratio = 1;
		}
	   }
	   else
	   	cells[ num ].n_ratio = 1;
	   cells[ num ].ng[0] =0;
	   cells[ num ].ng[1] =0;
	   cells[ num ].volume = (cells[num].dim[1]-cells[num].dim[0])*(cells[num].dim[3]-cells[num].dim[2])*(cells[num].dim[5]-cells[num].dim[4]);
	   cells[ num ].indices_inside = NULL;
           cells[ num ].num_inside = 0;
	   cells[ num ].weight = 0.0;
	   cells[ num ].F2 = 0.0;
	   cells[ num ].dn2 = 0.0;
	   cells[ num ].dn3 = 0.0;
	   cells[ num ].n_smooth = 0.0;
     cells[num].UU = 0.0;
     cells[num].UU0 = 0.0;
     cells[num].UUR = 0.0;
     cells[num].UU0R = 0.0;
     cells[num].p_m = 0.0;
           for(ii=0; ii<9; ii++)
		cells[ num ].alpha[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].beta[ii] = 0.0;
           for(ii=0; ii<3; ii++)
		cells[ num ].U_space[ii] = 0.0;
  	   cells[ num ].T = 0.0;
  	   cells[ num ].n = 0.0;
       cells[ num ].divU = 0.0;
        for(ii=0; ii<10; ii++)
            cells[ num ].sum_MMM[ii] = 0.0;
        for(ii=0; ii<6; ii++)
            cells[ num ].sum_MMMM[ii] = 0.0;
  	   for(ii=0; ii<3; ii++)
		cells[ num ].U[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].PIJ[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].M4[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].Q[ii] = 0.0;

  	   for(ii=0; ii<10; ii++)
		cells[ num ].M3[ii] = 0.0;

  	   for(ii=0; ii<6; ii++)
		cells[ num ].c[ii] = 0.0;

  	   for(ii=0; ii<3; ii++)
		cells[ num ].gamma[ii] = 0.0;

        cells[ num ].G_sum = 0.0;
  	  cells[ num ].DM2 = 0.0;
  	  cells[ num ].DM4 = 0.0;
  	  cells[ num ].crm = 0.0;
	  cells[ num ].omega_max = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M_n[ii] = -1.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M_n[ii] = 0.0;
    for(ii=0; ii<6; ii++)
  		cells[ num ].sum_M_np[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M1_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M1_f[ii] = 0.0;
    for(ii=0; ii<4; ii++)
      cells[num].sum_N[ii] = 0.0;
    for(ii=0; ii<4; ii++)
  		cells[ num ].sum_T_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
		cells[ num ].sum_M2_f[ii] = 0.0;
  for(ii=0; ii<6; ii++)
            cells[ num ].sum_M3_f[ii] = 0.0;
      for(ii=0; ii<6; ii++)
		cells[ num ].M2_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].M3_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].T_f[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		cells[ num ].sum_MM[ii] = 0.0;
      for(ii=0; ii<3; ii++)
		cells[ num ].sum_M[ii] = 0.0;

   cells[num].sum_weight = 0.0;
   cells[num].xx = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_0[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_1[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_2[ii] = 0.0;
   for(ii=0; ii<6; ii++)
		    cells[ num ].sum_MM_f_3[ii] = 0.0;


   for(ii=0; ii<3; ii++)
       cells[num].phi[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].Mp[ii]=0.0;
   for(ii=0; ii<3; ii++)
       cells[num].M5[ii]=0.0;
   cells[num].gamma_st = 0.0;
   if((*gas).model=="SPH" || (*gas).model=="Hybrid"){
     cells[num].model = sph;
   }
   else if( (*gas).model == "DSMC_VHS" ){
     cells[num].model = dsmc;
	}
}
}
}
void writting_position_to_file(double *x1, double *x2, double*x3, int N, double *t, int len_t, int step){
  int i, j;
  FILE *fp;
  char name_output[50];

  for (j=0; j < len_t; j++){
    snprintf(name_output, 50, "PostProc/Positions.vtk.%d", step);
    fp = fopen(name_output, "w");
    //fprintf(fp, "POINTS %d FLOAT\n",N);
    for (i=0; i < N; i++){
      fprintf(fp, "%20e%20e%20e\n", x1[i], x2[i], x3[i]);
    }
    fclose(fp);
  }
  snprintf(name_output, 50, "PostProc/tracking0.txt");
  fp = fopen(name_output, "a");
  fprintf(fp, "%20e%20e%20e%20d\n", x1[0], x2[0], x3[0], step);
  fclose(fp);
}
void writting_position_to_file2D(double *x1, double *x2, int N, double *t, int len_t, int step){
  int i, j;
  FILE *fp;
  char name_output[50];

  for (j=0; j < len_t; j++){
    snprintf(name_output, 50, "PostProc/Positions.vtk.%d", step);
    fp = fopen(name_output, "w");
    //fprintf(fp, "POINTS %d FLOAT\n",N);
    for (i=0; i < N; i++){
      fprintf(fp, "%20e%20e\n", x1[i], x2[i]);
    }
    fclose(fp);
  }
  snprintf(name_output, 50, "PostProc/tracking0.txt");
  fp = fopen(name_output, "a");
  fprintf(fp, "%20e%20e%20d\n", x1[0], x2[0], step);
  fclose(fp);
}
void writting_data_to_file(int num_steps, double delta_t, double *p, double* p_kin, double *T, struct GAS *gas, double* E_kin, int realization, double* p_mod){
  int j;
  FILE *fp;
  char name_output[50];
  snprintf(name_output, 50, "PostProc/data_real_%d.txt",realization);
  fp = fopen(name_output, "w");
  fprintf(fp, "%10s   %10s   %10s   %10s   %10s   %10s   %10s   %10s   %10s\n","#step", "time","p", "p_kin", "T", "p=n*k_b*T", "E_kin", "p_averaged", "p_mod");
  for (j=0; j < num_steps; j++){
    fprintf(fp, "%10d   %10lf   %10lf   %10lf   %10.3e   %10lf   %10lf   %10lf   %10lf\n",j+1, j*delta_t,( fabs(p[j*6+0])+fabs(p[j*6+1])+fabs(p[j*6+2])+fabs(p[j*6+3])+fabs(p[j*6+4])+fabs(p[j*6+0])  )/6, p_kin[j], T[j], double((*gas).n)*(*gas).kb*T[j], E_kin[j] , (*gas).p, p_mod[j] );
  }
  fclose(fp);
}
void writting_factor_to_file(double T, double n, double factor){
//   printf("T = %lf, n = %lf, factor = %lf",T, n, factor);
  FILE *fp;
  char name_output[50];
  snprintf(name_output, 50, "PostProc/data_real_%f.txt",n);
  fp = fopen(name_output, "a");
//   fprintf(fp, "%10s   %10s   %10s \n","#T", "n","factor");
//   for (j=0; j < num; j++){
    fprintf(fp, "%10lf   %10lf   %10lf \n",T, n, factor);
//   }
  fclose(fp);
}

void generating_Python_paraview_script(int num_step, struct BOX *Box){
  int  j;
  FILE *fp;
  char name_output[50];
  snprintf(name_output, 50, "ParaPython_Anim.py");
  fp = fopen(name_output, "w");
  fprintf(fp, "try: paraview.simple\nexcept: from paraview.simple import *\nparaview.simple._DisableFirstRenderCameraReset()\n\nPositions_vtk_ = ParticlesReader( FileName=[");
  for (j=1; j < num_step+1; j++){
    fprintf(fp, "'PostProc/Positions.vtk.%d'",j);
    if(j < num_step)
      fprintf(fp, ", ");
    if(j == num_step)
      fprintf(fp, "] )\n");
  }
fprintf(fp,"AnimationScene2 = GetAnimationScene()\n");
fprintf(fp,"AnimationScene1 = GetAnimationScene()\n");
  fprintf(fp,"AnimationScene2.EndTime = %d.0\n",num_step-1);
  fprintf(fp,"AnimationScene2.PlayMode = 'Snap To TimeSteps'\n");

  fprintf(fp,"AnimationScene1.EndTime = %d.0\n",num_step-1);
  fprintf(fp,"AnimationScene1.PlayMode = 'Snap To TimeSteps'\n");

  fprintf(fp,"RenderView2 = GetRenderView()\n");
  fprintf(fp,"a1_Scalar_PVLookupTable = GetLookupTableForArray( 'Scalar', 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )\n");

  fprintf(fp,"a1_Scalar_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )\n");

  fprintf(fp,"DataRepresentation2 = Show()\n");
  fprintf(fp,"DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]\n");
  fprintf(fp,"DataRepresentation2.SelectionPointFieldDataArrayName = 'Scalar'\n");
  fprintf(fp,"DataRepresentation2.ColorArrayName = ('POINT_DATA', 'Scalar')\n");
  fprintf(fp,"DataRepresentation2.LookupTable = a1_Scalar_PVLookupTable\n");
  fprintf(fp,"DataRepresentation2.ScaleFactor = 9.999584197998047\n");

  fprintf(fp,"RenderView2.CenterOfRotation = [49.997920989990234, 49.98411178588867, 49.82797622680664]\n");
  fprintf(fp,"a1_Scalar_PVLookupTable.ScalarOpacityFunction = a1_Scalar_PiecewiseFunction\n");
  fprintf(fp,"RenderView2.CameraViewUp = [0.058546596355159825, 0.9562916328121219, -0.2864936457737796]\n");

  fprintf(fp,"RenderView2.CameraPosition = [-22.86898207328301, 147.66798406367514, 360.99774167784255]\n");
//   fprintf(fp,"RenderView2.CameraPosition = [%lf, %lf, %lf]\n", 0.22*(*Box).Len[0], 1.5*(*Box).Len[1], 3.6*(*Box).Len[2]);

  fprintf(fp,"RenderView2.CameraClippingRange = [187.7422321793557, 519.3280615289876]\n");
  fprintf(fp,"RenderView2.CameraFocalPoint = [49.99792098999024, 49.98411178588867, 49.82797622680664]\n");
  fprintf(fp,"RenderView2.CameraParallelScale = 86.49295202040874\n");

  fprintf(fp,"DataRepresentation2.ColorArrayName = ('POINT_DATA', '')\n");

  fprintf(fp,"WriteAnimation('PostProc/anim.ogv', Magnification=1, Quality=2, FrameRate=15.000000)\n");


  fprintf(fp,"Render()");
  fclose(fp);
}
// #include <iostream>


void readSettingsFile(struct GAS *gas,struct BOX *box)
{
    double pi = acos(-1.0);
  (*box).Len[0] = -1.0;
  (*box).Len[1] = -1.0;
  (*box).Len[2] = -1.0;
  (*box).delta_dim[0] = -1.0;
  (*box).delta_dim[1] = -1.0;
  (*box).delta_dim[2] = -1.0;
  (*gas).kb = -1.0;
  (*gas).sigma = -1.0;
  (*gas).T = -1.0;
  (*gas).n = -1.0;
  (*gas).w = -1.0;
  (*gas).m = -1.0;
  (*gas).N = -1;
  (*gas).delta_t =  -1.0;
  (*gas).U0[0] = -1.0;
  (*gas).U0[1] = -1.0;
  (*gas).U0[2] = -1.0;
  (*box).num_steps = -1;
  (*box).num_realizations = -1;
  (*box).movie = -1;
  (*box).every = -1;
  (*box).after = -1;
  (*box).reset_every = -1;
  (*gas).frac_mft = -1;
  (*box).N[0] = -1;
  (*box).num_thr = -1;
  (*gas).Fn = -1;
  (*gas).mu = -1;
  (*gas).mu_corr = -1;
  (*gas).factor = -1.0;
  (*gas).phi0 = -1.0;
  (*box).problem = -1;

  (*box).U_wall_1 = -1.0;
  (*box).U_wall_2 = -1.0;
  (*box).U_wall_3 = -1.0;
  (*box).U_wall_4 = -1.0;
  (*box).U_wall_5 = -1.0;
  (*box).U_wall_6 = -1.0;
  (*box).U_wall_7 = -1.0;
  (*box).T_wall_1 = -1.0;
  (*box).T_wall_2 = -1.0;
  (*box).T_wall_3 = -1.0;
  (*box).T_wall_4 = -1.0;
  (*box).T_wall_5 = -1.0;
  (*box).T_wall_6 = -1.0;
  (*box).T_wall_7 = -1.0;
  (*box).n_wall_4 = -1.0;
  (*box).n_wall_5 = -1.0;
  (*box).n_wall_6 = -1.0;
  (*box).n_wall_7 = -1.0;
  (*box).reset = -1;
  (*box).direction[0] = 0;
  (*box).direction[1] = 0;
  (*box).direction[2] = 0;
  (*gas).output_name = "data.txt";
  (*gas).long_range =0;
  (*gas).epsilon = -1.0;//4.0*119.8*(*gas).kb;//4.47478200e-21;
  (*box).N_grid[0] = -1;
  (*box).N_grid[1] = -1;
  (*box).N_grid[2] = -1;

(*box).Nx[0] = -1;
(*box).Nx[1] = -1;
(*box).Nx[2] = -1;
(*box).Ny = -1;
(*box).s = 1.0;
(*box).Lx[0] = -1.0;
(*box).Lx[1] = -1.0;
(*box).Lx[2] = -1.0;
(*box).Ly = -1.0;
(*box).alpha[0] = -1.0;
(*box).alpha[1] = -1.0;
(*box).init_step = -1;
(*box).dcones_n = -1.0;
(*box).L1[0]=-1.0;
(*box).L1[1]=-1.0;
(*box).L2[0]=-1.0;
(*box).L2[1]=-1.0;
(*box).N1[0]=-1;
(*box).N1[1]=-1;
(*box).N2[0]=-1;
(*box).N2[1]=-1;
(*box).s1 = -1.0;
(*box).s2 = -1.0;
(*gas).N13 = 10;
(*gas).nvc = -1.0;
(*gas).nvh = -1.0;
(*gas).n_ratio = -1;
(*box).ghost = -1.0;
(*box).thermc = -1.0;
(*box).thermh = -1.0;
for(int i=0;i<4; i++)
  (*box).x2_prob[i]=-1.0;
(*gas).sigc = 0.9;
(*gas).alphac = 0.0;
(*gas).sigh = 0.9;
(*gas).alphah = 0.0;
(*gas).LvN1 = -1.0;
(*gas).LvN3 = -1.0;
(*gas).phi = -1.0;
(*gas).gp_nn =  GP;

for(int i=0;i<4;i++){
	(*box).JR[i] = 0.0;
	(*box).JL[i] = 0.0;
	(*box).NR[i] = 0.0;
	(*box).NL[i] = 0.0;
	(*box).Jref[i] = 0.0;
}
    FILE *inputFile;
    std::string dummyString;
    size_t len = 0;

    std::cout << std::endl;
    std::cout << "====== Settings =============================================================" << std::endl<<std::endl;

    char * cstr = new char [(*box).file_name.length()+1];
  std::strcpy (cstr, (*box).file_name.c_str());

    inputFile = fopen( cstr, "r");
//     inputFile = stdin;

//     int ch = getc(inputFile);

//    while((i = getchar()) != EOF)
          char *line = NULL;
     while ( getline(&line, &len, inputFile) != EOF)
    {
        // Get a line and store in line


      std::string s = line;
      if ( line[0] == '#')
	std::cout << s << "\n";
//    // If the first character of the line is not a '#'
      else // ( line[0] != '#')
      {
	    std::istringstream iss(s);
            iss >> dummyString;
	    if(dummyString == "N")
                 iss >> (*gas).N;
      else if(dummyString == "gp_nn")
                 iss >> (*gas).gp_nn;
	    else if(dummyString == "N13")
                 iss >> (*gas).N13;
      else if(dummyString == "sigc")
                 iss >> (*gas).sigc;
      else if(dummyString == "alphac")
                 iss >> (*gas).alphac;
      else if(dummyString == "sigh")
                iss >> (*gas).sigh;
      else if(dummyString == "alphah")
                iss >> (*gas).alphah;
	    else if(dummyString == "n")
                 iss >> (*gas).n;
	    else if(dummyString == "nc")
                 iss >> (*gas).nc;
	    else if(dummyString == "nh")
                 iss >> (*gas).nh;
      else if(dummyString == "nvc")
                 iss >> (*gas).nvc;
      else if(dummyString == "nvh")
                 iss >> (*gas).nvh;
	    else if(dummyString == "nv1")
                 iss >> (*gas).nv1;
	    else if(dummyString == "nv2")
                 iss >> (*gas).nv2;
	    else if(dummyString == "nv3")
                 iss >> (*gas).nv3;
	    else if(dummyString == "Tc")
                 iss >> (*gas).Tc;
	    else if(dummyString == "Th")
                 iss >> (*gas).Th;
	    else if(dummyString == "Tv1")
                 iss >> (*gas).Tv1;
	    else if(dummyString == "Tv2")
                 iss >> (*gas).Tv2;
	    else if(dummyString == "Tv3")
                 iss >> (*gas).Tv3;
	    else if(dummyString == "LcN")
                 iss >> (*gas).LcN;
	    else if(dummyString == "LhN")
                 iss >> (*gas).LhN;
	    else if(dummyString == "LvN")
                 iss >> (*gas).LvN;
	    else if(dummyString == "LvN1")
                 iss >> (*gas).LvN1;
	    else if(dummyString == "LvN3")
                 iss >> (*gas).LvN3;
                 else if(dummyString == "Lv1")
                            iss >> (*gas).Lv1;
           	    else if(dummyString == "Lv3")
                            iss >> (*gas).Lv3;
	    else if(dummyString == "Nh")
                 iss >> (*gas).Nh;
                 else if(dummyString == "Nv1")
                         iss >> (*gas).Nv1;
	    else if(dummyString == "Nv2")
                 iss >> (*gas).Nv2;
      else if(dummyString == "Nv3")
              iss >> (*gas).Nv3;
          else if(dummyString == "long_range")
                 iss >> (*gas).long_range;
	    else if(dummyString == "box_len"){
                 iss >> (*box).Len[0];
		 iss >> (*box).Len[1];
		 iss >> (*box).Len[2];}
	    else if(dummyString == "cone_len"){
                 iss >> (*box).Lx[0];
		 iss >> (*box).Lx[1];
		 iss >> (*box).Lx[2];
		 iss >> (*box).Ly;
		}
    else if(dummyString == "x2_prob"){
               iss >> (*box).x2_prob[0];
               iss >> (*box).x2_prob[1];
               iss >> (*box).x2_prob[2];
               iss >> (*box).x2_prob[3];
  }

	    else if(dummyString == "ghost")
                 iss >> (*box).ghost;
     else if(dummyString == "thermc")
                            iss >> (*box).thermc;
     else if(dummyString == "thermh")
           iss >> (*box).thermh;
	    else if(dummyString == "cone_num_cells"){
                 iss >> (*box).Nx[0];
		 iss >> (*box).Nx[1];
		 iss >> (*box).Nx[2];
		 iss >> (*box).Ny;
		 iss >> (*box).s;
		}
	    else if(dummyString == "cone_angles"){
                 iss >> (*box).alpha[0];
		 iss >> (*box).alpha[1];
		}
	    else if(dummyString == "flatnose"){
                 iss >> (*box).L1[0];
		 iss >> (*box).L1[1];
		 iss >> (*box).L2[0];
		 iss >> (*box).L2[1];
		 iss >> (*box).N1[0];
		 iss >> (*box).N1[1];
		 iss >> (*box).N2[0];
		 iss >> (*box).N2[1];
	    }
	   else if(dummyString == "flatratio"){
		 iss >> (*box).s1;
		 iss >> (*box).s2;
		}
	    else if(dummyString == "box_delta_dim"){
                 iss >> (*box).delta_dim[0];
		 iss >> (*box).delta_dim[1];
		 iss >> (*box).delta_dim[2];}
          else if(dummyString == "direction"){
                iss >> (*box).direction[0];
                iss >> (*box).direction[1];
                iss >> (*box).direction[2];}
	    else if(dummyString == "U0"){
                 iss >> (*gas).U0[0];
		 iss >> (*gas).U0[1];
		 iss >> (*gas).U0[2];}
	    else if(dummyString == "sigma")
		 iss >> (*gas).sigma;
	    else if(dummyString == "weight")
		 iss >> (*gas).w;
	    else if(dummyString == "dcones_n")
		 iss >> (*box).dcones_n;
	    else if(dummyString == "mass")
		 iss >> (*gas).m;
	    else if(dummyString == "steps")
		 iss >> (*box).num_steps;
	    else if(dummyString == "dt")
		iss >> (*gas).delta_t;
	    else if(dummyString == "init_step")
		iss >> (*box).init_step;
	    else if(dummyString == "model")
		iss >> (*gas).model;
          else if(dummyString == "output")
              iss >> (*gas).output_name;
	    else if(dummyString == "movie")
		iss >> (*box).movie;
	    else if(dummyString == "every")
		iss >> (*box).every;
	    else if(dummyString == "num_realizations")
		iss >> (*box).num_realizations;
	    else if(dummyString == "after")
		iss >> (*box).after;
    else if(dummyString == "reset_every")
    iss >> (*box).reset_every;
	    else if(dummyString == "frac_mft")
		iss >> (*gas).frac_mft;
	    else if(dummyString == "T")
		iss >> (*gas).T;
            else if(dummyString == "epsilon")
              iss >> (*gas).epsilon;
	    else if(dummyString == "kb")
		iss >> (*gas).kb;
	    else if(dummyString == "num_thr")
		iss >> (*box).num_thr;
	    else if(dummyString == "mu")
		iss >> (*gas).mu;
	    else if(dummyString == "visp")
		iss >> (*gas).visp;
	    else if(dummyString == "nualpha")
		iss >> (*gas).nualpha;
    else if(dummyString == "phi")
  iss >> (*gas).phi;
	    else if(dummyString == "phi0")
		iss >> (*gas).phi0;
	   else if(dummyString == "factor")
	        iss >> (*gas).factor;
	    else if(dummyString == "n_ratio")
                 iss >> (*gas).n_ratio;
	   else if(dummyString == "problem")
	        iss >> (*box).problem;
	   else if(dummyString == "U_wall_1")
	        iss >> (*box).U_wall_1;
	   else if(dummyString == "U_wall_2")
	        iss >> (*box).U_wall_2;
	   else if(dummyString == "U_wall_3")
	        iss >> (*box).U_wall_3;
	   else if(dummyString == "U_wall_4")
	        iss >> (*box).U_wall_4;
	   else if(dummyString == "U_wall_5")
	        iss >> (*box).U_wall_5;
	   else if(dummyString == "U_wall_6")
	        iss >> (*box).U_wall_6;
	   else if(dummyString == "U_wall_7")
	        iss >> (*box).U_wall_7;
	   else if(dummyString == "T_wall_1")
	        iss >> (*box).T_wall_1;
	   else if(dummyString == "T_wall_2")
	        iss >> (*box).T_wall_2;
	   else if(dummyString == "T_wall_3")
	        iss >> (*box).T_wall_3;
	   else if(dummyString == "T_wall_4")
	        iss >> (*box).T_wall_4;
	   else if(dummyString == "T_wall_5")
	        iss >> (*box).T_wall_5;
	   else if(dummyString == "T_wall_6")
	        iss >> (*box).T_wall_6;
	   else if(dummyString == "T_wall_7")
	        iss >> (*box).T_wall_7;
	   else if(dummyString == "n_wall_4")
	        iss >> (*box).n_wall_4;
	   else if(dummyString == "n_wall_5")
	        iss >> (*box).n_wall_5;
	   else if(dummyString == "n_wall_6")
	        iss >> (*box).n_wall_6;
	   else if(dummyString == "n_wall_7")
	        iss >> (*box).n_wall_7;
	   else if(dummyString == "reset")
	        iss >> (*box).reset;
    else if(dummyString == "drop_dia")
               iss >> (*box).drop_dia;

	    else if(dummyString == "num_cells"){
		iss >>  (*box).N[0];
		iss >>  (*box).N[1];
		iss >>  (*box).N[2];}
            else
            {
                std::cout << std::endl << "Unknown keyword in the settings file : "<<std::endl;
                std::cout << "Aborting...";
                exit(0);
            }
        }
    }
    fclose( inputFile);


   (*box).alpha[0] = pi/180.0*(*box).alpha[0];
   (*box).alpha[1] = pi/180.0*(*box).alpha[1];


    if ( (*box).delta_dim[0] < 0.0){
      (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
      (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
      (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];}
    if( (*box).movie < 0.0 )
      (*box).movie = 1;
    if( (*box).num_thr < 0 )
      (*box).num_thr = 1;
    if( (*box).after < 0.0 )
      (*box).after = 1;
      if( (*box).reset_every < 0.0 )
        (*box).reset_every = 1;
    if( (*box).num_realizations < 0.0 )
      (*box).num_realizations= 1;
    if( (*gas).kb < 0.0 )
      (*gas).kb= 1.38064852*1e-23;
    if((*gas).T < 0.0)
      (*gas).T = 1.0/(*gas).kb;
    if ((*gas).sigma < 0.0)
      (*gas).sigma = 1;
    if ( (*gas).n < 0.0)
      (*gas).n = 1e-3;
    if ( (*gas).w < 0.0 )
      (*gas).w = 1.0/20;
    if ( (*gas).m < 0.0 )
      (*gas).m = 1.0;
    if ( (*box).Len[0] < 0.0 ){
      (*box).Len[0] = 1.0;
      (*box).Len[1] = 1.0;
      (*box).Len[2] = 1.0;}
    if((*box).dcones_n<0.0)
	(*box).dcones_n = (*gas).n;
    if((*box).init_step<0)
	(*box).init_step = 0;
    if ( (*gas).U0[0] < 0.0 ){
      (*gas).U0[0] = 1.0;
      (*gas).U0[1] = 1.0;
      (*gas).U0[2] = 1.0;}
    if ( (*box).num_steps < 0)
      (*box).num_steps =  10;
    if ( (*box).every < 0)
      (*box).every = 1;

    if ( (*box).problem < 0)
      (*box).problem = 1;
    if ( (*box).T_wall_1 < 0)
      (*box).T_wall_1 = 1.0;
    if ( (*box).T_wall_2 < 0)
      (*box).T_wall_2 = 1.0;
    if ( (*box).T_wall_3 < 0)
      (*box).T_wall_3 = 1.0;
    if ( (*box).T_wall_4 < 0)
      (*box).T_wall_4 = 1.0;
/*
  if( (*gas).N < 0 ){
    (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
    (*gas).N = (int) ( (*gas).n*(*box).Volume )+1;

    (*box).N[0] = (int )( (*box).Len[0]/ (*box).delta_dim[0]);
    (*box).N[1] = (int )( (*box).Len[1]/ (*box).delta_dim[1]);
    (*box).N[2] = (int )( (*box).Len[2]/ (*box).delta_dim[2]);
  }
  else{
//      (*box).Volume = (*gas).N/(*gas).n;
//      (*box).Len[0] = cbrt((*box).Volume);
//      (*box).Len[1] = cbrt((*box).Volume);
//      (*box).Len[2] = cbrt((*box).Volume);

      if( (*box).N[0] < 0 ){
	(*box).N[0] = 1;
	(*box).N[1] = 1;
	(*box).N[2] = 1;}

	(*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
	(*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
	(*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

  }
 */

  (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
  (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
  (*box).Area[5] = (*box).Len[0]*(*box).Len[1];

  (*gas).Fn = (*gas).n*(*box).Volume/ (*gas).N;

double thermal;
thermal =  sqrt( (*gas).kb*(*gas).T/(*gas).m );

    //(*gas).nualpha = 1.4;

//if( (*gas).model == "DSMC_HS"){
      //(*gas).visp = 0.5;
(*gas).T0 = 273.0;
    if((*gas).epsilon < 0.0){
      (*gas).epsilon = 119.8*(*gas).kb;//4.47478200e-21;
     }


if(fabs((*gas).sigma-1.0)<1e-14){
      //do nothing
      (*gas).crref = thermal;
}
else if ( fabs( (*gas).visp - 0.5 )< 1e-13 ){
	(*gas).sigma = 1.016*5.0*sqrt( (*gas).m*(*gas).kb*(*gas).T0/pi );
	(*gas).sigma = (*gas).sigma/( 16.0*(*gas).mu );
	(*gas).sigma = sqrt((*gas).sigma);
	(*gas).crref = thermal;
}
//if( (*gas).model == "DSMC_VHS"){
else{
        //(*gas).visp = 0.81;
	(*gas).sigma = 1.016*5.0*((*gas).nualpha+1.0)*((*gas).nualpha+2.0)*sqrt( (*gas).m*(*gas).kb*(*gas).T0/pi );
	(*gas).sigma = (*gas).sigma/( (4.0*(*gas).nualpha)*(5.0-2.0*(*gas).visp)*(7.0-2.0*(*gas).visp)*(*gas).mu );
	(*gas).sigma = sqrt((*gas).sigma);

	double numer = (*gas).m*pow(4.0,(*gas).visp-1)*5.0*( (*gas).nualpha+1.0 )*( (*gas).nualpha+2.0 )*sqrt( pi );
        double denom =  pi*pow((*gas).sigma,2.0)*(*gas).nualpha*( 5.0-2.0*(*gas).visp )*( 7.0-2.0*(*gas).visp )*2.0*tgamma( 5.0/2.0-(*gas).visp );
        //double gamma_t = tgamma( 5.0/2.0-(*gas).visp );
        (*gas).crref = (numer/denom)*pow( (*gas).kb*273.0/(*gas).m,(*gas).visp )/(*gas).mu;
        (*gas).crref = pow((*gas).crref, 1.0/((*gas).visp*2.0-1.0));
}

    if ( (*gas).model == "DSMC_VHS" ){
//    (*gas).crref = (*gas).m*pow(4.0,(*gas).visp-1)*5.0*( (*gas).nualpha+1.0 )*( (*gas).nualpha+2.0 )*sqrt( pi );

    }

  if(  (*gas).model == "MD" && (*box).problem == vacuum ){
	(*gas).sigma = 3.405e-10;
        double L;
	double eqdist = 1.0*pow(2.0, 1.0/6.0)*(*gas).sigma;
	(*gas).eqdist = eqdist;
//	double eqdist = pow(1.0/(*gas).n,1.0/3.0);

	int N = 30;
	int N13 = (*gas).N13;
	(*box).N_grid[0] = N13;
	(*box).N_grid[1] = N;
	(*box).N_grid[2] = N13;
	(*gas).n = 1.0/(eqdist*eqdist*eqdist);

	(*box).Len[0] = N13*eqdist;
	(*box).Len[1] = 4.0*N*eqdist;
	(*box).Len[2] = N13*eqdist;
        (*gas).Fn = 1.0;
        //(*gas).N = 5000.0*(*box).N[1];
/*
        L = pow((*gas).N/(*box).Len[1]/(*gas).n/0.5,1.0/2.0);
        (*box).Len[0] = L;
        (*box).Len[1] = (*box).Len[1];
        (*box).Len[2] = L;
*/
        (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
        (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
        (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
        (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
        (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
        (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

/*
	(*box).N_grid[1] = floor( pow((*gas).N/(L*L)*(*box).Len[1]*(*box).Len[1] ,1.0/3.0) );
	(*box).N_grid[0] = floor( (*box).N_grid[1]*L/(*box).Len[1] );
	(*box).N_grid[2] = (*box).N_grid[0];
*/
	(*gas).N = (*box).N_grid[0]*(*box).N_grid[1]*(*box).N_grid[2];
	//(*gas).n = (*gas).N/(*box).Volume;
	//(*gas).Fn = (*gas).n*(*box).Volume/ (*gas).N;
	printf("\n\n grid= %d %d %d\n", (*box).N_grid[0], (*box).N_grid[1], (*box).N_grid[2]);
	printf("N = %d\n",  (*box).N_grid[0]*(*box).N_grid[1]*(*box).N_grid[2]);
	printf("n = %e\n", (*gas).n);
	printf("sigma = %e\n",(*gas).sigma);
	printf("eqdist = %e\n",eqdist);
	printf("gas.T = %e\n",(*gas).T);
  }


else if ( (*box).problem == dcones){
////////////////		Double cone
	double Lx = (*box).Lx[0] + (*box).Lx[1] + (*box).Lx[2];
	int Nx = (*box).Nx[0] + (*box).Nx[1] + (*box).Nx[2];
	//int Nc = ceil( (1.0*(*gas).N)/(1.0*(*box).Ny*Nx) );
	//(*gas).N = Nc*(*box).Ny*Nx;
	//(*box).Nc = Nc;

        (*box).Volume = acos(-1.0)*(*box).Ly*(*box).Ly*(*box).Lx[0];
	double R1U = (*box).Ly, R2U = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1], R1L = 0.0, R2L = tan((*box).alpha[0])*(*box).Lx[1];
	(*box).Volume += volume_conical_frustum(R1U, R2U, (*box).Lx[1]) - volume_conical_frustum(R1L, R2L, (*box).Lx[1]);
	R1U = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1];
	R2U = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1]+tan((*box).alpha[1])*(*box).Lx[2];
	R1L = tan((*box).alpha[0])*(*box).Lx[1];
	R2L = tan((*box).alpha[0])*(*box).Lx[1]+tan((*box).alpha[1])*(*box).Lx[2];
	(*box).Volume += volume_conical_frustum(R1U, R2U,(*box).Lx[2])
		 	-volume_conical_frustum(R1L, R2L,(*box).Lx[2]);
	(*gas).Fn = (*gas).n*(*box).Volume/(*gas).N;
}
else if( (*box).problem == flatnose ){
	double Lx = (*box).L1[0]+(*box).L1[1];
	double Ly = (*box).L2[0]+(*box).L2[1];
	(*box).Volume = pi*Ly*Ly*Lx - pi*(*box).L2[0]*(*box).L2[0]*(*box).L1[1];
	(*gas).Fn = 0.085e14;
	(*gas).n = (*gas).Fn*(*gas).N/(*box).Volume;
	//(*gas).Fn = (*gas).n*(*box).Volume/(*gas).N;
}

else if ( (*box).problem == vacuum && (*gas).model != "MD"){

	(*gas).sigma = 3.405e-10;
        double L;
	double eqdist = 1.0*pow(2.0, 1.0/6.0)*(*gas).sigma;
	(*gas).eqdist = (*gas).sigma;//eqdist;//(*gas).sigma*( 1.068+0.3837*(*gas).kb*(*gas).T/(*gas).epsilon )/( 1.0+0.4293*(*gas).kb*(*gas).T/(*gas).epsilon );

	int N = 30;
	(*gas).n = 1.0/(eqdist*eqdist*eqdist);


	(*box).Len[1] = 4.0*N*eqdist;
	(*gas).Fn = (*gas).n*(*box).Len[0]*(*box).Len[2]*(*box).Len[1]/4.0/(*gas).N;

	(*gas).Nv1 = 0.0;//floor( 0.01*(*gas).n*(*box).Len[0]*(*box).Len[2]*(*box).Len[1]*3.0/4.0/(*gas).Fn );

	(*gas).N = (*gas).Nv1 + (*gas).N;
        (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
        (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
        (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
        (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
        (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
        (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

	printf("N = %ld\n",  (*gas).N);
	printf("n = %e\n", (*gas).n);
	printf("sigma = %e\n",(*gas).sigma);
	printf("eqdist = %e\n",(*gas).eqdist);
	printf("gas.T = %e\n",(*gas).T);
}
if(  (*gas).model == "MD" && (*box).problem == inverted ){
	(*gas).sigma = 3.405e-10;
	double d;
	int N, N13;
	double eqdist = pow(2.0, 1.0/6.0)*(*gas).sigma;
	//(*gas).eqdist = eqdist;

	printf("\n\n   :::   Initial setting   :::\n");
	printf("   nc %lf, nh %lf, nv1 %lf, nv2 %lf, nv3 %lf\n\n\n", (*gas).nc, (*gas).nh, (*gas).nv1, (*gas).nv2, (*gas).nv3);


	(*gas).nc = (*gas).nc/pow(eqdist, 3.0);
	(*gas).nh = (*gas).nh/pow(eqdist, 3.0);
	(*gas).nv1 = (*gas).nv1/pow(eqdist, 3.0);
	(*gas).nv2 = (*gas).nv2/pow(eqdist, 3.0);
	(*gas).nv3 = (*gas).nv3/pow(eqdist, 3.0);

	// cold droplet:
	d = pow(1.0/(*gas).nc, 1.0/3.0);
	(*gas).eqdist = d;
	N = floor((*gas).LcN*eqdist/d);
	(*gas).LcN = N;
	(*gas).Lc = d*N;
	(*gas).Lh = (*gas).Lc;

	(*box).Len[0] = N13*d;
	(*box).Len[2] = N13*d;

	N13 = (*gas).N13;

	(*gas).n = 0.5*((*gas).nc+(*gas).nh);

	d = pow(1.0/(*gas).nv2, 1.0/3.0);
	(*gas).Lv = (*gas).LvN*d;
	(*gas).Lv1 = 0.0;
	(*gas).Lv3 = 0.0;

	(*gas).Lc_center = (*gas).Lv1+0.5*(*gas).Lc;
	(*gas).Lh_center = (*gas).Lv+0.5*(*gas).Lh+1.0*(*gas).Lc+(*gas).Lv1;
	(*gas).Lc_ratio = 1.0*(*gas).Lc;
	(*gas).Lh_ratio = 1.0*(*gas).Lh;
	(*gas).n_ratio = 1;
	(*gas).Fn = 1.0;

	(*box).N_grid[0] = N13;
	(*box).N_grid[1] = N;
	(*box).N_grid[2] = N13;



	(*gas).n = 0.5*((*gas).nc+(*gas).nh);

	(*gas).Nc = (*gas).N13*(*gas).N13*N;
	(*gas).Nh = (*gas).N13*(*gas).N13*N;

	//(*gas).nc = (*gas).Nc*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lc );

	(*gas).Nv1 = 0;//floor( (*gas).nv1*(*box).Len[0]*(*box).Len[2]*0.5*(*gas).Lc/(*gas).Fn );
	(*gas).nv1 = 0.0;//(*gas).Nv1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*0.5*(*gas).Lc );

	(*gas).Nv2 = floor( (*gas).nv2*(*box).Len[0]*(*box).Len[2]*(*gas).Lv/(*gas).Fn );
	(*gas).nv2 = (*gas).Nv2*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv );

	(*gas).Nv3 = 0;//floor( (*gas).nv3*(*box).Len[0]*(*box).Len[2]*0.5*(*gas).Lh/(*gas).Fn );
	(*gas).nv3 = 0.0;//(*gas).Nv3*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*0.5*(*gas).Lh );

	(*gas).N = (*gas).Nc + (*gas).Nh + (*gas).Nv1 + (*gas).Nv2 + (*gas).Nv3;
	(*box).Len[1] = (*gas).Lv1+(*gas).Lh+(*gas).Lv+(*gas).Lc+(*gas).Lv3;

        (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
        (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
        (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
        (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
        (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
        (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

	(*gas).T = 0.5*((*gas).Tc+(*gas).Th);
	printf("N = %ld\n",  (*gas).N);
	printf("Nc = %d\n",  (*gas).Nc);
	printf("Nh = %d\n",  (*gas).Nh);
	printf("nc = %e\n", (*gas).nc);
	printf("nh = %e\n", (*gas).nh);
	printf("sigma = %e\n",(*gas).sigma);
	printf("eqdist = %e\n",(*gas).eqdist);
	printf("gas.Tc = %e\n",(*gas).Tc);
	printf("gas.Th = %e\n",(*gas).Th);

  }
else if(  (*box).problem == zoomed_inverted ){
  // in the input, we should give n, Lvn, length of other dimensions, N, T
  (*gas).sigma = 3.405e-10;
  (*gas).eqdist = (*gas).sigma;
  (*gas).n = (*gas).n/pow((*gas).sigma, 3.0);
  (*gas).nvc = (*gas).nvc/pow((*gas).sigma, 3.0);
  (*gas).nvh = (*gas).nvh/pow((*gas).sigma, 3.0);

  (*gas).Lv = (*gas).LvN*(*gas).eqdist;
  (*gas).Fn = (*gas).n*(*box).Len[0]*(*box).Len[2]*(*gas).Lv/( 1.0*(*gas).N );

  (*box).Len[1] = (*gas).Lv;

  (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
  (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
  (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
  (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
  (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
  (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

  printf("N = %ld\n",  (*gas).N);

  printf("sigma = %e\n",(*gas).sigma);
  printf("eqdist = %e\n",(*gas).eqdist);
  printf("vapour T: %e\n",(*gas).T);
  printf("n=%e  \n", (*gas).n);

  printf("\n+++ sigc = %e\n",(*gas).sigc);
  printf("+++ alphac = %e\n\n",(*gas).alphac);
  printf("\n+++ sigh = %e\n",(*gas).sigh);
  printf("+++ alphah = %e\n\n",(*gas).alphah);

}
else if ( (*box).problem == inverted && (*gas).model != "MD"){

	(*gas).sigma = 3.405e-10;
	double eqdist = (*gas).sigma;//1.0*pow(2.0, 1.0/6.0)*(*gas).sigma;
	printf("eqdist    =     %e\n",eqdist);
	(*gas).eqdist = (*gas).sigma;

	printf("\n\n   :::   Initial setting   :::\n");
	printf("   nc %lf, nh %lf, nv1 %lf, nv2 %lf, nv3 %lf\n\n\n", (*gas).nc, (*gas).nh, (*gas).nv1, (*gas).nv2, (*gas).nv3);

	(*gas).nc = (*gas).nc/pow(eqdist, 3.0);
	(*gas).nh = (*gas).nh/pow(eqdist, 3.0);
	(*gas).nv1 = (*gas).nv1/pow(eqdist, 3.0);
	(*gas).nv2 = (*gas).nv2/pow(eqdist, 3.0);
	(*gas).nv3 = (*gas).nv3/pow(eqdist, 3.0);

	(*gas).n = 0.5*((*gas).nc+(*gas).nh);
	(*gas).Lc = (*gas).LcN*eqdist;
	(*gas).Lh = (*gas).LhN*eqdist;
	(*gas).Lv = (*gas).LvN*eqdist;
	(*gas).Lv1 = (*gas).LvN1*eqdist;
	(*gas).Lv3 = (*gas).LvN3*eqdist;

	(*gas).Lc_center = (*gas).Lv1+0.5*(*gas).Lc;
	(*gas).Lh_center = (*gas).Lv+0.5*(*gas).Lh+1.0*(*gas).Lc+(*gas).Lv1;
	(*gas).Lc_ratio = (*gas).Lc;//0.05*(*gas).Lc;
	(*gas).Lh_ratio = (*gas).Lh;//0.05*(*gas).Lh;

	if( (*gas).n_ratio < 0 )
		(*gas).n_ratio = floor( (*gas).nh/(*gas).nv2 );

	(*gas).Fn = (*gas).nv2*(*box).Len[0]*(*box).Len[2]*(*gas).Lv/( 1.0*(*gas).Nv2 );

	int sg1, sg3;
	if((*gas).LvN1>0)
		sg1 = 1;
	else
		sg1 = 0;

	if((*gas).LvN3>0)
		sg3 = 1;
	else
		sg3 = 0;


	(*gas).Nc1 = 0;//floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) )/(1.0*(*gas).Fn) )*sg1;
	(*gas).Nc2 = floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( (*gas).Lc_ratio  )/(1.0*(*gas).n_ratio*(*gas).Fn) );
	(*gas).Nc3 = 0;//floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) )/(1.0*(*gas).Fn) );
	(*gas).Nc = (*gas).Nc1+(*gas).Nc2+(*gas).Nc3;

	//(*gas).nc = (*gas).Nc1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) ) );
	(*gas).nc = (*gas).Nc2*(*gas).n_ratio*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( (*gas).Lc_ratio ) );


	(*gas).Nh1 = 0;//floor( (*gas).nh*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lh-(*gas).Lh_ratio ) )/(1.0*(*gas).Fn) );
	(*gas).Nh2 = floor( (*gas).nh*(*box).Len[0]*(*box).Len[2]*( (*gas).Lh_ratio )/(1.0*(*gas).n_ratio*(*gas).Fn) );
	(*gas).Nh3 = 0;//floor( (*gas).nh*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lh-(*gas).Lh_ratio ) )/(1.0*(*gas).Fn) )*sg3;
	(*gas).Nh = (*gas).Nh1+(*gas).Nh2+(*gas).Nh3;

	//(*gas).nh = (*gas).Nh1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( (*gas).Lh-(*gas).Lh_ratio ) );
	(*gas).nh = (*gas).Nh2*(*gas).n_ratio*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( (*gas).Lh_ratio ) );



	(*gas).Nv1 = floor( (*gas).nv1*(*box).Len[0]*(*box).Len[2]*(*gas).Lv1/(*gas).Fn );
	if((*gas).Lv1 > 0.0)
		(*gas).nv1 = (*gas).Nv1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv1 );

	(*gas).Nv3 = floor( (*gas).nv3*(*box).Len[0]*(*box).Len[2]*(*gas).Lv3/(*gas).Fn );
	if((*gas).Lv3 > 0.0)
		(*gas).nv3 = (*gas).Nv3*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv3 );

	(*gas).N = (*gas).Nc1 + (*gas).Nc2 + (*gas).Nc3 + (*gas).Nh1 + (*gas).Nh2+ (*gas).Nh3 + (*gas).Nv1 + (*gas).Nv2 + (*gas).Nv3;
	(*box).Len[1] = (*gas).Lv+(*gas).Lh+(*gas).Lc+(*gas).Lv1+(*gas).Lv3;

        (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
        (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
        (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
        (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
        (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
        (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

	(*gas).T = 0.5*((*gas).Tc+(*gas).Th);

  if((*box).x2_prob[0]<0.0){
    (*box).x2_prob[0] = (*gas).Lc-3.0*(*gas).sigma;
    (*box).x2_prob[1] = (*gas).Lc+6.0*(*gas).sigma;
    (*box).x2_prob[2] = (*box).Len[1]-(*gas).Lh - 6.0*(*gas).sigma;
    (*box).x2_prob[3] = (*box).Len[1]-(*gas).Lh + 3.0*(*gas).sigma;
  }
  else{
    for(int i=0; i<4; i++){
      (*box).x2_prob[i] =(*box).x2_prob[i]*(*gas).sigma;
    }
  }

	printf("N = %ld\n",  (*gas).N);
	printf("Nc1:3 = %d  %d  %d\n",  (*gas).Nc1, (*gas).Nc2, (*gas).Nc3);
	printf("Nh1:3 = %d  %d  %d\n",  (*gas).Nh1, (*gas).Nh2, (*gas).Nh3);
	printf("Nv1 = %d\n",  (*gas).Nv1);
	printf("Nv2 = %d\n",  (*gas).Nv2);
	printf("Nv3 = %d\n",  (*gas).Nv3);

	printf("sigma = %e\n",(*gas).sigma);
	printf("eqdist = %e\n",(*gas).eqdist);
	printf("Liquid T = %e, %e\n",(*gas).Tc, (*gas).Th);
	printf("vapour T: %e, %e, %e\n",(*gas).Tv1, (*gas).Tv2, (*gas).Tv3);
	printf("nv1=%e  nv2=%e  nv3=%e\n", (*gas).nv1, (*gas).nv2, (*gas).nv3);
	printf("nc=%e  nh=%e\n", (*gas).nc, (*gas).nh);
  printf("\n ++ x2_prob = %e %e %e %e\n",(*box).x2_prob[0]/(*gas).sigma,(*box).x2_prob[1]/(*gas).sigma,(*box).x2_prob[2]/(*gas).sigma,(*box).x2_prob[3]/(*gas).sigma);
}
else if ( (*box).problem == evaporation || (*box).problem == wall){
	//double sc = 1e-10;
	//double gcm3tom3 = 1e3/((*gas).m);
	(*gas).sigma = 3.405e-10;
	(*gas).eqdist = (*gas).sigma;

	printf("\n\n   :::   Initial setting   :::\n");
	printf("   nc %lf, nv1 %lf, nv2 %lf\n\n\n", (*gas).nc, (*gas).nv1, (*gas).nv2);

	(*gas).nc = (*gas).nc/pow((*gas).sigma, 3.0);
	(*gas).nh = 0.0;
	(*gas).nv1 = (*gas).nv1/pow((*gas).sigma, 3.0);
	(*gas).nv2 = (*gas).nv2/pow((*gas).sigma, 3.0);
	(*gas).nv3 = 0.0;
	(*gas).LvN3 = 0;

	(*gas).n = (*gas).nc;
	(*gas).Lc = (*gas).LcN*(*gas).sigma;
	(*gas).Lh = 0.0;
	(*gas).Lv = (*gas).LvN*(*gas).sigma;
	(*gas).Lv1 = (*gas).LvN1*(*gas).sigma;
	(*gas).Lv3 = 0.0;

	(*gas).Lc_center = (*gas).Lv1+0.5*(*gas).Lc;
	(*gas).Lh_center = 0.0;

	// for 100K
	//(*gas).Lc_ratio = 0.5*(*gas).Lc;
	// for 80 K
	(*gas).Lc_ratio = (*gas).Lc;

	(*gas).Lh_ratio = 0.0;

	if( (*gas).n_ratio < 0 )
		(*gas).n_ratio = floor( (*gas).nc/(*gas).nv2 );

	//(*gas).Fn = (*gas).nv2*(*box).Len[0]*(*box).Len[2]*(*gas).Lv/( 1.0*(*gas).Nv2 );
  // here, we fixed the FN
  (*gas).Fn = 0.1e+14;
	int sg1;
	if((*gas).LvN1>0)
		sg1 = 1;
	else
		sg1 = 0;


	(*gas).Nc1 = 0;//floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) )/(1.0*(*gas).Fn) )*sg1;
	(*gas).Nc2 = floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( (*gas).Lc_ratio + 0.5*( (*gas).Lc-(*gas).Lc_ratio )*(1-sg1) )/(1.0*(*gas).n_ratio*(*gas).Fn) );
	(*gas).Nc3 = 0;//floor( (*gas).nc*(*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) )/(1.0*(*gas).Fn) );
	(*gas).Nc = (*gas).Nc1+(*gas).Nc2+(*gas).Nc3;

	//(*gas).nc = (*gas).Nc1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( 0.5*( (*gas).Lc-(*gas).Lc_ratio ) ) );
	(*gas).nc = (*gas).Nc2*(*gas).n_ratio*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*( (*gas).Lc_ratio ) );

	(*gas).Nh1 = 0;
	(*gas).Nh2 = 0;
	(*gas).Nh3 = 0;
	(*gas).Nh = 0;

	(*gas).nh = 0.0;
	(*gas).Nv3 = 0;

	(*gas).Nv1 = floor( (*gas).nv1*(*box).Len[0]*(*box).Len[2]*(*gas).Lv1/(*gas).Fn );
  (*gas).Nv2 = floor( (*gas).nv2*(*box).Len[0]*(*box).Len[2]*(*gas).Lv/(*gas).Fn );
	if((*gas).Lv1 > 0.0)
		(*gas).nv1 = (*gas).Nv1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv1 );

	(*gas).N = (*gas).Nc1 + (*gas).Nc2 + (*gas).Nc3 + (*gas).Nv1 + (*gas).Nv2;
	(*box).Len[1] = (*gas).Lv+(*gas).Lc+(*gas).Lv1;

        (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
        (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
        (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
        (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
        (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
        (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
        (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
        (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

	(*gas).T = (*gas).Tc;

  if((*box).x2_prob[0]<0.0){

    //(*box).x2_prob[2] = 30*1e-10;
		//(*box).x2_prob[0] = 70*1e-10;
    //(*box).x2_prob[3] = 50*1e-10;
		//(*box).x2_prob[1] = 90*1e-10;

    (*box).x2_prob[0] = (*gas).Lv1 + (*gas).Lc -6.0*(*gas).sigma;
    (*box).x2_prob[1] = (*gas).Lv1 + (*gas).Lc +3.0*(*gas).sigma;
    (*box).x2_prob[2] = (*gas).Lv1 - 3.0*(*gas).sigma;
    (*box).x2_prob[3] = (*gas).Lv1 + 6.0*(*gas).sigma;

    printf("\n ++ x2_prob = %lf %lf %lf %lf\n",(*box).x2_prob[0]/(*gas).sigma,(*box).x2_prob[1]/(*gas).sigma,(*box).x2_prob[2]/(*gas).sigma,(*box).x2_prob[3]/(*gas).sigma);

    printf("nv1=%e  nv2=%e  nc=%e \n", (*gas).nv1, (*gas).nv2, (*gas).nc);

  }
}
else if ( (*box).problem == injection){
	//double sc = 1e-10;
	//double gcm3tom3 = 1e3/((*gas).m);
	(*gas).sigma = 3.405e-10;
	(*gas).eqdist = (*gas).sigma;

	printf("\n\n   :::   Initial setting   :::\n");
	printf("   nc,Tc (state of chamber) (%e,%e)\n   nh,Th (state of injected liquid) (%e,%e)\n\n", (*gas).nc, (*gas).Tc, (*gas).nh, (*gas).Th);

	(*gas).nc = (*gas).nc/pow((*gas).sigma, 3.0);
	(*gas).nh =  (*gas).nh/pow((*gas).sigma, 3.0);
  (*gas).N = (*gas).N*(*box).N[0]*(*box).N[1];
  (*box).Len[0] = (*box).Len[0]*(*gas).sigma;
  (*box).Len[1] = (*box).Len[1]*(*gas).sigma;
  (*gas).T = (*gas).Tc;
  (*gas).n = (*gas).nc;
  (*gas).Fn = (*gas).n*(*box).Len[0]*(*box).Len[1]*(*box).Len[2]/( 1.0*(*gas).N );

  (*box).drop_dia = (*box).drop_dia*(*gas).sigma;
  (*gas).Nh = floor(((*gas).nh-(*gas).nc)*pi*(*box).drop_dia*(*box).drop_dia/4.0/(*gas).Fn);
  (*gas).N = (*gas).N + (*gas).Nh;
  printf("Nh=%d   N=%ld\n", (*gas).Nh, (*gas).N - (*gas).Nh);

  (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
  (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
  (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
  (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
  (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
  (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
  (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
  (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];
  }
  else if ( (*box).problem == shock){
  	//double sc = 1e-10;
  	//double gcm3tom3 = 1e3/((*gas).m);
  	(*gas).sigma = 3.405e-10;
  	(*gas).eqdist = (*gas).sigma;

  	printf("\n\n   :::   Initial setting   :::\n");
  	printf("    nv1 %lf, nv3 %lf\n\n\n",  (*gas).nv1, (*gas).nv3);

  	(*gas).nv1 = (*gas).nv1/(*gas).m*1e-6;
  	(*gas).nv3 = (*gas).nv3/(*gas).m*1e-6;
  	(*gas).Lv1 = (*gas).Lv1*1.0;
  	(*gas).Lv3 = (*gas).Lv3*1.0;
    (*gas).T = 0.5*((*gas).Tv1+(*gas).Tv3);
    (*gas).T0 = 273.0;

    (*gas).crref = sqrt( (*gas).kb*(*gas).T/(*gas).m );


  	(*gas).Fn =(*gas).nv1*(*box).Len[0]*(*box).Len[2]*(*gas).Lv1/( 1.0*(*gas).Nv1 );


  	(*gas).Nv1 = floor( (*gas).nv1*(*box).Len[0]*(*box).Len[2]*(*gas).Lv1/(*gas).Fn );
  	(*gas).nv1 = (*gas).Nv1*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv1 );

    (*gas).Nv3 = floor( (*gas).nv3*(*box).Len[0]*(*box).Len[2]*(*gas).Lv3/(*gas).Fn );
  	(*gas).nv3 = (*gas).Nv3*(*gas).Fn/( (*box).Len[0]*(*box).Len[2]*(*gas).Lv3 );


    if((*gas).model=="SPH" || (*gas).model=="Hybrid"){
      (*gas).Fn = 0.005*(*gas).Fn;
    }
  	(*gas).N =  (*gas).Nv1 + (*gas).Nv3;
  	(*box).Len[1] = (*gas).Lv1 + (*gas).Lv3;

          (*box).Volume = (*box).Len[0]*(*box).Len[1]*(*box).Len[2];
          (*box).Area[0] = (*box).Len[1]*(*box).Len[2];
          (*box).Area[1] = (*box).Len[1]*(*box).Len[2];
          (*box).Area[2] = (*box).Len[0]*(*box).Len[2];
          (*box).Area[3] = (*box).Len[0]*(*box).Len[2];
          (*box).Area[4] = (*box).Len[0]*(*box).Len[1];
          (*box).Area[5] = (*box).Len[0]*(*box).Len[1];
          (*box).delta_dim[0] = (*box).Len[0]/(*box).N[0];
          (*box).delta_dim[1] = (*box).Len[1]/(*box).N[1];
          (*box).delta_dim[2] = (*box).Len[2]/(*box).N[2];

  (*gas).n = 0.5*( (*gas).nv1+(*gas).nv3 );
	printf("N = %ld\n",  (*gas).N);
	printf("Nv1,3 = %d,  %d\n",  (*gas).Nv1, (*gas).Nv3);

	printf("sigma = %e\n",(*gas).sigma);
	printf("T = %e \n",(*gas).T );
	printf("nv1=%e  nv2=%e \n", (*gas).nv1, (*gas).nv3);



}
else{
	(*gas).eqdist = pow(2.0, 1.0/6.0)*(*gas).sigma;
}


//double rho = (*gas).n*(*gas).m;
//double b = 2*pi*pow((*gas).sigma,3)/(3.0);

(*gas).b = 2*pi*pow((*gas).sigma,3)/(3.0);
(*gas).mu = 5.0/( 16.0*(*gas).sigma*(*gas).sigma )*sqrt( (*gas).m*(*gas).kb*(*gas).T0/acos(-1.0) );
(*gas).k = 15.0*(*gas).kb*(*gas).mu/(4.0*(*gas).m);

  double Y = increase_collision_rate((*gas).n, gas);//(1-(*gas).b*rho/8.0)/ (pow((1.0 - (*gas).b*rho/4.0),3));

  //Y = 1+ 0.05556782*(*gas).n*(*gas).b+0.01394451*(*gas).n*(*gas).n*(*gas).b*(*gas).b- 0.0013396*(*gas).n*(*gas).n*(*gas).n*(*gas).b*(*gas).b*(*gas).b;
  //Y = Y/( 1.0 -  0.56943218*(*gas).n*(*gas).b + 0.08289011*(*gas).n*(*gas).n*(*gas).b*(*gas).b);

      printf("Y = %e \n",Y);

  //(*gas).lambda = 1.0/( Y* sqrt(2.0)*pi*(*gas).sigma*(*gas).sigma*(*gas).n);
  double lambda;
 lambda = 4.0*(*gas).nualpha*(5.0-2.0*(*gas).visp)*(7.0-2.0*(*gas).visp);
 lambda = lambda/( 5.0*((*gas).nualpha+1.0)*((*gas).nualpha+2.0) );
 lambda = lambda*sqrt( (*gas).m/(2.0*pi*(*gas).kb*(*gas).T) )*(*gas).mu/( (*gas).m*(*gas).n );
 (*gas).lambda = lambda;
 if((*box).problem == dcones){
	(*gas).Kn = (*gas).lambda/(*box).Lx[1];
 }
else if((*box).problem == flatnose ){
	(*gas).Kn = (*gas).lambda/(*box).L2[0];
}
else{
  	(*gas).Kn = (*gas).lambda/(*box).Len[1];
 }
  if ( (*gas).delta_t < 0.0 ){
    //if( (*gas).model == "CBA_HS" )
      //(*gas).delta_t = 0.04*( (*gas).lambda  ) /  (  Y* max(  sqrt( (*gas).kb*(*box).T_wall/(*gas).m ) , (*box).U_wall ) );
      //(*gas).delta_t = 0.04*( (*gas).lambda  ) /  (  Y* (  sqrt( (*gas).kb*(*box).T_wall/(*gas).m )  ) );
	//(*gas).delta_t = 0.5*(  (*box).Len[1]/(*box).N[1] ) / max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) , max( (*box).U_wall_1, (*box).U_wall_2) );
    if( (*gas).model == "DSMC_VHS" || (*gas).model == "CBA_VHS" || (*gas).model == "ESMC"||(*gas).model == "CBA_HS" || (*gas).model == "MD" || (*gas).model == "FP_dense"){
         if( (*gas).frac_mft < 0.0 )
	         (*gas).delta_t = 0.5*min(  (*gas).lambda, (*box).Len[1]/(*box).N[1])/max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) ,max( (*box).U_wall_1, (*box).U_wall_2) );
         else
		(*gas).delta_t = (*gas).frac_mft*min(  (*gas).lambda, (*box).Len[1]/(*box).N[1])/max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) , max( (*box).U_wall_1, (*box).U_wall_2) );
    }
    else
      if( (*gas).frac_mft < 0.0 )
	  (*gas).delta_t = 0.5*(*box).Len[1]/(*box).N[1] / max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) , max( (*box).U_wall_1, (*box).U_wall_2) );
      else
        (*gas).delta_t = (*gas).frac_mft*(*box).Len[1]/(*box).N[1] / max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) , max( (*box).U_wall_1, (*box).U_wall_2) );
}
if((*box).problem == dcones)
	(*gas).delta_t = (*gas).frac_mft*min(  (*gas).lambda, (*box).Lx[1]/(*box).Nx[1])/max(  sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m ) , (*box).U_wall_7 );
else if((*box).problem == flatnose)
	(*gas).delta_t = (*gas).frac_mft*min(  (*gas).lambda, (*box).L1[0]/(*box).N1[0])/sqrt( (*gas).kb*(*box).T_wall_1/(*gas).m );
//   (*gas).delta_t =  0.01*((*box).Len[1] / (*box).N[1] ) /sqrt((*gas).kb*(*gas).T/(*gas).m);


    if ( (*gas).phi0 < 0.0 )
      (*gas).phi0 = 0.0;

      if( (*gas).mu_corr < 0.0 )
	  (*gas).mu_corr = 3.0*(*gas).n*(*gas).b/4.0;


    //(*gas).crref = pow( (*gas).crref, 1.0/ (2.0*(*gas).visp-1.0) );

    (*gas).avgeraging_time_till_now = 0.0;


// old values: working
    (*gas).ast = -1.64835851e-28;
    (*gas).bst = 6.91304716e+09;

// values for testing with Marcel
//  (*gas).ast = -0.4E-28;
//  (*gas).bst = 6.91304716e+09;

  //  (*gas).ast = -1.0420475583159304e-27;
    //(*gas).bst = 11685753098.056;

 // (*gas).ast = -4.833784252e-28;
 // (*gas).bst = 9919700000.0;
//(*gas).ast = -1.0420541947931818e-27;
//(*gas).bst = 11685768576.264046;
// new values: minimizing pressure
//     (*gas).ast = -1.0252953129075e-27;
//    (*gas).bst = 11682533132.705481;
// new values: matching Sutherland
//    (*gas).ast =-4.7623362e-28;
//    (*gas).bst = 9.9195611e+9;

// old, which was sort of working
//    (*gas).phi = 2.23770358e-21;
// new
   if( (*gas).phi<0.0 ){
      (*gas).phi = 4.0*(*gas).epsilon;//2.2711776e-21
    }
   else{
     (*gas).phi = (*gas).phi*(*gas).epsilon;
   }
    if( (*box).ghost < 0.0 )
	(*box).ghost = 0.0;
    else
	(*box).ghost = (*box).ghost*(*gas).sigma;

  if( (*box).thermc < 0.0 )
  (*box).thermc = 0.0;
  else
  (*box).thermc = (*box).thermc*(*gas).sigma;

  if( (*box).thermh < 0.0 )
  (*box).thermh = 0.0;
  else
  (*box).thermh = (*box).thermh*(*gas).sigma;

if((*box).problem == dcones){
	printf("Double Cone dimensions dim.: Lx = [%e, %e, %e], Ly = %e\n",(*box).Lx[0],(*box).Lx[1],(*box).Lx[2], (*box).Ly);
	printf("Double Cone: num. of cells: Nx[%d, %d, %d] and Ny=%d\n",(*box).Nx[0],(*box).Nx[1],(*box).Nx[2], (*box).Ny);
	printf("Double Cone angles: alpha = %lf, %lf\n",(*box).alpha[0],(*box).alpha[1]);
	printf("Volume = %e\n", (*box).Volume);
}
else if((*box).problem == flatnose){
	printf("Flat nose dimensions dim.: L1 = [%e, %e], L2 = [%e, %e]\n",(*box).L1[0],(*box).L1[1],(*box).L2[0], (*box).L2[1]);
	printf("Flat nose: num. of cells: N1 = [%d, %d] and N2=[%d, %d]\n",(*box).N1[0],(*box).N1[1],(*box).N2[0], (*box).N2[1]);
	printf("Flatnose: s1 = %lf and s2 = %lf", (*box).s1, (*box).s2);
	printf("Volume = %e\n", (*box).Volume);
}
else{
    printf("Box: dim. = [%e, %e, %e]\n",(*box).Len[0],(*box).Len[1],(*box).Len[2]);
    printf("Box: num. of cells = [%d, %d, %d]\n",(*box).N[0],(*box).N[1],(*box).N[2]);
}
    printf("Box: dim. of cells = [%e, %e, %e]\n",(*box).delta_dim[0],(*box).delta_dim[1],(*box).delta_dim[2]);
    printf("Gas: sigma = %e\n", (*gas).sigma);
    printf("Gas: number of computational particles = %ld\n", (*gas).N);
    printf("Gas: number density = %e\n",(*gas).n);
    printf("Gas: density = %e\n",(*gas).n*(*gas).m);
    printf("Gas: weight = %e\n",(*gas).w);
    printf("Gas: mass = %e\n",(*gas).m);
    printf("Gas: T = %e\n",(*gas).T);
    printf("Gas: kb = %e\n",(*gas).kb);
    printf("Gas: std of initial Vel = [%lf, %lf, %lf]\n",(*gas).U0[0],(*gas).U0[1],(*gas).U0[2]);
    printf("dt = %e\n",(*gas).delta_t);
    printf("number of steps = %d\n",(*box).num_steps);
    printf("number of realizations = %d\n",(*box).num_realizations);
    printf("Knudsen number = %e\n",(*gas).Kn);
    printf("Mean Free Path = %e\n",(*gas).lambda);
    printf("Mean Free Time = %e\n",(*gas).lambda/sqrt((*gas).kb*(*gas).T/(*gas).m));
//     printf("Collision Model = %s\n",(*gas).model);
    std::cout << "Collision Model = "<< (*gas).model << "\n";
    std::cout << "print on console every = "<< (*box).every << " steps \n";
    std::cout << "calculate pressure after "<< (*box).after << "\n";
    std::cout << "Reset Every "<< (*box).reset_every << "\n";
    std::cout << "Fn = "<< (*gas).Fn << "\n";
    std::cout << "mu = "<< (*gas).mu << "\n";
    std::cout << "factor = "<< (*gas).factor << "\n";
    printf("U_wall_1 = %lf\n", (*box).U_wall_1);
    printf("U_wall_2 = %lf\n", (*box).U_wall_2);
    printf("U_wall_3 = %lf\n", (*box).U_wall_3);
    printf("U_wall_4 = %lf\n", (*box).U_wall_4);
    printf("T_wall_1 = %lf\n", (*box).T_wall_1);
    printf("T_wall_2 = %lf\n", (*box).T_wall_2);
    printf("T_wall_3 = %lf\n", (*box).T_wall_3);
    printf("T_wall_4 = %lf\n", (*box).T_wall_4);

    if((*box).problem == closed_box)
      printf("problem is the simple closed box\n");
    else if ((*box).problem == couette_flow)
      printf("problem is the couette flow\n");
    printf("mu^corr = %e \n",(*gas).mu_corr);

   printf("long_range = %d \n",(*gas).long_range);

    printf("\n\n");

    printf("direction= %d %d %d\n", (*box).direction[0], (*box).direction[1], (*box).direction[2]);


printf("ghost/sigma = %lf\n", (*box).ghost/(*gas).sigma);
printf("thermc/sigma = %lf\n", (*box).thermc/(*gas).sigma);
printf("thermh/sigma = %lf\n", (*box).thermh/(*gas).sigma);

printf(" (*box).ghost = %e\n", (*box).ghost );

printf(" (*gas).ast = %e,   (*gas).bst = %e\n", (*gas).ast, (*gas).bst );



    (*gas).U[0] = (*gas).U0[0];
    (*gas).U[1] = (*gas).U0[1];
    (*gas).U[2] = (*gas).U0[2];

    (*box).step = 0;
    (*gas).t_post = 0.0;
    delete [] cstr;
    return;
}
