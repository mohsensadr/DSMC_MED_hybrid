#include "cpp_headers.h"

void dcones_vlasov_integral(double *U1, double *U2, struct GAS *gas, struct BOX *box, struct CELLS *cells, int *index){

      int i, num;

      double pi = acos(-1);



      int j,id;
      double r0 = 1.122*(*gas).sigma;
      double cut_off = 3.0*r0;


      int NNr = 10;
      int NNphi = 10;
      int NNtheta = 10;
      double *rr = (double *) malloc( NNr * sizeof(double) );
      double *phi = (double *) malloc( NNphi * sizeof(double) );
      double *theta = (double *) malloc( NNtheta * sizeof(double) );

      double *delt_r = (double *) malloc( NNr * sizeof(double) );
      for (i=0; i<NNr; ++i){
		rr[i] = r0 + ( (cut_off-r0)*(i) )/( 1.0* NNr);
		delt_r[i] = (cut_off-r0)/( 1.0* NNr);
      }
      for (i=0; i<NNphi; ++i)
		phi[i] = 0.0 + (pi*(i+0.5))/(1.0*NNphi);
      for (i=0; i<NNtheta; ++i)
		theta[i] = 0.0 + (2.0*pi*(i+0.5))/(1.0*NNtheta);



double ndummy = 0.0;


int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
int Ny = (*box).Ny;
double x1,x2,x3,X,Y;
int num2, indx, indy;
double dphi, dtheta;
int ir, iphi, itheta;
dphi = pi/(1.0*NNphi);
dtheta = 2.0*pi/(1.0*NNtheta);

for(num=0; num<Nx*Ny; num++)
	for(i=0; i<4; i++)
		cells[num].V[i] = 0.0;
/// cells[].V[4]
//	 _3_
//	|   |
//	0   1
//	|_2_|
//
///// cells[].dim
//
//	2--3
//	|  |
//	0--1
for(j=0; j<Ny; j++){
	for(i=0; i<Nx; i++){
		num = j*Nx+i;
		// find V[0]
		if(i==0){
			for(ir=0; ir<NNr; ir++){
				for(iphi=0; iphi<NNphi; iphi++){
					for(itheta=0; itheta<NNtheta; itheta++){
						x1 = cells[num].cell_center[0]-cells[num].dx/2.0+rr[ir]*sin(phi[iphi])*sin(theta[itheta]);
						x2 = (cells[num].dim[0]+cells[num].dim[2])/2.0+rr[ir]*cos(theta[itheta]);
						x3 = rr[ir]*sin(phi[iphi])*cos(theta[itheta]);
						X = x1;
						Y = sqrt(x2*x2+x3*x3);
						num2 = dcone_index(X, Y, cells, box, &indx, &indy);
						if(indx>=0 && indx<Nx && indy>=0 && indy<Ny)
							ndummy = cells[num2].n;
						else if( ( X>(*box).Lx[0] && X<(*box).Lx[0]+(*box).Lx[1] && Y<tan((*box).alpha[0])*(X-(*box).Lx[0]) ) || ( X>(*box).Lx[0]+(*box).Lx[1] && Y<(*box).Ly+tan((*box).alpha[1])*(X-(*box).Lx[0]-(*box).Lx[1]) ))
							ndummy = (*box).dcones_n;
						else
							ndummy = (*gas).n;
						cells[num].V[0] += exp(-(*gas).bst*rr[ir])/(4.0*pi)*rr[ir]*sin(phi[iphi])*delt_r[ir]*dphi*dtheta*ndummy;
					}
				}
			}
		}
		else{
			cells[num].V[0] = cells[j*Nx+i-1].V[1];
		}
		// find V[1]
		for(ir=0; ir<NNr; ir++){
			for(iphi=0; iphi<NNphi; iphi++){
				for(itheta=0; itheta<NNtheta; itheta++){
					x1 = cells[num].cell_center[0]+cells[num].dx/2.0+rr[ir]*sin(phi[iphi])*sin(theta[itheta]);
					x2 = (cells[num].dim[1]+cells[num].dim[3])/2.0+rr[ir]*cos(theta[itheta]);
					x3 = rr[ir]*sin(phi[iphi])*cos(theta[itheta]);
					X = x1;
					Y = sqrt(x2*x2+x3*x3);
					num2 = dcone_index(X, Y, cells, box, &indx, &indy);
					if(indx>=0 && indx<Nx && indy>=0 && indy<Ny)
						ndummy = cells[num2].n;
					else if( ( X>(*box).Lx[0] && X<(*box).Lx[0]+(*box).Lx[1] && Y<tan((*box).alpha[0])*(X-(*box).Lx[0]) ) || ( X>(*box).Lx[0]+(*box).Lx[1] && Y<(*box).Ly+tan((*box).alpha[1])*(X-(*box).Lx[0]-(*box).Lx[1]) ))
						ndummy = (*box).dcones_n;
					else
						ndummy = (*gas).n;
					cells[num].V[1] += exp(-(*gas).bst*rr[ir])/(4.0*pi)*rr[ir]*sin(phi[iphi])*delt_r[ir]*dphi*dtheta*ndummy;
				}
			}
		}
		// find V[2]
		if(j==0){
			for(ir=0; ir<NNr; ir++){
				for(iphi=0; iphi<NNphi; iphi++){
					for(itheta=0; itheta<NNtheta; itheta++){
						x1 = cells[num].cell_center[0]+rr[ir]*sin(phi[iphi])*sin(theta[itheta]);
						x2 = (cells[num].dim[0]+cells[num].dim[1])/2.0+rr[ir]*cos(theta[itheta]);
						x3 = rr[ir]*sin(phi[iphi])*cos(theta[itheta]);
						X = x1;
						Y = sqrt(x2*x2+x3*x3);
						num2 = dcone_index(X, Y, cells, box, &indx, &indy);
						if(indx>=0 && indx<Nx && indy>=0 && indy<Ny)
							ndummy = cells[num2].n;
						else if( ( X>(*box).Lx[0] && X<(*box).Lx[0]+(*box).Lx[1] && Y<tan((*box).alpha[0])*(X-(*box).Lx[0]) ) || ( X>(*box).Lx[0]+(*box).Lx[1] && Y<(*box).Ly+tan((*box).alpha[1])*(X-(*box).Lx[0]-(*box).Lx[1]) ))
							ndummy = (*box).dcones_n;
						else
							ndummy = (*gas).n;
						cells[num].V[2] += exp(-(*gas).bst*rr[ir])/(4.0*pi)*rr[ir]*sin(phi[iphi])*delt_r[ir]*dphi*dtheta*ndummy;
					}
				}
			}
		}
		else{
			cells[num].V[2] = cells[(j-1)*Nx+i].V[3];
		}
		// find V[3]
		for(ir=0; ir<NNr; ir++){
			for(iphi=0; iphi<NNphi; iphi++){
				for(itheta=0; itheta<NNtheta; itheta++){
					x1 = cells[num].cell_center[0]+rr[ir]*sin(phi[iphi])*sin(theta[itheta]);
					x2 = (cells[num].dim[2]+cells[num].dim[3])/2.0+rr[ir]*cos(theta[itheta]);
					x3 = rr[ir]*sin(phi[iphi])*cos(theta[itheta]);
					X = x1;
					Y = sqrt(x2*x2+x3*x3);
					num2 = dcone_index(X, Y, cells, box, &indx, &indy);
					if(indx>=0 && indx<Nx && indy>=0 && indy<Ny)
						ndummy = cells[num2].n;
					else if( ( X>(*box).Lx[0] && X<(*box).Lx[0]+(*box).Lx[1] && Y<tan((*box).alpha[0])*(X-(*box).Lx[0]) ) || ( X>(*box).Lx[0]+(*box).Lx[1] && Y<(*box).Ly+tan((*box).alpha[1])*(X-(*box).Lx[0]-(*box).Lx[1]) ))
						ndummy = (*box).dcones_n;
					else
						ndummy = (*gas).n;
					cells[num].V[3] += exp(-(*gas).bst*rr[ir])/(4.0*pi)*rr[ir]*sin(phi[iphi])*delt_r[ir]*dphi*dtheta*ndummy;
				}
			}
		}
	}
}
for(i=0; i<Nx; i++){
	for(j=0; j<Ny; j++){
		num = j*Nx+i;
		if(i<(*box).Nx[0]){
			cells[num].F1 = (cells[num].V[1]-cells[num].V[0])/cells[num].dx;
			cells[num].F2 = (cells[num].V[3]-cells[num].V[2])/cells[num].dy;
		}
		else if(i<(*box).Nx[0]+(*box).Nx[1]){
			cells[num].F1 = (cells[num].V[1]-cells[num].V[0])/cells[num].dx + (cells[num].V[2]-cells[num].V[3])*sin((*box).alpha[0])/cells[num].dy;
			cells[num].F2 = (cells[num].V[3]-cells[num].V[2])/cells[num].dy*cos((*box).alpha[0]);
		}
		else{
			cells[num].F1 = (cells[num].V[1]-cells[num].V[0])/cells[num].dx + (cells[num].V[2]-cells[num].V[3])*sin((*box).alpha[1])/cells[num].dy;
			cells[num].F2 = (cells[num].V[3]-cells[num].V[2])/cells[num].dy*cos((*box).alpha[1]);
		}
		cells[num].F1 = cells[num].F1/(*gas).m*(*gas).ast;
		cells[num].F2 = cells[num].F2/(*gas).m*(*gas).ast;
	}
}
      for(id=0; id<(*gas).N; id++){
	    num = index[id];
	    U1[id] = U1[id] + cells[num].F1*(*gas).delta_t;
	    U2[id] = U2[id] + cells[num].F2*(*gas).delta_t;
      }
}

void dcones_intersection(double xi, double yi, double x, double y, double *xs, double *ys, double *len, int *w_id, struct BOX *box){
// finding the first intersection with the walls wall
// input: xi, yi, x, y
// output: xs, ys, len, w_id

/*
// id of walls:
	       @
	      /|
	     / |
	  4 /  |3
	   /   |
	  /    @
	 @    /
	/    /
      5/    /2
      /    /
  6  /    @
 @--@	 /
7|	/1
 |     /
 @----@
   0
*/
// the particles is originally on a segment of the line y=ax+b from (xi,yi) to (x,y)
double a,b;
double xw,yw;
double xd,yd;
double dlen;
int i;
double eps = 1e-15;
int check;
a = x-xi;
b = y-yi;
*len = sqrt((xi-x)*(xi-x) + (yi-y)*(yi-y));
*w_id = -1;

  for(i=0; i<8; i++){
	check = 0;
	//check out the wall
	if(i==0){
	///////////////////  wall = 0
		yw = 0.0;
		if(fabs(b)>1e-15){
			xw = a*(yw-yi)/b+xi;
			if( (b>0.0 && yi<yw && yw<y) || (b<0.0 && y<yw && yw<yi) ){
				if(xw>eps && xw<(*box).Lx[0]){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(fabs(y-yw)<eps && x>0.0 && x<(*box).Lx[0]){
				xd = x;
				yd = yw+eps;
				check = 1;
			}
		}
	}
	else if(i==1){
	///////////////////  wall = 1
		if(fabs(a)>1e-15){
			xw = (-(*box).Lx[0]*tan((*box).alpha[0]) - yi +(b/a)*xi)/(b/a-tan((*box).alpha[0]));
			yw = tan((*box).alpha[0])*(xw-(*box).Lx[0]);
			if( (fabs(b)>1e-15 && b>0.0 && yi<yw && yw<y) || (fabs(b)>1e-15 && b<0.0 && y<yw && yw<yi) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15)){
				if(xw>(*box).Lx[0] && xw<(*box).Lx[0]+(*box).Lx[1] && yw>0.0 && yw<tan((*box).alpha[0])*(*box).Lx[1] ){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(x>(*box).Lx[0] && x<(*box).Lx[0]+(*box).Lx[1]){
				xw = x;
				yw = tan((*box).alpha[0])*(xw-(*box).Lx[0]);
				if(yw>0.0 && yw<tan((*box).alpha[0])*(*box).Lx[1]){
					if( (fabs(b)>1e-15 && b>0.0 && yi-yw<1e-15 && yw-y<1e-15) || (fabs(b)>1e-15 && b<0.0 && y-yw<1e-15 && yw-yi<1e-15) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15) ){
						xd = xw;
						yd = yw+eps;
						check = 1;
					}
				}
			}
		}
	}
	else if(i==2){
	///////////////////  wall = 2
		if(fabs(a)>1e-15){
			xw = ( (*box).Lx[1]*tan((*box).alpha[0])-((*box).Lx[0]+(*box).Lx[1])*tan((*box).alpha[1]) - yi +(b/a)*xi )/( b/a-tan((*box).alpha[1]) );
			yw = tan((*box).alpha[1])*(xw-(*box).Lx[0]-(*box).Lx[1])+(*box).Lx[1]*tan((*box).alpha[0]);
			if( (fabs(b)>1e-15 && b>0.0 && yi-yw<1e-15 && yw-y<1e-15) || (fabs(b)>1e-15 && b<0.0 && y-yw<1e-15 && yw-yi<1e-15) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15) ){
				if(xw>(*box).Lx[0]+(*box).Lx[1] && xw<(*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2] && yw>tan((*box).alpha[0])*(*box).Lx[1] && yw<tan((*box).alpha[0])*(*box).Lx[1]+tan((*box).alpha[1])*(*box).Lx[2]){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(x>(*box).Lx[0]+(*box).Lx[1] && x<(*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2]){
				xw = x;
				yw = tan((*box).alpha[1])*(xw-(*box).Lx[0]-(*box).Lx[1]) + (*box).Lx[1]*tan((*box).alpha[0]);
				if( (fabs(b)>1e-15 && b>0.0 && yi-yw<1e-15 && yw-y<1e-15) || (fabs(b)>1e-15 && b<0.0 && y-yw<1e-15 && yw-yi<1e-15) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15) ){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
	}
	else if(i==3){
	///////////////////  wall = 3
		xw = (*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2];
		if(fabs(a)>1e-15){
			if( (a>0.0 && xi-xw<1e-15 && xw-x<1e-15) || (a<0.0 && x-xw<1e-15 && xw-xi<1e-15)){
				yw = (b/a)*(xw-xi)+yi;
				if( (*box).Lx[1]*tan((*box).alpha[0])+(*box).Lx[2]*tan((*box).alpha[1])<yw && yw <(*box).Ly +(*box).Lx[1]*tan((*box).alpha[0])+(*box).Lx[2]*tan((*box).alpha[1]) ){
					xd = xw-eps;
					yd = yw;
					check = 1;
				}
			}
		}
		else{
			if(fabs(x-xw)<eps){

				xd = xw-eps;
				yd = y;
				check = 1;
			}
		}
	}
	else if(i==4){
	///////////////////  wall = 4
		if(fabs(a)>1e-15){
			xw = ((*box).Lx[1]*tan((*box).alpha[0])+(*box).Ly-((*box).Lx[0]+(*box).Lx[1])*tan((*box).alpha[1]) - yi +(b/a)*xi)/(b/a-tan((*box).alpha[1]));
			yw = tan((*box).alpha[1])*(xw-(*box).Lx[0]-(*box).Lx[1])+(*box).Lx[1]*tan((*box).alpha[0])+(*box).Ly;
			if( (fabs(b)>1e-15 && b>0.0 && yi<yw && yw<y) || (fabs(b)>1e-15 && b<0.0 && y<yw && yw<yi) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15)){
				if(xw>(*box).Lx[0]+(*box).Lx[1] && xw<(*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2]){
					xd = xw;
					yd = yw-eps;
					check = 1;
				}
			}
		}
		else{
			if(x>(*box).Lx[0]+(*box).Lx[1] && x<(*box).Lx[0]+(*box).Lx[1]+(*box).Lx[2]){
				xw = x;
				yw = tan((*box).alpha[1])*(xw-(*box).Lx[0]-(*box).Lx[1]) + (*box).Lx[1]*tan((*box).alpha[0])+(*box).Ly;
				if( (fabs(b)>1e-15 && b>0.0 && yi-yw<2e-15 && yw-y<2e-15) || (fabs(b)>1e-15 && b<0.0 && y-yw<2e-15 && yw-yi<2e-15)  || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15)){
					xd = xw;
					yd = yw-eps;
					check = 1;
				}
			}
		}
	}
	else if(i==5){
	///////////////////  wall = 5
		if(fabs(a)>1e-15){
			xw = ( (*box).Ly-(*box).Lx[0]*tan((*box).alpha[0]) - yi +(b/a)*xi)/(b/a-tan((*box).alpha[0]));
			yw = tan((*box).alpha[0])*(xw-(*box).Lx[0])+(*box).Ly;
			if( (fabs(b)>1e-15 && b>0.0 && yi<yw && yw<y) || (fabs(b)>1e-15 && b<0.0 && y<yw && yw<yi) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15)){
				if(xw>(*box).Lx[0] && xw<(*box).Lx[0]+(*box).Lx[1]  && yw>(*box).Ly && yw<(*box).Ly+tan((*box).alpha[0])*(*box).Lx[1]){
					xd = xw;
					yd = yw-eps;
					check = 1;
				}
			}
		}
		else{
			if(x>(*box).Lx[0] && x<(*box).Lx[0]+(*box).Lx[1]){
				xw = x;
				yw = tan((*box).alpha[0])*(xw-(*box).Lx[0])+(*box).Ly;
				if( (fabs(b)>1e-15 && b>0.0 && yi<yw && yw<y) || (fabs(b)>1e-15 && b<0.0 && y<yw && yw<yi) || (fabs(b)<1e-15 && fabs(yw-yi)<1e-15)){
					xd = xw;
					yd = yw-eps;
					check = 1;
				}
			}
		}
	}
	else if(i==6){
	///////////////////  wall = 6
		yw = (*box).Ly;
		if(fabs(b)>1e-15){
			xw = a*(yw-yi)/b+xi;
			if( (b>0.0 && yi<yw && yw<y) || (b<0.0 && y<yw && yw<yi) ){
				if(xw>0.0 && xw<(*box).Lx[0]){
					xd = xw;
					yd = yw-eps;
					check = 1;
				}
			}
		}
		else{
			if(fabs(y-yw)<eps && x>0.0 && x<(*box).Lx[0]){
				xd = x;
				yd = yw-eps;
				check = 1;
			}
		}
	}
	else if(i==7){
	///////////////////  wall = 7
		xw = 0.0;
		if(fabs(a)>1e-15){
			if( (a>0.0 && xi<xw && xw<x) || (a<0.0 && x<xw && xw<xi) ){
				yw = (b/a)*(xw-xi)+yi;
				if( 0.0<yw && yw-(*box).Ly<1e-15 ){
					xd = xw+eps;
					yd = yw;
					check = 1;
				}
			}
		}
		else{
			if(fabs(x-xw)<eps){
				xd = xw+eps;
				yd = y;
				check = 1;
			}
		}
	}
	if(check ==1){
		dlen = sqrt( (xi-xd)*(xi-xd)+(yi-yd)*(yi-yd) );
		if(dlen - *len < 2e-15){
			*xs = xd;
			*ys = yd;
			*len = dlen;
			*w_id = i;
		}
	}
  }
}
void flatnose_intersection(double xi, double yi, double x, double y, double *xs, double *ys, double *len, int *w_id, struct BOX *box){
// finding the first intersection with the walls wall
// input: xi, yi, x, y
// output: xs, ys, len, w_id

/*
// id of walls:
_____4____
|         |
|         3
5    __2__|
|    |
|    1
|_0__|
*/
// the particles is originally on a segment of the line y=ax+b from (xi,yi) to (x,y)
double a,b;
double xw,yw;
double xd,yd;
double dlen;
int i;
double eps = 1e-15;
int check;
a = x-xi;
b = y-yi;
*len = sqrt((xi-x)*(xi-x) + (yi-y)*(yi-y));
*w_id = -1;

  for(i=0; i<6; i++){
	check = 0;
	//check out the wall
	if(i==0){
	///////////////////  wall = 0
		yw = 0.0;
		if(fabs(b)>1e-15){
			xw = a*(yw-yi)/b+xi;
			if( (b>0.0 && yi<yw && yw<y) || (b<0.0 && y<yw && yw<yi) ){
				if(xw>eps && xw<(*box).L1[0]){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(fabs(y-yw)<eps && x>0.0 && x<(*box).L1[0]){
				xd = x;
				yd = yw+eps;
				check = 1;
			}
		}
	}
	else if(i==1){
	///////////////////  wall = 1
		xw = (*box).L1[0];
		if(fabs(a)>1e-15){
			if( (a>0.0 && xi-xw<1e-15 && xw-x<1e-15) || (a<0.0 && x-xw<1e-15 && xw-xi<1e-15)){
				yw = (b/a)*(xw-xi)+yi;
				if( yw>0.0 && yw < (*box).L2[0] ){
					xd = xw-eps;
					yd = yw;
					check = 1;
				}
			}
		}
		else{
			if(fabs(x-xw)<eps){

				xd = xw-eps;
				yd = y;
				check = 1;
			}
		}
	}
	else if(i==2){
	///////////////////  wall = 2
		yw = (*box).L2[0];
		if(fabs(b)>1e-15){
			xw = a*(yw-yi)/b+xi;
			if( (b>0.0 && yi<yw && yw<y) || (b<0.0 && y<yw && yw<yi) ){
				if(xw-(*box).L1[0]>eps && xw<(*box).L1[0]+(*box).L1[1]){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(fabs(y-yw)<eps && x>(*box).L1[0] && x<(*box).L1[0]+(*box).L1[1]){
				xd = x;
				yd = yw+eps;
				check = 1;
			}
		}
	}
	else if(i==3){
	///////////////////  wall = 3
		xw = (*box).L1[0]+(*box).L1[1];
		if(fabs(a)>1e-15){
			if( (a>0.0 && xi-xw<1e-15 && xw-x<1e-15) || (a<0.0 && x-xw<1e-15 && xw-xi<1e-15)){
				yw = (b/a)*(xw-xi)+yi;
				if( yw>(*box).L2[0] && yw < (*box).L2[0]+(*box).L2[1] ){
					xd = xw-eps;
					yd = yw;
					check = 1;
				}
			}
		}
		else{
			if(fabs(x-xw)<eps){

				xd = xw-eps;
				yd = y;
				check = 1;
			}
		}
	}
	else if(i==4){
	///////////////////  wall = 4
		yw = (*box).L2[0]+(*box).L2[1];
		if(fabs(b)>1e-15){
			xw = a*(yw-yi)/b+xi;
			if( (b>0.0 && yi<yw && yw<y) || (b<0.0 && y<yw && yw<yi) ){
				if(xw>eps && xw<(*box).L1[0]+(*box).L1[1]){
					xd = xw;
					yd = yw+eps;
					check = 1;
				}
			}
		}
		else{
			if(fabs(y-yw)<eps && x>0.0 && x<(*box).L1[0]+(*box).L1[1]){
				xd = x;
				yd = yw+eps;
				check = 1;
			}
		}
	}
	else if(i==5){
	///////////////////  wall = 5
		xw = 0.0;
		if(fabs(a)>1e-15){
			if( (a>0.0 && xi-xw<1e-15 && xw-x<1e-15) || (a<0.0 && x-xw<1e-15 && xw-xi<1e-15)){
				yw = (b/a)*(xw-xi)+yi;
				if( yw>0.0 && yw < (*box).L2[0]+(*box).L2[1] ){
					xd = xw-eps;
					yd = yw;
					check = 1;
				}
			}
		}
		else{
			if(fabs(x-xw)<eps){

				xd = xw-eps;
				yd = y;
				check = 1;
			}
		}
	}
	if(check ==1){
		dlen = sqrt( (xi-xd)*(xi-xd)+(yi-yd)*(yi-yd) );
		if(dlen - *len < 2e-15){
			*xs = xd;
			*ys = yd;
			*len = dlen;
			*w_id = i;
		}
	}
  }
}
void flatnose_wall_collision(double *U1,double *U2,double *U3,double *x1,double *x2, struct BOX *box, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, int *w_id, double xs, double ys, double len, double *dt, struct GAS *gas){
// Update velocity and position after collision with the corresponding wall
/*
// id of walls:
_____4____
|         |
|         3
5    __2__|
|    |
|    1
|_0__|
*/
//  int i;
if(*w_id == -1){
	*dt = 0.0;
}
else{
double U10, U20, U30, U2D0;
double dt_rm;
double DX, DY, DZ, UI,VI,WI;
double eps = 1e-15;
  std::random_device rdnew;
  std::mt19937 gennew(rdnew());

  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);

U10 = *U1;
U20 = *U2;
U2D0 = sqrt(U10*U10 + U20*U20);
	if(*w_id==0){
	///////////////////  wall = 0
		*x1_old = xs;
		*x2_old = ys;
		*x2 = ys + (ys-*x2);
		*U2 = -*U2;
		*dt = len/U2D0;
	}
	else if(*w_id==1){
	///////////////////  wall = 1
		if(*dt>1e-15){
			*x1_old = xs;
			*x2_old = ys;
			dt_rm = *dt - len/U2D0;
			U10 = -sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gennew)));
			U20 =  sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal ( gennew );//diffusive
			U30 =  sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal ( gennew );//diffusive
			UI = U10;
			VI = U20;
			WI = U30;
			DX = UI*dt_rm;
			DY = VI*dt_rm;
			DZ = WI*dt_rm;
			*x1 = xs + DX;
			*x2 = sqrt( (ys+DY)*(ys+DY)+DZ*DZ );
			*U1 = UI;
			*U2 = (VI*(ys+DY)+WI*DZ)/(*x2);
			*U3 = (WI*(ys+DY)-VI*DZ)/(*x2);
			*dt = *dt - dt_rm;
		}
		else{
			*x1 = (*box).L1[0] - eps;
			*w_id = -1;
		}
	}
	else if(*w_id==2){
	///////////////////  wall = 2
		if(*dt>1e-15){
			*x1_old = xs;
			*x2_old = ys;
			dt_rm = *dt - len/U2D0;
			U10 = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gennew );//diffusive
			U20 = sqrt(2.0* (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * sqrt(-log(uniform(gennew)));//diffusive
			U30 = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gennew );//diffusive
			UI = U10;
			VI = U20;
			WI = U30;
			DX = UI*dt_rm;
			DY = VI*dt_rm;
			DZ = WI*dt_rm;
			*x1 = xs + DX;
			*x2 = sqrt( (ys+DY)*(ys+DY)+DZ*DZ );
			*U1 = UI;
			*U2 = (VI*(ys+DY)+WI*DZ)/(*x2);
			*U3 = (WI*(ys+DY)-VI*DZ)/(*x2);
			*dt = *dt - dt_rm;
		}
		else{
			*x2 = (*box).L2[0] + eps;
			//*x1 = xs;
			*w_id=-1;
		}
	}
	else if(*w_id==3){
	///////////////////  wall = 3
		*flag = 0;
	}
	else if(*w_id==4){
	///////////////////  wall = 4
		*flag = 0;
	}
	else if(*w_id==5){
	///////////////////  wall = 5
		*flag = 0;
	}
}
// end of function
}
void dcones_wall_collision(double *U1,double *U2,double *U3,double *x1,double *x2, struct BOX *box, struct CELLS *cells, int* flag, double *x1_old, double *x2_old, int *w_id, double xs, double ys, double len, double *dt, struct GAS *gas){
// Update velocity and position after collision with the corresponding wall
/*
// id of walls:
	       @
	      /|
	     / |
	  4 /  |3
	   /   |
	  /    @
	 @    /
	/    /
      5/    /2
      /    /
  6  /    @
 @--@	 /
7|	/1
 |     /
 @----@
   0
*/
//  int i;
if(*w_id == -1){
	*dt = 0.0;
}
else{
double U10, U20, U30, U2D0;
double dt_rm;
double DX, DY, DZ, UI,VI,WI;
  std::random_device rdnew;
  std::mt19937 gennew(rdnew());

  std::uniform_real_distribution<> uniform(0.0, 1.0);
  std::normal_distribution<> Normal(0.0,1.0);

U10 = *U1;
U20 = *U2;
U2D0 = sqrt(U10*U10 + U20*U20);
	if(*w_id==0){
	///////////////////  wall = 0
		*x1_old = xs;
		*x2_old = ys;
		*x2 = ys + (ys-*x2);
		*U2 = -*U2;
		*dt = len/U2D0;
	}
	else if(*w_id==1){
	///////////////////  wall = 1
		if(*dt>1e-15){
			*x1_old = xs;
			*x2_old = ys;
			dt_rm = *dt - len/U2D0;
			U10 = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal ( gennew );//diffusive
			U20 = sqrt(2.0* (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * sqrt(-log(uniform(gennew)));//diffusive
			U30 = sqrt( (*gas).kb*(*box).T_wall_1/( (*gas).m) ) * Normal ( gennew );//diffusive
			// rotating access R(theta)
			rotate(&U10, &U20, (*box).alpha[0], &UI, &VI);
			WI = U30;
			DX = UI*dt_rm;
			DY = VI*dt_rm;
			DZ = WI*dt_rm;
			*x1 = xs + DX;
			*x2 = sqrt( (ys+DY)*(ys+DY)+DZ*DZ );
			*U1 = UI;
			*U2 = (VI*(ys+DY)+WI*DZ)/(*x2);
			*U3 = (WI*(ys+DY)-VI*DZ)/(*x2);
			*dt = *dt - dt_rm;
		}
		else{
			*x2 = (xs-(*box).Lx[0])*tan( (*box).alpha[0] ) + 1e-15;
			*x1 = xs;
			*w_id = -1;
		}
	}
	else if(*w_id==2){
	///////////////////  wall = 2
		if(*dt>1e-15){
			*x1_old = xs;
			*x2_old = ys;
			dt_rm = *dt - len/U2D0;
			U10 = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gennew );//diffusive
			U20 = sqrt(2.0* (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * sqrt(-log(uniform(gennew)));//diffusive
			U30 = sqrt( (*gas).kb*(*box).T_wall_2/( (*gas).m) ) * Normal ( gennew );//diffusive
			// rotating access R(theta)
			rotate(&U10, &U20, (*box).alpha[1], &UI, &VI);
			WI = U30;
			DX = UI*dt_rm;
			DY = VI*dt_rm;
			DZ = WI*dt_rm;
			*x1 = xs + DX;
			*x2 = sqrt( (ys+DY)*(ys+DY)+DZ*DZ );
			*U1 = UI;
			*U2 = (VI*(ys+DY)+WI*DZ)/(*x2);
			*U3 = (WI*(ys+DY)-VI*DZ)/(*x2);
			*dt = *dt - dt_rm;
		}
		else{
			*x2 = (*box).Lx[1]*tan( (*box).alpha[0] ) + (xs-(*box).Lx[0]-(*box).Lx[1])*tan( (*box).alpha[1] ) + 1e-15;
			*x1 = xs;
			*w_id=-1;
		}
	}
	else if(*w_id==3){
	///////////////////  wall = 3
		*flag = 0;
	}
	else if(*w_id==4){
	///////////////////  wall = 4
		*flag = 0;
	}
	else if(*w_id==5){
	///////////////////  wall = 5
		*flag = 0;
	}
	else if(*w_id==6){
	///////////////////  wall = 6
		*flag = 0;
	}
	else if(*w_id==7){
	///////////////////  wall = 7
		*flag = 0;
	}
}
// end of function
}
int new_particles(struct GAS *gas, struct BOX *box,  struct CELLS *cells){
double beta, s;

double dummy;
double area;
double pi = acos(-1.0);
double R1, R2, h;
// generate new partiles for for walls 4,5,6,7
// wall 7 gets input temperature, density and velocity according to the Mach
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform(0.0, 1.0);
int num, ng, j;
double Un;
int Nx, Ny;
if( (*box).problem == dcones ){
	Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
	Ny = (*box).Ny;
}
else if ( (*box).problem == flatnose ){
	Nx = (*box).N1[0] + (*box).N1[1];
	Ny = (*box).N2[0] + (*box).N2[1];
}
for(num=0; num<Ny*Nx; num++){
	cells[num].ng[0] = 0;
	cells[num].ng[1] = 0;
}
if ( (*box).problem == flatnose ){
	for(int w_id=5; w_id>=4; w_id--){
		if(w_id==5){
			(*box).Ndot_w5 = 0;
			(*box).N_tot_genenrate = 0;

//			area = pi*pow((*box).L2[0]+(*box).L2[1],2.0);
			area = 2.0*pi*( (*box).L2[0]+(*box).L2[1] );
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_5) );
			s = (*box).U_wall_5*beta;
			dummy = (*box).n_wall_5/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/(2.0*sqrt(pi)+uniform(gen));
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
			(*box).Ndot_w5 = ng;
			(*box).N_tot_genenrate = ng;
/*
			for(j=0; j<Ny; j++){
				num = j*Nx + 0;
				area = pi*( cells[num].dim[3]*cells[num].dim[3] - cells[num].dim[0]*cells[num].dim[0]);
//				area = 2.0*pi*cells[num].dy*( cells[num].dim[2] + cells[num].dim[0])/2.0;
//				area = 2.0*pi*cells[num].dy;
				beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_5) );
				s = (*box).U_wall_5*beta;
				dummy = (*box).n_wall_5/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/(2.0*sqrt(pi)+uniform(gen));
				dummy = dummy*(*gas).delta_t*area;
				ng = floor(dummy/(*gas).Fn);
				cells[ num ].ng[0] =ng;
				(*box).Ndot_w5 += ng;
				(*box).N_tot_genenrate += ng;
			}
*/
		}
		else if(w_id==4){
			(*box).Ndot_w4 = 0;
			for(j=0; j<Nx; j++){
				num = (Ny-1)*Nx + j;
				if(cells[num].num_inside>1){
					area = 2.0*pi*( (*box).L2[0]+(*box).L2[1] )*cells[num].dx;
					beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
					s = -cells[num].U_space[1]*beta;
					dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
					dummy = dummy*(*gas).delta_t*area;
					ng = floor(dummy/(*gas).Fn);
					cells[ num ].ng[1] =ng;
					(*box).Ndot_w4 += ng;
					(*box).N_tot_genenrate += ng;
				}
			}
		}
	}
}
else if( (*box).problem == dcones ){
for(int w_id=7; w_id>=4; w_id--){
	if(w_id==7){
		(*box).Ndot_w7 = 0;
		(*box).N_tot_genenrate = 0;
		for(j=0; j<(*box).Ny; j++){
			num = j*Nx + 0;
			area = pi*(cells[num].dim[0]+cells[num].dim[2])*cells[num].dy;
			//area = pi*( cells[num].dim[2]*cells[num].dim[2]-cells[num].dim[0]*cells[num].dim[0] );
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_7) );
			s = (*box).U_wall_7*beta;
			dummy = (*box).n_wall_7/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
			cells[ num ].ng[0] =ng;
			(*box).Ndot_w7 += ng;
			(*box).N_tot_genenrate += ng;
			//printf("ng=%d in j=%d\n",ng,j);
		}
	}
	else if(w_id==6){
		(*box).Ndot_w6 = 0;
		for(j=0; j<(*box).Nx[0]; j++){
			num = ((*box).Ny-1)*Nx + j;

			area = 2.0*pi*(*box).Ly*cells[num].dx;
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			s = -cells[num].U_space[1]*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
			cells[ num ].ng[1] =ng;
			(*box).Ndot_w6 += ng;
			(*box).N_tot_genenrate += ng;
		}
	}
	else if(w_id==5){
		(*box).Ndot_w5 = 0;
		for(j=0; j<(*box).Nx[1]; j++){
			num = ((*box).Ny-1)*Nx + (*box).Nx[0]+j;
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			Un = cells[num].U_space[1]*cos((*box).alpha[0])-cells[num].U_space[0]*sin((*box).alpha[0]);
			s = -Un*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			R1 = cells[num].dim[3];
			R2 = cells[num].dim[2];
			h =  cells[num].dx;
			area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
			cells[ num ].ng[1] =ng;
			(*box).Ndot_w5 += ng;
			(*box).N_tot_genenrate += ng;

		}
	}
	else if(w_id==4){
		(*box).Ndot_w4 = 0;
		for(j=0; j<(*box).Nx[2]; j++){
			num = ((*box).Ny-1)*Nx + (*box).Nx[0]+(*box).Nx[1]+j;
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			Un = cells[num].U_space[1]*cos((*box).alpha[1])-cells[num].U_space[0]*sin((*box).alpha[1]);
			s = -Un*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			R1 = cells[num].dim[3];
			R2 = cells[num].dim[2];
			h =  cells[num].dx;
			area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
			cells[ num ].ng[1] =ng;
			(*box).Ndot_w4 += ng;
			(*box).N_tot_genenrate += ng;
		}
	}
}
}



/*
beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_7) );
s = (*box).U_wall_7*beta;
dummy = (*box).n_wall_7/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
area = pi*(*box).Ly*(*box).Ly;
dummy = dummy*(*gas).delta_t*area;
(*box).Ndot_w7 = floor(dummy/(*gas).Fn);

beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_6) );
s = 0.0;//(*box).U_wall_6*beta;
dummy = (*box).n_wall_6/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
area = 2.0*pi*(*box).Ly*(*box).Lx[0];
dummy = dummy*(*gas).delta_t*area;
(*box).Ndot_w6 = floor(dummy/(*gas).Fn);

beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_5) );
s = 0.0;//(*box).U_wall_5*beta;
dummy = (*box).n_wall_5/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
R1 = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1];
R2 = (*box).Ly;
h =  (*box).Lx[1];
area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
//area = pi/sin((*box).alpha[0])*( pow((*box).Ly+tan((*box).alpha[0])*(*box).Lx[1],2.0)-pow((*box).Ly,2.0) );
dummy = dummy*(*gas).delta_t*area;
(*box).Ndot_w5 = floor(dummy/(*gas).Fn);

beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_4) );
s = 0.0;//(*box).U_wall_4*beta;
dummy = (*box).n_wall_4/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
R1 = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1]+tan((*box).alpha[1])*(*box).Lx[2];
R2 = (*box).Ly+tan((*box).alpha[0])*(*box).Lx[1];
h =  (*box).Lx[2];
area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
//area = pi/sin((*box).alpha[1])*( pow((*box).Ly+tan((*box).alpha[0])*(*box).Lx[1]+tan((*box).alpha[1])*(*box).Lx[2],2.0) - pow((*box).Ly+tan((*box).alpha[0])*(*box).Lx[1],2.0) );
dummy = dummy*(*gas).delta_t*area;
(*box).Ndot_w4 = floor(dummy/(*gas).Fn);

(*box).N_tot_genenrate = floor( (*box).Ndot_w4 + (*box).Ndot_w5 + (*box).Ndot_w6 + (*box).Ndot_w7);
*/
return (*box).N_tot_genenrate;
}
void  generate_new_particles(double *x1n,double *x2n,double *U1n,double *U2n,double *U3n,double *x1_oldn,double *x2_oldn, struct GAS *gas, struct BOX *box, double *dtn,  struct CELLS *cells){
/*
for(i=0; i<(*gas).N; i++){
	x1n[i] = x1[i];
	x2n[i] = x2[i];
	U1n[i] = U1[i];
	U2n[i] = U2[i];
	U3n[i] = U3[i];
	Fnn[i] = Fn[i];
	x1_oldn[i] = x1_old[i];
	x2_oldn[i] = x2_old[i];
}
*/
double std;
double UI, VI;
double XI, YI;
double eps = 1e-15;
int w_id,i, num;
double dummy;
double uw;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis_x(0.0, 1.0);
std::uniform_real_distribution<> uniform(0.0, 1.0);

 int j=0;

//int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
int Nx, Ny;
if( (*box).problem == dcones ){
	Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
	Ny = (*box).Ny;
}
else if ( (*box).problem == flatnose ){
	Nx = (*box).N1[0] + (*box).N1[1];
	Ny = (*box).N2[0] + (*box).N2[1];
}

int ng;
double QA, FS2;
double VIr, UIr;
int n_tot=0;
double VMP, SC, FS1,  U,  FS1m;
int   n_inflow;
double pratio;
if ( (*box).problem == flatnose ){
  for(w_id=5; w_id>=4; w_id--){
	if(w_id==5){
		// generate new velocites from a
		if((*box).step > (*box).init_step)
			uw = (*box).U_wall_5;//(*box).U_wall_5;
		else
			uw = (*box).U_wall_5;

		VMP = sqrt(2.0*(*gas).kb*(*box).T_wall_5/(*gas).m);
		SC = (*box).U_wall_5/VMP;
		QA = 3.0;
		if(SC<-3.0)
			QA = fabs(SC)+1.0;
		FS1 =  SC+sqrt(SC*SC+2.0);//SC + sqrt(2.0+SC*SC);
		FS1m = SC-sqrt(SC*SC+2.0);
		FS2 = 0.5*(1.0+SC*(2.0*SC-FS1));//0.5*(1.0+SC*(2.0*SC-FS1));
		n_tot=0;
		std = sqrt((*gas).kb*(*box).T_wall_5/(*gas).m);
		std::normal_distribution<> dis_u5(0.0,std);
		ng = (*box).Ndot_w5;
		for(i=n_tot; i<n_tot+ng; i++){
			U2n[i] = dis_u5(gen);
			U3n[i] = dis_u5(gen);
			U = (10.0*uniform(gen)*FS1-SC);
			pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
			while(pratio<uniform(gen)){
				U = (10.0*uniform(gen)*FS1-SC);
				pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
			}
			U1n[i] = (U+SC)*VMP;
			x1n[i] = eps;
			x2n[i] = ((*box).L2[0]+(*box).L2[1])*sqrt( uniform(gen) );
		}
/*
		for(j=0; j<Ny; j++){
			num = j*Nx + 0;
			ng = cells[num].ng[0];
			for(i=n_tot; i<n_tot+ng; i++){
				U2n[i] = dis_u5(gen);
				U3n[i] = dis_u5(gen);
				U = (10.0*uniform(gen)*FS1-SC);
				pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
				while(pratio<uniform(gen)){
					U = (10.0*uniform(gen)*FS1-SC);
					pratio = (2.0*(U+SC))*exp( 0.5+(SC/2.0)*FS1m-U*U )/FS1;
				}
				U1n[i] = (U+SC)*VMP;

//				yes = 0;
//				while(yes == 0){
//					U = -QA+2.0*QA*uniform(gen);
//					UN=U+SC;
//					if(UN>0.0){
//						A=(2.*UN/FS1)*exp(FS2-U*U);
//						if(A>uniform(gen)){
//							U1n[i] = UN*VMP;//+(*box).U_wall_5;
//							yes = 1;
//						}
//					}
//				}
//
				x1n[i] = eps;
//				x2n[i] = cells[num].dim[0] + dis_x(gen)*cells[num].dy;
				x2n[i] = sqrt( pow(cells[num].dim[0],2)+uniform(gen)*(pow(cells[num].dim[3],2)-pow(cells[num].dim[0],2)) );
			}

			n_tot += ng;
		}
*/
/*
		std = sqrt((*gas).kb*(*box).T_wall_5/(*gas).m);
		std::normal_distribution<> dis_u5(0.0,std);
		n_tot=0;
		for(j=0; j<Ny; j++){
			num = j*Nx + 0;
			ng = cells[num].ng[0];
			for(i=n_tot; i<n_tot+ng; i++){
				U2n[i] = dis_u5(gen);
				U3n[i] = dis_u5(gen);
				dummy = -1.0;
				while(dummy<0.0)
					dummy = dis_u5(gen) + uw;
				U1n[i] = dummy;
				//U1n[i] =dis_u5(gen);
				//U1n[i] = sqrt(2.0* (*gas).kb*(*box).T_wall_5/( (*gas).m) ) * sqrt(-log(uniform(gen)))+uw;
				x1n[i] = eps;
				x2n[i] = cells[num].dim[0] + dis_x(gen)*cells[num].dy;
			}
			n_tot += ng;
		}
*/
		n_tot += ng;
		n_inflow = n_tot;
	}
	else if(w_id==4){
		// generate new velocites from a
		for(j=0; j<Nx; j++){
			num = (Ny-1)*Nx + j;
			VMP = sqrt(2.0*(*gas).kb*(*box).T_wall_5/(*gas).m);
			//std = sqrt((*gas).kb*cells[num].T/(*gas).m);
			std = sqrt((*gas).kb*(*box).T_wall_5/(*gas).m);
			std::normal_distribution<> dis_u4(0.0,std);
			ng = cells[num].ng[1];
			if(ng>0){
				for(i=n_tot; i<n_tot+ng; i++){

					U1n[i] = dis_u4(gen)+cells[num].U_space[0];
					dummy = 1.0;
					while(dummy>0.0)
						dummy = dis_u4(gen)+cells[num].U_space[1];
					U2n[i] = dummy;
					U3n[i] = dis_u4(gen)+cells[num].U_space[2];

/*
					U1n[i] = dis_u4(gen);
					U2n[i] = -sqrt(-log(uniform(gen)))*VMP;
					U3n[i] = dis_u4(gen);
*/
					x1n[i] = cells[num].cell_center[0]-cells[num].dx/2.0+dis_x(gen)*cells[num].dx;
					x2n[i] = (*box).L2[0]+(*box).L2[1]-eps;
				}
				n_tot += ng;
			}
		}
	}
  }
}
else if ( (*box).problem == dcones ){
  for(w_id=7; w_id>=4; w_id--){
	if(w_id==7){
		// generate new velocites from a
		if((*box).step > (*box).init_step)
			uw = (*box).U_wall_7;
		else
			uw = 0.0;
		std = sqrt((*gas).kb*(*box).T_wall_7/(*gas).m);
		std::normal_distribution<> dis_u7(0.0,std);
//		std::normal_distribution<> dis_uw(uw,std);
//		(*box).Ndot_w7 = 0;
//		(*box).N_tot_genenrate = 0;
		n_tot=0;
		for(j=0; j<(*box).Ny; j++){
			num = j*Nx + 0;
/*
			area = pi*( cells[num].dim[2]*cells[num].dim[2]-cells[num].dim[0]*cells[num].dim[0] );
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*(*box).T_wall_7) );
			s = (*box).U_wall_7*beta;
			dummy = (*box).n_wall_7/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
*/
			ng = cells[num].ng[0];
			for(i=n_tot; i<n_tot+ng; i++){
				U2n[i] = dis_u7(gen);
				U3n[i] = dis_u7(gen);
				dummy = -1.0;
				while(dummy<0.0)
					dummy = dis_u7(gen) + uw;//dis_uw(gen);
				U1n[i] = dummy;
				x1n[i] = eps;
				x2n[i] = cells[num].dim[0] + dis_x(gen)*(cells[num].dim[2]-cells[num].dim[0]);
			}
			n_tot += ng;
//			(*box).Ndot_w7 += ng;
//			(*box).N_tot_genenrate += ng;
		}
	}
	else if(w_id==6){
		// generate new velocites from a
		// generate new velocites from a
		if((*box).step > (*box).init_step)
			uw = (*box).U_wall_6;
		else
			uw = 0.0;
//		(*box).Ndot_w6 = 0;
		for(j=0; j<(*box).Nx[0]; j++){
			num = ((*box).Ny-1)*Nx + j;
			std = sqrt((*gas).kb*cells[num].T/(*gas).m);
			std::normal_distribution<> dis_u6(0.0,std);
/*
			area = 2.0*pi*(*box).Ly*cells[num].dx;
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			s = -cells[num].U_space[1]*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
*/
			ng = cells[num].ng[1];
			if(ng>0){
				for(i=n_tot; i<n_tot+ng; i++){
					U1n[i] = dis_u6(gen)+cells[num].U_space[0];//+uw;
					dummy = 1.0;
					while(dummy>0.0)
						dummy = dis_u6(gen)-cells[num].U_space[1];
//					U2n[i] = -sqrt(2.0* (*gas).kb*(*box).T_wall_6/( (*gas).m) ) * sqrt(-log(uniform(gen)));
					U2n[i] = dummy;
					U3n[i] = dis_u6(gen)+cells[num].U_space[2];
					x1n[i] = j*cells[num].dx+dis_x(gen)*cells[num].dx;
					x2n[i] = (*box).Ly-eps;
				}
				n_tot += ng;
//				(*box).Ndot_w6 += ng;
//				(*box).N_tot_genenrate += ng;
			}
		}
	}
	else if(w_id==5){
		// generate new velocites from a
		if((*box).step > (*box).init_step)
			uw = (*box).U_wall_5;
		else
			uw = 0.0;
//		(*box).Ndot_w5 = 0;
		for(j=0; j<(*box).Nx[1]; j++){
			num = ((*box).Ny-1)*Nx + (*box).Nx[0]+j;
			std = sqrt((*gas).kb*cells[num].T/(*gas).m);
			std::normal_distribution<> dis_u5(0.0,std);
/*
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			Un = cells[num].U_space[1]*cos((*box).alpha[0])-cells[num].U_space[0]*sin((*box).alpha[0]);
			s = -Un*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			R1 = cells[num].dim[3];
			R2 = cells[num].dim[2];
			h =  cells[num].dx;
			area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
//			area = pi/sin((*box).alpha[0])*( pow((*box).Ly+tan((*box).alpha[0])*(*box).Lx[1],2.0)-pow((*box).Ly,2.0) );
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
*/
			ng = cells[num].ng[1];
			if(ng>0){
				for(i=n_tot; i<n_tot+ng; i++){
/*
					dummy = 1.0;
					while (dummy>1e-15){
						UI = dis_u5(gen)+cells[num].U_space[0];
						VI = dis_u5(gen)+cells[num].U_space[1];
						dummy = VI*cos((*box).alpha[0])-UI*sin((*box).alpha[0]);
					}
					U1n[i] = UI;
					U2n[i] = VI;
					U3n[i] = dis_u5(gen)+cells[num].U_space[0];;
*/
					rotate(&cells[num].U_space[0], &cells[num].U_space[1], -(*box).alpha[0], &UIr, &VIr);
					UI = dis_u5(gen)+UIr;
					dummy = 1.0;
					while(dummy>0.0)
						dummy = dis_u5(gen)-VIr;
					VI = dummy;
//					VI = -sqrt(2.0* (*gas).kb*(*box).T_wall_5/( (*gas).m) ) * sqrt(-log(uniform(gen)));
				// rotating access R(theta)
					rotate(&UI, &VI, (*box).alpha[0], &U1n[i], &U2n[i]);
					U3n[i] = dis_u5(gen)+cells[num].U_space[2];

					XI = dis_x(gen)*cells[num].dx/cos((*box).alpha[0]);
					YI = 0.0;
					// rotating points R(theta)
					rotate(&XI, &YI, (*box).alpha[0], &x1n[i], &x2n[i]);

					x1n[i] = x1n[i] + (*box).Lx[0] + j*(*box).Lx[1]/(1.0*(*box).Nx[1]);
					x2n[i] = x2n[i] + cells[num].dim[2]-eps;
				}
				n_tot += ng;
//				(*box).Ndot_w5 += ng;
//				(*box).N_tot_genenrate += ng;
			}
		}
	}
	else if(w_id==4){
		// generate new velocites from a
		if((*box).step > (*box).init_step)
			uw = (*box).U_wall_4;
		else
			uw = 0.0;
//		(*box).Ndot_w4 = 0;
		for(j=0; j<(*box).Nx[2]; j++){
			num = ((*box).Ny-1)*Nx + (*box).Nx[0]+(*box).Nx[1]+j;
			std = sqrt((*gas).kb*cells[num].T/(*gas).m);
			std::normal_distribution<> dis_u4(0.0,std);
/*
			beta = sqrt( (*gas).m/(2.0*(*gas).kb*cells[num].T) );
			Un = cells[num].U_space[1]*cos((*box).alpha[1])-cells[num].U_space[0]*sin((*box).alpha[1]);
			s = -Un*beta;
			dummy = cells[num].n/beta*( exp(-s*s)+sqrt(pi)*s*(1.0+erf(s)) )/2.0/sqrt(pi);
			R1 = cells[num].dim[3];
			R2 = cells[num].dim[2];
			h =  cells[num].dx;
			area = pi*(R1+R2)*sqrt((R1-R2)*(R1-R2)+h*h);
			dummy = dummy*(*gas).delta_t*area;
			ng = floor(dummy/(*gas).Fn);
*/
			ng = cells[num].ng[1];
			if(ng>0){
				for(i=n_tot; i<n_tot+ng; i++){
/*
					dummy = 1.0;
					while (dummy>1e-15){
						UI = dis_u4(gen)+uw;
						VI = dis_u4(gen);
						dummy = VI*cos((*box).alpha[1])-UI*sin((*box).alpha[1]);
					}
					U1n[i] = UI;
					U2n[i] = VI;
					U3n[i] = dis_u4(gen);
*/
					rotate(&cells[num].U_space[0], &cells[num].U_space[1], -(*box).alpha[1], &UIr, &VIr);
					UI = dis_u4(gen)+UIr;
					dummy = 1.0;
					while(dummy>0.0)
						dummy = dis_u4(gen)-VIr;
					VI = dummy;
//					VI = -sqrt(2.0* (*gas).kb*(*box).T_wall_4/( (*gas).m) ) * sqrt(-log(uniform(gen)));
					// rotating access R(theta)
					rotate(&UI, &VI, (*box).alpha[1], &U1n[i], &U2n[i]);
					U3n[i] = dis_u4(gen)+cells[num].U_space[2];;

					XI = dis_x(gen)*cells[num].dx/cos((*box).alpha[1]);
					YI = 0.0;
					// rotating points R(theta)
					rotate(&XI, &YI, (*box).alpha[1], &x1n[i], &x2n[i]);

					x1n[i] = x1n[i] + (*box).Lx[0]+(*box).Lx[1]+j*(*box).Lx[2]/(1.0*(*box).Nx[2]);
					x2n[i] = x2n[i] + cells[num].dim[2]-eps;
				}
				n_tot += ng;
//				(*box).Ndot_w4 += ng;
//				(*box).N_tot_genenrate += ng;
			}
		}
	}
  }
}

// check the mean of input particles
(*gas).n_inflow = n_inflow;
/*
double mn[3], st[3], mnx2n;
mn[0] = mean(U1n, n_inflow);
mn[1] = mean(U2n, n_inflow);
mn[2] = mean(U3n, n_inflow);
st[0] = (*gas).m*standard_deviation(U1n, n_inflow, &mn[0])*standard_deviation(U1n, n_inflow, &mn[0])/(*gas).kb;
st[1] = (*gas).m*standard_deviation(U2n, n_inflow, &mn[1])*standard_deviation(U2n, n_inflow, &mn[1])/(*gas).kb;
st[2] = (*gas).m*standard_deviation(U3n, n_inflow, &mn[2])*standard_deviation(U3n, n_inflow, &mn[2])/(*gas).kb;
*/
// evolve new particles
//double XI, YI;
double  WI;
double X,Y;
double DX,DY,DZ;

//double dx[3];
//  dx[0] = (*box).Lx[0]/(1.0*(*box).Nx[0]);
//  dx[1] = (*box).Lx[1]/(1.0*(*box).Nx[1]);
//  dx[2] = (*box).Lx[2]/(1.0*(*box).Nx[2]);
/*
int ndummy=0;
for(num=0; num<(*box).Ny*Nx; num++){
		ndummy += cells[num].ng[0];
		ndummy += cells[num].ng[1];
}
printf("n_tot=%d and (*box).N_tot_genenrate=%d\n, n_dummy=%d",n_tot,(*box).N_tot_genenrate,ndummy);
*/
for(i=0; i<(*box).N_tot_genenrate; i++){
	x1_oldn[i] = x1n[i];
	x2_oldn[i] = x2n[i];
/*
	if(x1n[i]<(*box).Lx[0]){
		j = floor(x1n[i]/dx[0]);
    		k = floor( log(1-(1-(*box).s)*x2n[i]/(*box).dy0) / log((*box).s) );
	}
	else if(x1n[i]<(*box).Lx[0]+(*box).Lx[1]){
		j = (*box).Nx[0]+floor((x1n[i]-(*box).Lx[0])/dx[1]);
		k = floor( log(1-(1-(*box).s)*(x2n[i]-(x1n[i]-(*box).Lx[0])*tan((*box).alpha[0]))/(*box).dy0) / log((*box).s) );
	}
	else{
		j = (*box).Nx[0]+(*box).Nx[1]+floor((x1n[i]-(*box).Lx[0]-(*box).Lx[1])/dx[2]);
		k = floor( log(1-(1-(*box).s)*(x2n[i]-(*box).Lx[1]*tan((*box).alpha[0])-(x1n[i]-(*box).Lx[0]-(*box).Lx[1])*tan((*box).alpha[1]))/(*box).dy0) / log((*box).s) );
	}
	num =  k*Nx + j;
    	if(j<0  )
      		printf("inside: OOPS, particle is out from left, x1=%e x2=%e, j = %d\n",x1n[i], x2n[i],j);
    	if(j>(*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2]-1 )
      		printf("inside: OOPS, particle is out from right, x1=%e x2=%e, j = %d\n",x1n[i], x2n[i],j);
    	if(k<0 )
      		printf("inside: OOPS, particle is out from down, x1=%e x2=%e, k = %d\n",x1n[i], x2n[i], k);
    	if(k>(*box).Ny-1 )
      		printf("inside: OOPS, particle is out from up, x1=%e x2=%e, k = %d\n",x1n[i], x2n[i],k);
*/
	dtn[i] = uniform(gen)*(*gas).delta_t;
	XI = x1n[i];
	YI = x2n[i];
	VI = U2n[i];
	WI = U3n[i];
	DX = U1n[i]*dtn[i];
	DY = VI*dtn[i];
	DZ = WI*dtn[i];
	X = XI + DX;
	Y = sqrt( (YI+DY)*(YI+DY)+DZ*DZ );
	U2n[i] = (VI*(YI+DY)+WI*DZ)/Y;
	U3n[i] = (WI*(YI+DY)-VI*DZ)/Y;
	x1n[i] = X;
	x2n[i] = Y;
}
/*
mn[0] = mean(U1n, n_inflow);
mn[1] = mean(U2n, n_inflow);
mn[2] = mean(U3n, n_inflow);
st[0] = (*gas).m*standard_deviation(U1n, n_inflow, &mn[0])*standard_deviation(U1n, n_inflow, &mn[0])/(*gas).kb;
st[1] = (*gas).m*standard_deviation(U2n, n_inflow, &mn[1])*standard_deviation(U2n, n_inflow, &mn[1])/(*gas).kb;
st[2] = (*gas).m*standard_deviation(U3n, n_inflow, &mn[2])*standard_deviation(U3n, n_inflow, &mn[2])/(*gas).kb;
*/
}
int flatnose_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x1_old,double *x2_old, int* flag, struct GAS *gas, struct BOX *box, struct CELLS *cells,
double *U1n,double *U2n,double *U3n,double *x1n,double *x2n,double *x1_oldn,double *x2_oldn, int* flagn, double *dtn){
// applying BC on double cone problem
/*
// id of walls:
	       @
	      /|
	     / |
	  4 /  |3
	   /   |
	  /    @
	 @    /
	/    /
      5/    /2
      /    /
  6  /    @
 @--@	 /
7|	/1
 |     /
 @----@
   0
*/
int i,w_id,out;
double len = 0.0;
double xs=0.0, ys=0.0;
double dt;
int N_in = 0;
//int num;

N_in = 0;
int Nx = (*box).N1[0] + (*box).N1[1];
int Ny = (*box).N2[0] + (*box).N2[1];
int j,k ;
//double U10, U20, U30, U2D0, dt_rm;
//double DX, DY, DZ, UI,VI,WI;
// BC check for old particles
//printf("check old particles\n");
for(i=0; i<(*gas).N; i++){
//	x1[i] = 0.0207128288195467;
//	x2[i] = 0.00165614204549755;
//	x1_old[i] = 0.0200435717620908;
//	x2_old[i] = 0.00141158023330257;
//	U1[i] = 0.0;
//	U2[i] = 66088.019;
	out = 1;
	dt = (*gas).delta_t;
	while(out == 1){
		w_id = -1;
		// find the first wall hit by the particle
		flatnose_intersection(x1_old[i], x2_old[i], x1[i], x2[i], &xs, &ys, &len, &w_id, box);
		// collision with that walls
		flatnose_wall_collision(&U1[i], &U2[i], &U3[i], &x1[i], &x2[i], box, cells, &flag[i], &x1_old[i], &x2_old[i], &w_id, xs, ys, len, &dt, gas);
		if(w_id == -1 || w_id > 2)
			out = 0;
	}
	if(flag[i] == 1){
		N_in ++;
	}
	// Check
		//num =  flatnose_index(x1[i], x2[i], cells, box, &j, &k);
    		if(j<0  && flag[i]==1)
      			printf("OOPS, particle is out from left, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(j>Nx-1 && flag[i]==1)
      			printf("OOPS, particle is out from right, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(k<0 && flag[i]==1)
      			printf("OOPS, particle is out from down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(k>Ny-1 && flag[i]==1)
      			printf("OOPS, particle is out from up, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(k<(*box).N2[0] && j >= (*box).N1[0] && flag[i]==1)
      			printf("OOPS, particle is out from right down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
}
// BC check for new particles
//printf("check new particles\n");
for(i=0; i<(*box).N_tot_genenrate; i++){
//	x1_oldn[i] = 9.486832e-08;
//	x1n[i] = 3.657813e-10;
//	x2_oldn[i] = 3.360686e-09;
//	x2n[i] = 3.362265e-09;
	out = 1;
	dt = dtn[i];
	flagn[i] = 1;
	while(out == 1){
		w_id = -1;
		// find the first wall hit by the particle
		flatnose_intersection(x1_oldn[i], x2_oldn[i], x1n[i], x2n[i], &xs, &ys, &len, &w_id, box);
		// collision with that walls
		flatnose_wall_collision(&U1n[i], &U2n[i], &U3n[i], &x1n[i], &x2n[i], box, cells, &flagn[i], &x1_oldn[i], &x2_oldn[i], &w_id, xs, ys, len, &dt, gas);
		if(w_id == -1 || w_id > 2)
			out = 0;
	}
	if(flagn[i] == 1){
		N_in ++;
	}
	// Check
		//num =  flatnose_index(x1n[i], x2n[i], cells, box, &j, &k);
    		if(j<0  && flagn[i]==1)
      			printf("OOPS, particle is out from left, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(j>Nx-1 && flagn[i]==1)
      			printf("OOPS, particle is out from right, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(k<0 && flagn[i]==1)
      			printf("OOPS, particle is out from down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(k>Ny-1 && flagn[i]==1)
      			printf("OOPS, particle is out from up, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(k<(*box).N2[0] && j >= (*box).N1[0] && flagn[i]==1)
      			printf("OOPS, particle is out from right down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
}

return N_in;
}

int dcones_BC(double *U1,double *U2,double *U3,double *x1,double *x2, double *x1_old,double *x2_old, int* flag, struct GAS *gas, struct BOX *box, struct CELLS *cells,
double *U1n,double *U2n,double *U3n,double *x1n,double *x2n,double *x1_oldn,double *x2_oldn, int* flagn, double *dtn){
// applying BC on double cone problem
/*
// id of walls:
	       @
	      /|
	     / |
	  4 /  |3
	   /   |
	  /    @
	 @    /
	/    /
      5/    /2
      /    /
  6  /    @
 @--@	 /
7|	/1
 |     /
 @----@
   0
*/
int i,w_id,out;
double len = 0.0;
double xs=0.0, ys=0.0;
double dt;
int N_in = 0;

N_in = 0;
int Nx = (*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2];
int j,k ;
  double dx[3];
//double U10, U20, U30, U2D0, dt_rm;
//double DX, DY, DZ, UI,VI,WI;


  dx[0] = (*box).Lx[0]/(1.0*(*box).Nx[0]);
  dx[1] = (*box).Lx[1]/(1.0*(*box).Nx[1]);
  dx[2] = (*box).Lx[2]/(1.0*(*box).Nx[2]);
// BC check for old particles
//printf("check old particles\n");
for(i=0; i<(*gas).N; i++){
//	x1_old[i] = 3.160829e-09;
//	x1[i] = 3.160828e-09;
//	x2_old[i] = 2.964188e-09;
//	x2[i] = 2.964240e-09;
//	U1[i] = 0.0;
//	U2[i] = 66088.019;
	out = 1;
	dt = (*gas).delta_t;
	while(out == 1){
		w_id = -1;
		// find the first wall hit by the particle
		dcones_intersection(x1_old[i], x2_old[i], x1[i], x2[i], &xs, &ys, &len, &w_id, box);
		// collision with that walls
		dcones_wall_collision(&U1[i], &U2[i], &U3[i], &x1[i], &x2[i], box, cells, &flag[i], &x1_old[i], &x2_old[i], &w_id, xs, ys, len, &dt, gas);
		if(w_id == -1 || w_id > 2)
			out = 0;
	}
	if(flag[i] == 1){
		N_in ++;
	}
	// Check
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
		//num =  k*Nx + j;
    		if(j<0  && flag[i]==1)
      			printf("OOPS, particle is out from left, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(j>(*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2]-1 && flag[i]==1)
      			printf("OOPS, particle is out from right, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(k<0 && flag[i]==1)
      			printf("OOPS, particle is out from down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
    		if(k>(*box).Ny-1 && flag[i]==1)
      			printf("OOPS, particle is out from up, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1[i], x2[i],x1_old[i], x2_old[i],flag[i],w_id);
}
// BC check for new particles
//printf("check new particles\n");
for(i=0; i<(*box).N_tot_genenrate; i++){
//	x1_oldn[i] = 9.486832e-08;
//	x1n[i] = 3.657813e-10;
//	x2_oldn[i] = 3.360686e-09;
//	x2n[i] = 3.362265e-09;
	out = 1;
	dt = dtn[i];
	flagn[i] = 1;
	while(out == 1){
		w_id = -1;
		// find the first wall hit by the particle
		dcones_intersection(x1_oldn[i], x2_oldn[i], x1n[i], x2n[i], &xs, &ys, &len, &w_id, box);
		// collision with that walls
		dcones_wall_collision(&U1n[i], &U2n[i], &U3n[i], &x1n[i], &x2n[i], box, cells, &flagn[i], &x1_oldn[i], &x2_oldn[i], &w_id, xs, ys, len, &dt, gas);
		if(w_id == -1 || w_id > 2)
			out = 0;
	}
	if(flagn[i] == 1){
		N_in ++;
	}
	// Check
		if(x1n[i]<(*box).Lx[0]){
		    	j = floor(x1n[i]/dx[0]);
    			k = floor( log(1-(1-(*box).s)*x2n[i]/(*box).dy0) / log((*box).s) );
			//if(k>0)
			//	printf("x2[%d] = %e\n", i, x2[i]);
		}
		else if(x1n[i]<(*box).Lx[0]+(*box).Lx[1]){
			j = (*box).Nx[0]+floor((x1n[i]-(*box).Lx[0])/dx[1]);
			k = floor( log(1-(1-(*box).s)*(x2n[i]-(x1n[i]-(*box).Lx[0])*tan((*box).alpha[0]))/(*box).dy0) / log((*box).s) );
		}
		else{
			j = (*box).Nx[0]+(*box).Nx[1]+floor((x1n[i]-(*box).Lx[0]-(*box).Lx[1])/dx[2]);
			k = floor( log(1-(1-(*box).s)*(x2n[i]-(*box).Lx[1]*tan((*box).alpha[0])-(x1n[i]-(*box).Lx[0]-(*box).Lx[1])*tan((*box).alpha[1]))/(*box).dy0) / log((*box).s) );
		}
		//num =  k*Nx + j;
    		if(j<0  && flagn[i]==1)
      			printf("OOPS, particle is out from left, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(j>(*box).Nx[0]+(*box).Nx[1]+(*box).Nx[2]-1 && flagn[i]==1)
      			printf("OOPS, particle is out from right, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(k<0 && flagn[i]==1)
      			printf("OOPS, particle is out from down, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
    		if(k>(*box).Ny-1 && flagn[i]==1)
      			printf("OOPS, particle is out from up, x1=%e x2=%e,x1_old=%e x2_old=%e, flag = %d, w_id =%d\n",x1n[i], x2n[i],x1_oldn[i], x2_oldn[i],flagn[i],w_id);
}

return N_in;

}
