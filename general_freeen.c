#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#include "in_proto.h"
#include "in_vdefs.h"
#include "io_interface.c"

#define NDIM 2			// number of dimensions. if == 2, then Ly is is small

typedef struct{
	VecI rs;
} substrateNode;

substrateNode *subNode;

int nSub;

/* names and folders */
char name[120];
char foldername[120];
char fullname_iter[120];
char fullname_Wfields[120];
char fullname_Wsub[120];
char fullname_free_en[120];
char fullname_conv[120];
char fullname_snap[120];
char fullname_param[120];
char fullname_press[120];
char fullname_profiling[120];
char fullname_test[120];
char path[80];
char folder[80];
int write_snapshot;

/* global properties of the system*/
double kbT;								// kbT :-)
double rho_liq;						// density in the bulk (at the top of the box)
double lambda;							// relaxation parameter of the SC-scheme

/* properties of the run */
int	restart;							// restart flag
int	iteration_num;					// if restarted, the last existing iteration
int	iterations;						// desired iterations limit
int	filled_init;					// initial state of the interface
double convergence;					// convergence

/* size and grid related properties */
VecR 	L;									// box length
VecI 	nRep;								// number of replicas of the initial box
VecI 	grid;								// number of grid points
int subGridX;							// MPI-related to keep grid.x value for the substrate
int globStartIdx;						// MPI-related global starting index in a CPU
VecR 	corr;								// dimensions of the corrugation
VecR 	rho_fd, rho_bd, rho_cd;		// forw, backw and central density differences
double dx,dy,dz;						// size of the grid cell
double dx2, dy2, dz2;				// size of the grid cell squared
double volume;							// volume of the computation domain
double nMol;							// number of water molecules for NVT-ensemble
int 	nDims;							// number of dimensions
int 	NVT;								// canonical ensemble

/* substrate' parameters */
double rCut;							// cut off distance for potential
double rrCut;							// squared cut off distance for potential
double sigma_sub;						// sigma of the substrate potential
double sigma_sub2;					// sigma of the substrate potential squared
double eps;								// amplitude of the potential

/* statistical properties of the ensemble */
double total_N;						// total number of molecules in the system
double Q;								// statistical sum

/* time variables for simple profiling */
time_t begin, end;					// two marks of the time-measurement
double seconds;						// length of time interval between measurements

/* integration variables */
double *koeff_x_semi, *koeff_y_semi, *koeff_z_open;
double *XYa,*YYa;
double *XYat,*YYat;

/* substrate potential variables */
double *subPotZ;

/* relevant fields */
double ***wa;
double ***sub;
double ***waout;
double ***rho;

/* free energy terms (may be omitted in principal) */
double ***df_drho;
double ***grad;

#define TWRITE 50000

/* ======================================================== */

#define	V11	  -1.086641506142464
#define	W111		0.023102120829070
#define	KAPPA	2 * 0.015739171

/*
 c        On periodic boundary conditions:
 c       |<------------- L.x ------------->|
 c       |                                 |
 c    0__|__1__2__.....__grid.x-1__grid.x__|__grid.x+1
 c       |                                 |
 c       |                                 |
 c
 c       0 = grid.x, grid.x+1 = 1   _____ = L.x / grid.x
 c
 c    Make sure since the integrations will be performed by
 c    Simpson method grid.x grid.y grid.z must be odd numbers
 */

int myRank;
int MPIsize;

int main (int argc, char **argv){
	time(&begin);
	
	/* Start up MPI */
	MPI_Init(&argc, &argv);
	
	/* Find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	
	PassConsoleParams (argc, argv);

	double r_sub, ri2, ri6;
	double distance2;
	VecR dist2;
	int i1,j1,k1,kk;
	int is, js, ks, sn;
	int subDepth;
	double maxSubPot;
	VecR cm0;
	double rrToCm0;
	
	DivideSpace ();
	AllocArrays ();
	
	FILE  *Wtest, *Wrestart, *Wfields, *Wsub;
	FILE	*snapshot, *iterkeeper, *converge, *wPressure;
		
	/* Lengths are measured in nm */
	iterations = 60000000;
	kbT = 4.116404397;
	rCut = 5.;
	sigma_sub = .3;

	/* Squared lengths */
	sigma_sub2 = SQR(sigma_sub);
	rrCut = SQR(rCut);

	strcpy (foldername, path);
	strcat (foldername, folder);
	
	if (myRank == 0) {
		InitParameters ();				// read parameters from input line
		StoreParameters ();				// save all parameters to the file
	}
	
	DefineSimpson ();
	
	/* set initial values of field, gradients, etc. to zero. 
	 Create nodes of flat substrate */
	nSub = 0;
	subDepth = (int) (rCut / dz);
	maxSubPot = 0.;

	// calculate total number of box replicas in x,y-directions,
  // including the original box
	V_SET (nRep, (int) (2 * ((int) (rCut / L.x) + 1) + 1),
							 (int) (2 * ((int) (rCut / L.y) + 1) + 1),
							 0);
	printf ("Use %d parent cell replicas in x- and %8d in y-dir.\n", nRep.x, nRep.y);
	
	// set all the fields to zero
	for (i1 = 0; i1 < grid.x+1; i1++){
		for (j1 = 0; j1 < grid.y+1; j1++){
			for (k1 = 0; k1 < grid.z+1; k1++){
				/* set free energy terms to zero*/
				SET_TO_ZERO (df_drho,i1,j1,k1);
				SET_TO_ZERO (grad,i1,j1,k1);
				SET_TO_ZERO (sub,i1,j1,k1);
				
				/* set initial field to zero */
				SET_TO_ZERO (wa,i1,j1,k1);
				
				subPotZ[k1] = 0.;
			}
		}
	}
	printf ("set up fields to zero\n");

	AllocSubstrate ();
	
	printf ("Finished substrate' allocation\n");
	
	if (restart == 0) {
		/* set up the flat substrate */
		for (i1 = 1; i1 < nRep.x * subGridX + 1; i1++){
			for (j1 = 1; j1 < nRep.y * grid.y + 1; j1++){
				// create "0"th substrate plane
				V_SET(subNode[nSub].rs, i1, j1, 0);
				++nSub;
			}
		}
		printf ("set up the flat substrate\n");
	
		/* calculate z-dependence of the potential created by a flat substrate */
		for (k1 = 1; k1 < grid.z+1; k1++){
			for(sn = 0; sn < nSub; sn++){
				dist2.x = SQR(subNode[sn].rs.x - ((int) ((nRep.x - 1) / 2) + 0.5) * subGridX - 1);
				dist2.y = SQR(subNode[sn].rs.y - ((int) ((nRep.y - 1) / 2) + 0.5) * grid.y - 1);
			
				/* take into account the "0"th plane and the planes up to rCut beneath */
				for (ks = 0; ks < subDepth; ks++){
					dist2.z = SQR(ks + k1);
					distance2 = dist2.x * dx2 + dist2.y * dy2 + dist2.z * dz2;
					if (distance2 < rrCut) {
						ri2 = sigma_sub2/distance2;
						ri6 = CUBE(ri2);
						subPotZ[k1] +=	ri6 * (ri6 - 1.);
					}
				}
			}
			maxSubPot = MAX(maxSubPot, subPotZ[k1]);
			/* in the end we have one dimensional (w.r.t. z) potential */
		}
		printf ("Made the 0th substrate layer. It has %d nodes\n", nSub);
	
		/* creation of corrugation nodes (if any) initial fields */
		nSub = 0;	// now it count only corrugation nodes
		
		iteration_num = 0;
		total_N = 0.;
		int filled_x_grids = 0;

		printf ("Creating substrate nodes and liquid field\n");
		printf ("corr.z is %8.4f\n", corr.z);
		
		// create corrugation sites (inkl. replicas of the parent cell)
		for (i1 = 1; i1 < nRep.x * subGridX + 1; i1++){
			for (j1 = 1; j1 < nRep.y * grid.y + 1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create substrate nodes taking into account corrugation */
					if ( (((i1 - 1) % subGridX) * dx < corr.x) &&
										 (((j1 - 1) % grid.y) * dy < corr.y) && 
										 ((k1 - .5)*dz < corr.z) ) {
						V_SET(subNode[nSub].rs, i1, j1, k1);
						++nSub;
					} else {
					}
				}
			}
		}
		
		/* create spherical (nDims == 3) or cylindrical droplet (nDims == 2) */
		if (nDims == 3) {
			switch (filled_init) {
				case 0: 	rrToCm0 = pow(.75 * nMol / (rho_liq * M_PI),1./3.);
							rrToCm0 += .02 * dx;				// increase the radius by 2%
							break;
				case 1:	rrToCm0 = pow(1.5 * nMol / (rho_liq * M_PI),1./3.); 
							rrToCm0 -= .02 * rrToCm0;		// decrease the radius by 2%
							break;
			}
		} else if (nDims == 2) {
			switch (filled_init) {
				case 0:	rrToCm0 = sqrt(nMol / (rho_liq * M_PI * L.y));
							rrToCm0 += .02 * dx;				// increase the radius by 2%
							break;
				case 1:	rrToCm0 = sqrt(2. * nMol / (rho_liq * M_PI * L.y)); 
							rrToCm0 -= .02 * rrToCm0;		// decrease the radius by 2%
							break;
			}
		}

		/* put drop's CoM into the xy-plane's middle and adjust its z-position */
		switch (filled_init) {
			case 0:	V_SET(cm0, .5 * subGridX, .5 * L.y, rrToCm0 + dz);
						break;
			case 1:	V_SET(cm0, .5 * subGridX, .5 * L.y, dz);
						break;
		}

		rrToCm0 *= rrToCm0;				// make a square out of drop's radius

		printf ("rrToCm0 is %8.4f\n", rrToCm0);
		printf ("cm0 is %8.4f %8.4f %8.4f \n", cm0.x, cm0.y, cm0.z);
		
		/* fill the system with the liquid */
		for (i1 = 1; i1 < grid.x+1; i1++){
			for (j1 = 1; j1 < grid.y+1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create initial liquid field */
					if (NVT == 1) {
						if (nDims == 3 && 
								(SQR((globStartIdx + i1 - 0.5)*dx - cm0.x) +
								 SQR((j1 - 0.5)*dy - cm0.y) + 
								 SQR((k1 - 0.5 - 1.)*dz - corr.z - cm0.z) < rrToCm0) &&
							  	(filled_x_grids < nMol / (rho_liq * dx * dy * dz))) {
							wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
							total_N = total_N + 1.;
							++filled_x_grids;
#warning: for NVT the total_N and filled_x_grids should be known to other CPUs
						} else if (nDims == 2 &&
								(SQR((globStartIdx + i1 - 0.5)*dx-cm0.x) +
								 SQR((k1 - 0.5 - 1.)*dz - corr.z - cm0.z) < rrToCm0) && 
								(filled_x_grids < nMol / (rho_liq * dx * dy * dz))) {
							wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
							total_N = total_N + 1.;
							++filled_x_grids;
#warning: for NVT the total_N and filled_x_grids should be known to other CPUs
						} else {
						}
					} else {
						// grand canonical ensemble
						wa[i1][j1][k1] = 0.;
						// if initial liquid state is the one with the flat interface
						if (filled_init == 0) {
							if ( (k1 - 0.5) * dz > 1.4 * corr.z){
								wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
								total_N = total_N + 1.;
							}
						} else {
						// if initial liquid state is the filled one
							if ( (i1 - 1) * dx < corr.x && (j1 - 1) * dy < corr.y && (k1 - 0.5) * dz < corr.z ) {
							} else {
								wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
								total_N = total_N + 1.;
							}
						}
					}
				}
			}
		}
		printf ("filled_x_grids %d\n", filled_x_grids);
		total_N = total_N * rho_liq * dx * dy * dz;
		printf ("total_N is %8.4f\n", total_N);
		
		printf ("Made the corrugations. Now the substrate has %d nodes\n", nSub);

		VecI currNode;
		V_ZERO (currNode);
		subDepth = (int) ((corr.z + rCut) / dz) + 1;
		
		/* assign the potential provided by the corrugations */
		for (i1 = 1; i1 < grid.x+1; i1++){
			for (j1 = 1; j1 < grid.y+1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create substrate potential */
					sub[i1][j1][k1] = subPotZ[k1]; // set 1d potential of a flat substrate
					if ( (globStartIdx + i1 - 1) * dx < corr.x && (j1 - 1) * dy < corr.y && (k1 - 0.5) * dz < corr.z ) {
						// if inside a corrugation, set large potential:
						sub[i1][j1][k1] =	maxSubPot;
					} else if (k1 >= subDepth) {
					} else {
						// if the node is accesible to fluid, calculate sub.potential there
						for(sn = 0; sn < nSub; sn++){
							V_SET (currNode, (int)((nRep.x-1) / 2) * subGridX + globStartIdx + i1, (int)((nRep.y-1) / 2) * grid.y + j1, k1);
							V_SET (dist2, SQR(subNode[sn].rs.x - currNode.x), SQR(subNode[sn].rs.y - currNode.y), SQR(subNode[sn].rs.z - currNode.z));
							distance2 = dist2.x * dx2 + dist2.y * dy2 + dist2.z * dz2;
							if (distance2 < rrCut) {
								ri2 = sigma_sub2/distance2;
								ri6 = CUBE(ri2);
								sub[i1][j1][k1] += ri6 * (ri6 - 1.);
							}
						}
					}
					sub[i1][j1][k1] *=	4. * eps;

					if (isinf(sub[i1][j1][k1]) == 1 || isnan(sub[i1][j1][k1]) == 1) {
						printf("Error when calculating substrate potential at i=[%4d], j=[%4d], k=[%4d]\n", i1, j1, k1);
						exit(EXIT_FAILURE);
					}
				}
			}
			printf ("Finished subsrtrate creation for i1 = %d \n", i1);
		}
		
		PrintSubFile ();			// write down the sub potential into the output file

		printf ("Finished with initial fields creation! total_N is %10.8f\n",total_N);
		
	} else {	// if we are restarting an old calculation
		ReadSubFile ();
#warning: have to check Wfields work at restart when its IO is accomplished
		//PrintWfieldsToCompare ();
	}
	
	// main iteration loop
	for (kk = 0; kk < iterations; kk++){
		if (kk % TWRITE != 0) {
			write_snapshot = 0;
		} else {
			write_snapshot = 1;
		}
		
		if (kk != 0) {
			for(i1 = 1; i1 < grid.x+1; i1++){
				for(j1 = 1; j1 < grid.y+1; j1++){
					for(k1 = 1; k1 < grid.z+1; k1++){
						wa[i1][j1][k1] = (1.-lambda) * wa[i1][j1][k1] + lambda * waout[i1][j1][k1];
					}
				}
			}

			if ((int)(kk*10) % TWRITE == 0) {
				CheckConverge (kk);
			}

		} else {
			if (restart == 1) write_snapshot = 0;
			MAKE_FILENAME(fullname_conv,"converge_test.dat");
		}
		
		CalcPartSumQ (NVT);
		
		/* density calculation */
		// inside the simulation region
		if (NVT == 1){
			for(i1 = 1; i1 < grid.x+1; i1++){
				for(j1 = 1; j1 < grid.y+1; j1++){
					for(k1 = 1; k1 < grid.z+1; k1++){
						rho[i1][j1][k1] = total_N * exp(- wa[i1][j1][k1]) / Q;
					}
				}
			}
		}	else {
			for(i1 = 1; i1 < grid.x+1; i1++){
				for(j1 = 1; j1 < grid.y+1; j1++){
					for(k1 = 1; k1 < grid.z+1; k1++){
						rho[i1][j1][k1] = rho_liq * exp(V11*rho_liq + W111*SQR(rho_liq) - wa[i1][j1][k1]);
					}
				}
			}
		}
		
		/* account for periodic boundaries
		 c
		 c In order to implement easier
		 c the periodic boundary conditions
		 c we have the
		 c   q(0,..,..) and the q(grid.x+1,..,..) cells
		 c   q(..,0,..) and the q(..,grid.y+1,..) cells
		 c   q(..,..,0) and the q(..,..,grid.z+1) cells
		 c
		 */
		
		for(i1 = 1; i1 < grid.x+1; i1++){
			// z-direction
			for(j1 = 1; j1 < grid.y+1; j1++){
				rho[i1][j1][0] = 0.; //rho[i1][j1][grid.z+1] = 0.;
				// an attempt to implement reflecting boundary
				rho[i1][j1][grid.z+1] = rho[i1][j1][grid.z-1];
			}
			// x-direction
			for(k1 = 1; k1 < grid.z+1; k1++){
				rho[i1][0][k1] = rho[i1][grid.y][k1]; rho[i1][grid.y+1][k1] = rho[i1][1][k1];
			}
		}
		
/*		for(k1 = 1; k1 < grid.z+1; k1++){
			// x-direction
//			for(i1 = 1; i1 < grid.x+1; i1++){
//				rho[i1][0][k1] = rho[i1][grid.y][k1]; rho[i1][grid.y+1][k1] = rho[i1][1][k1];
//			}
			// y-direction
			for(j1 = 1; j1 < grid.y+1; j1++){
				rho[0][j1][k1] = rho[grid.x][j1][k1]; rho[grid.x+1][j1][k1] = rho[1][j1][k1];
			}
		}
*/
		CommHaloYZ ();
		
		/* FINISHED density calculation */
		
		PrintSnapField (0, kk);
		/* finding new fields */
		for(i1 = 1; i1 < grid.x+1; i1++){
			for(j1 = 1; j1 < grid.y+1; j1++){
				for(k1 = 1; k1 < grid.z+1; k1++){
					df_drho[i1][j1][k1] = rho[i1][j1][k1] * (V11 + W111 * rho[i1][j1][k1]);

					grad[i1][j1][k1] = -KAPPA * ( ((rho[i1+1][j1][k1] - 2.*rho[i1][j1][k1] + rho[i1-1][j1][k1])/dx2) +
												 ((rho[i1][j1+1][k1] - 2.*rho[i1][j1][k1] + rho[i1][j1-1][k1])/dy2));
					if (k1 == 1) {
						grad[i1][j1][k1] -= KAPPA * ( (rho[i1][j1][k1+1] - 3.*rho[i1][j1][k1] + 2.*rho[i1][j1][k1-1]) / dz2 );
//					} else if (k1 == grid.z) {
//						grad[i1][j1][k1] -= KAPPA * ( (2.*rho[i1][j1][k1+1] - 3.*rho[i1][j1][k1] + rho[i1][j1][k1-1]) / dz2 );
					} else {
						grad[i1][j1][k1] -= KAPPA * ( (rho[i1][j1][k1+1] - 2.*rho[i1][j1][k1] + rho[i1][j1][k1-1]) / dz2 );
					}

					waout[i1][j1][k1] = df_drho[i1][j1][k1] + grad[i1][j1][k1] + sub[i1][j1][k1];

					
				}
			}
		}
		PrintSnapField (1, kk);

		if (kk % TWRITE == 0) {
			if (write_snapshot == 1 && myRank == 0) {
				// evaluates pressure and calculates free energy
				FILE *wPressure;
				
				// save the pressure at the node [2][5][5]
				wPressure = fopen(fullname_press,"a");
				fprintf (wPressure,	"step %d. The gas pressure is %6.10f \n",
								 kk + iteration_num,
								 kbT * 10. * (rho[2][5][5] +
															(0.5 * V11 * SQR(rho[2][5][5])) +
															(2. * W111 * CUBE(rho[2][5][5])/3.)) );
				fclose(wPressure);
			} else {
			}
			
			// calculate free energy
			CalcFreeEn(iteration_num, kk + iteration_num, NVT);
		}
	}  /* FINISHED SCF iteration loop */
	
	/* Shut down MPI */
	MPI_Finalize();
	
	return(0);
} /* FINISHED main function */

/* calculation of the Q coefficient!!! */
void CalcPartSumQ (int _NVT) {
	int i, j, k;
	MPI_Status status;
	double *sbuf=NULL, *rbuf=NULL;
	
	sbuf = (double*) malloc(2*sizeof(double));
	
	sbuf[0] = sbuf[1] = Q = 0.;
	
	if (_NVT == 1) {
		for(i = 1; i < grid.x+1; i++){
			XYa[i] = 0.;
			for(j = 1; j < grid.y+1; j++){
				YYa[j] = 0.;
				for(k = 1; k < grid.z+1; k++){
					YYa[j] = YYa[j] + koeff_z_open[k]*exp(-wa[i][j][k])*dz;
				}
				XYa[i] = XYa[i] + koeff_y_semi[j]*YYa[j]*dy;
			}
			Q = Q + koeff_x_semi[i]*XYa[i]*dx;
		}
		total_N = total_N / MPIsize;		// just needed to do this. Re-summation after MPI_Allgather
	} else {
		total_N = 0.;
		
		for(i = 1; i < grid.x+1; i++){
			XYa[i]  = 0.;
			XYat[i] = 0.;
			for(j = 1; j < grid.y+1; j++){
				YYa[j]  = 0.;
				YYat[j] = 0.;
				for(k = 1; k < grid.z+1; k++){
					YYa[j]  = YYa[j]  + koeff_z_open[k]*exp(-wa[i][j][k])*dz;
					YYat[j] = YYat[j] + koeff_z_open[k]*rho[i][j][k]*dz;
				}
				XYa[i]  = XYa[i]  + koeff_y_semi[j]*YYa[j]*dy;
				XYat[i] = XYat[i] + koeff_y_semi[j]*YYat[j]*dy;
			}
			Q = Q + koeff_x_semi[i]*XYa[i]*dx;
			total_N = total_N + koeff_x_semi[i]*XYat[i]*dx;
		}
		sbuf[1]=total_N;
	}
	sbuf[0]=Q;
	
	rbuf = (double*) malloc(MPIsize*2*sizeof(double));
	
	MPI_Allgather (sbuf, 2, MPI_DOUBLE, rbuf, 2, MPI_DOUBLE, MPI_COMM_WORLD);
	
	double _Q = 0., _total_N = 0.;
	for (int i = 0; i < MPIsize; i++) {
		_Q += rbuf[2*i];
		_total_N += rbuf[2*i+1];
//		if (myRank == 0) {
//			printf ("Q is %18.14f and total_N is %18.14f\n", rbuf[2*i], rbuf[2*i+1]);
//		}
	}
	Q = _Q;
	total_N = _total_N;
	
	// free buffers
	free(rbuf);
	free(sbuf);
	
//	printf("finished CalcPartSumQ. Value of Q is %8.4f, total_N is %8.4f\n", Q, total_N);
}

void CalcFreeEn (int temp_iteration_num, int temp_kk, int _NVT) {
	V_ZERO (rho_fd);
	V_ZERO (rho_bd);
	V_ZERO (rho_cd);
	
	int i, j, k;
	
	double  term1,term2,term3,term4,term_sub;
	double 	sum_final;
	extern double ***sub;
	
	FILE *free_energy_out;
	MPI_Status status;
	double *sbuf=NULL, *rbuf=NULL;
	
	sbuf = (double*) malloc(sizeof(double));
	
	if (myRank == 0) {
		sprintf(name,"free_energy_%d.dat",temp_iteration_num);
		MAKE_FILENAME(fullname_free_en,name);
		free_energy_out = fopen(fullname_free_en,"a");
		if (temp_kk-iteration_num == 0) {
			fprintf (free_energy_out, "#   step          Q     FreeEn   total_N\n");
		} else {
			// do nothing
		}
	}
	
	if (_NVT == 1) {
		term1 = -total_N * (log(Q / total_N) + 1.);
	} else {
		term1 =	-rho_liq * Q * exp(rho_liq * (V11 + W111 * rho_liq));
	}
	term2 = 0.;
	term3 = 0.;
	term4 = 0.;
	term_sub = 0.;
	
	sum_final = 0.;

	/* calculation of the free energy*/
	for(i = 1; i < grid.x+1; i++){
		XYa[i] = 0.;
		for(j = 1; j < grid.y+1; j++){
			YYa[j] = 0.;
			for(k = 1; k < grid.z+1; k++){
				/* central differences * 2. in 3 dimensions */
				rho_cd.x = (rho[i+1][j][k] - rho[i-1][j][k])/dx;
				rho_cd.y = (rho[i][j+1][k] - rho[i][j-1][k])/dy;
				
				if (k == 1) {
					rho_cd.z = (rho[i][j][k+1] + rho[i][j][k] - 2.*rho[i][j][k-1])/dz;
					//				} else if (k == grid.z) {
					//					rho_cd.z = (2.*rho[i][j][k+1] - rho[i][j][k] - rho[i][j][k-1])/dz;
				} else {
					rho_cd.z = (rho[i][j][k+1] - rho[i][j][k-1])/dz;
				}
				
				term2 = V11 * 0.5 * SQR(rho[i][j][k]) +	W111 * (1./3.) * CUBE(rho[i][j][k]);
				term3 = 0.125 * KAPPA * (SQR(rho_cd.x) + SQR(rho_cd.y) + SQR(rho_cd.z));
				term4 = -wa[i][j][k] * rho[i][j][k];
				term_sub = sub[i][j][k] * rho[i][j][k];
				
				YYa[j] = YYa[j] +  koeff_z_open[k]*(term2 + term3 + term4 + term_sub)*dz;

			}
			XYa[i] = XYa[i] + koeff_y_semi[j]*YYa[j]*dy;
		}
		sum_final  = sum_final  + koeff_x_semi[i]*XYa[i]*dx;
	}
	sbuf[0]=sum_final;
	
	rbuf = (double*) malloc(MPIsize*sizeof(double));
	
	MPI_Gather (sbuf, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	double _sum_final = 0.;
	for (int i = 0; i < MPIsize; i++) {
		_sum_final += rbuf[i];
	}
	
	if (myRank == 0) {
		sum_final =  _sum_final + term1;

		/* print into the file */
		fprintf(free_energy_out,"%d %10.8f %10.8f %10.8f\n", temp_kk, Q, sum_final * kbT, total_N);
		
		/* print onto the screen */
		printf("%d %10.8f %10.8f %g %10.8f\n", temp_kk, sum_final * kbT, Q, convergence, total_N);
		
		fclose(free_energy_out);
	}
	
	// free buffers
	free(rbuf);
	free(sbuf);
}

void DefineSimpson () {
	int i;
	/* ///////////////////////////////////////
	 c      Define the coefficients for the
	 c    Simpson OPEN formula of integration
	 c/////////////////////////////////////// */
	
	/*     koeff_x_open[1] 	  = koeff_y_open[1]    = koeff_z_open[1] 	= 25./24.;
	 koeff_x_open[2] 	  = koeff_y_open[2]    = koeff_z_open[2] 	= 25./24.;
	 koeff_x_open[3]	  = koeff_y_open[3]    = koeff_z_open[3] 	= 3./4.;
	 koeff_x_open[grid.x-2] = koeff_y_open[grid.y-2] = koeff_z_open[grid.z-2] = 3./4.;
	 koeff_x_open[grid.x-1] = koeff_y_open[grid.y-1] = koeff_z_open[grid.z-1] = 25./24.;
	 koeff_x_open[grid.x]   = koeff_y_open[grid.y]   = koeff_z_open[grid.z] 	= 25./24.;
	 
	 for(i = 4; i <= grid.x - 3; i++){
	 if (i % 2 == 0) koeff_x_open[i] = 4./3.;
	 else koeff_x_open[i] = 2./3.;
	 }
	 
	 for(i = 4; i <= grid.y - 3; i++){
	 if (i % 2 == 0) koeff_y_open[i] = 4./3.;
	 else koeff_y_open[i] = 2./3.;
	 }
	 
	 for(i = 4; i <= grid.z - 3; i++){
	 if (i % 2 == 0) koeff_z_open[i] = 4./3.;
	 else koeff_z_open[i] = 2./3.;
	 }
	 */
	
	/* now we use semi-open formulae in x and y directions
	 * with 100% they preserve a correct volume through
	 * periodic boundary conditions in this directions
	 */
	if (myRank == 0) {
		koeff_x_semi[1] = 3./8.;
		koeff_x_semi[2] = 7./6.;
		koeff_x_semi[3] = 23./24.;
	} else {
		koeff_x_semi[1] = koeff_x_semi[2] = koeff_x_semi[3] = 1.;
	}
	
	koeff_y_semi[1] = 3./8.;
	koeff_y_semi[2] = 7./6.;
	koeff_y_semi[3] = 23./24.;
	
	if (myRank == MPIsize - 1) {
		koeff_x_semi[grid.x-2] = 55./24.;
		koeff_x_semi[grid.x-1] = -1./6.;
		koeff_x_semi[grid.x]	= 11./8.;
	} else {
		koeff_x_semi[grid.x-2] = koeff_x_semi[grid.x-1] = koeff_x_semi[grid.x] = 1.;
	}
	
	koeff_y_semi[grid.y-2]	= 55./24.;
	koeff_y_semi[grid.y-1]	= -1./6.;
	koeff_y_semi[grid.y]		= 11./8.;
	/*	koeff_x_close[grid.x-2]	= koeff_y_close[grid.y-2]	= 23./24.;
	 koeff_x_close[grid.x-1]	= koeff_y_close[grid.y-1]	= 7./6.;
	 koeff_x_close[grid.x]	= koeff_y_close[grid.y]		= 3./8.;
	 */
	
	/* in z-direction we compose our own formulae due to
	 * the location of the fixed values 1/2 step outside the
	 * computational domain (the corresponded derivatives had
	 * also been changed)
	 */
	koeff_z_open[1] = koeff_z_open[grid.z]		= 11./ 24.;
	koeff_z_open[2] = koeff_z_open[grid.z-1]	= 7./ 8.;
	koeff_z_open[3] = koeff_z_open[grid.z-2]	= 5./ 3.;
	
	for(i = 4; i <= grid.x - 3; i++){
		koeff_x_semi[i] = 1.;
	}
	
	for(i = 4; i <= grid.y - 3; i++){
		koeff_y_semi[i] = 1.;
	}
	
	for(i = 4; i <= grid.z - 3; i++){
		koeff_z_open[i] = 1.;
	}
}

void AllocArrays () {
	
	ALLOC_MEM (subPotZ, grid.z+1, double);

	ALLOC_MEM (koeff_x_semi, grid.x+1, double);
	ALLOC_MEM (koeff_y_semi, grid.y+1, double);
	ALLOC_MEM (koeff_z_open, grid.z+1, double);
	
	ALLOC_MEM (XYa, grid.x+1, double);
	ALLOC_MEM (YYa, grid.y+1, double);
	ALLOC_MEM (XYat, grid.x+1, double);
	ALLOC_MEM (YYat, grid.y+1, double);
	
	ALLOC_MEM3 (wa, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sub, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (waout, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (rho, grid.x+2, grid.y+2, grid.z+2, double);
	
	ALLOC_MEM3 (df_drho, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (grad, grid.x+1, grid.y+1, grid.z+1, double);
}

void AllocSubstrate () {
	ALLOC_MEM (subNode, (int)(nRep.x * subGridX * nRep.y * grid.y * ((int)(corr.z / dz) + 1)), substrateNode);
}

void DivideSpace () {
	subGridX = grid.x;
	if (MPIsize != 0) {
		globStartIdx = (int) (floor(myRank * subGridX / MPIsize));
		grid.x = (int) (floor((myRank + 1) * subGridX / MPIsize) - floor(myRank * subGridX / MPIsize));
		printf ("Proc %d: my grid.x is %d, my index starts with %d in global sence\n", myRank, grid.x, globStartIdx);
	} else {
//		do nothing
	}
}

/* Check for convergence */
void CheckConverge (int _kk) {
	FILE *converge;
	MPI_Status status;
	double *sbuf=NULL, *rbuf=NULL;
	
	sbuf = (double*) malloc(2*sizeof(double));
	
	double W2A, tempA;
	W2A = tempA = 0.;
	
	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				W2A = W2A +  wa[i][j][k] * wa[i][j][k];
				tempA = tempA + SQR(wa[i][j][k] - waout[i][j][k]);
			}
		}
	}
	sbuf[0]=W2A; sbuf[1]=tempA;
	
	if (myRank == 0) {
		rbuf = (double*) malloc(MPIsize*2*sizeof(double));
	}

	MPI_Gather (sbuf, 2, MPI_DOUBLE, rbuf, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myRank == 0) {
		W2A = tempA = 0.;			// we will start collecting their values from rbuf
		for (int i = 0; i < MPIsize; i++) {
			W2A += rbuf[2*i];
			tempA += rbuf[2*i+1];
		}
		
		convergence = tempA/W2A;
		printf ("convergence is %gf\n", convergence);

//		if((int)(_kk*10) % TWRITE == 0){
			converge = fopen(fullname_conv,"a");
			fprintf(converge, "%d %g \n", _kk + iteration_num, convergence);
			fclose(converge);
//		}
		// free buffers
		free(rbuf);
	}
	
	// free buffers
	free(sbuf);
}

/* Communicate YZ planes between processors */
void CommHaloYZ () {
	int i, j, k, index;					// running indices and index of the node to be copied
	int count;								// number of data to transfer
	int rnode, snode;						// rank of the node to receive from and to send to
	double *sbuf=NULL, *rbuf=NULL;
	MPI_Status status;
	
	//////////////////////
	//// X-direction /////
	//////////////////////
	
	count = grid.y * grid.z;
	sbuf = (double*) malloc(count*sizeof(double));
	rbuf = (double*) malloc(count*sizeof(double));

	/* send to right, recv from left */
	snode = myRank + 1;
	rnode = myRank - 1;
	if (myRank == 0) rnode = MPIsize - 1;
	if (myRank == MPIsize - 1) snode = 0;
	
	// prepare message for sending
	i = grid.x;
	for (k = 1; k < grid.z + 1; k++) {
		for (j = 1; j < grid.y + 1; j++) {
			index = grid.y*(k-1) + (j-1);		// rightmost yz-plane with real data
			sbuf[index] = rho[i][j][k];
		}
	}
	
	// send and receive data or use memcpy if number of CPU in x-dir is 1
	if (MPIsize > 1) {
//		MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, MPI_COMM_WORLD, rbuf, count, MPI_DOUBLE, rnode, MPI_COMM_WORLD, comm_cart, &status);
		if (myRank % 2 == 0) {
			MPI_Send (sbuf, count, MPI_DOUBLE, snode, 302, MPI_COMM_WORLD);
			MPI_Recv (rbuf, count, MPI_DOUBLE, rnode, 301, MPI_COMM_WORLD, &status);
		} else {
			MPI_Recv (rbuf, count, MPI_DOUBLE, rnode, 302, MPI_COMM_WORLD, &status);
			MPI_Send (sbuf, count, MPI_DOUBLE, snode, 301, MPI_COMM_WORLD);
		}
	} else {
		memcpy(rbuf,sbuf,count*sizeof(double));
	}
	
	// unpack message
	i = 0;
	for (k = 1; k < grid.z + 1; k++) {
		for (j = 1; j < grid.y + 1; j++) {
			index = grid.y*(k-1) + (j-1);		// rightmost yz-plane with real data
			rho[i][j][k] = rbuf[index];
		}
	}
	
	/* send to left, recv from right i = 2, 8, 10, 12, 14 */
	snode = myRank - 1;
	rnode = myRank + 1;
	if (myRank == 0) snode = MPIsize - 1;
	if (myRank == MPIsize - 1) rnode = 0;
	
	// prepare message for sending
	i = 1;
	for (k = 1; k < grid.z + 1; k++) {
		for (j = 1; j < grid.y + 1; j++) {
			index = grid.y*(k-1) + (j-1);		// rightmost yz-plane with real data
			sbuf[index] = rho[i][j][k];
		}
	}
	
	// send and receive data or use memcpy if number of CPU in x-dir is 1
	if (MPIsize > 1) {
		//		MPI_Sendrecv(sbuf, count, MPI_DOUBLE, snode, MPI_COMM_WORLD, rbuf, count, MPI_DOUBLE, rnode, MPI_COMM_WORLD, comm_cart, &status);
		if (myRank % 2 == 0) {
			MPI_Send (sbuf, count, MPI_DOUBLE, snode, 402, MPI_COMM_WORLD);
			MPI_Recv (rbuf, count, MPI_DOUBLE, rnode, 401, MPI_COMM_WORLD, &status);
		} else {
			MPI_Recv (rbuf, count, MPI_DOUBLE, rnode, 402, MPI_COMM_WORLD, &status);
			MPI_Send (sbuf, count, MPI_DOUBLE, snode, 401, MPI_COMM_WORLD);
		}
	} else {
		memcpy(rbuf,sbuf,count*sizeof(double));
	}
	
	// unpack message
	i = grid.x+1;
	for (k = 1; k < grid.z + 1; k++) {
		for (j = 1; j < grid.y + 1; j++) {
			index = grid.y*(k-1) + (j-1);		// rightmost yz-plane with real data
			rho[i][j][k] = rbuf[index];
		}
	}
	
	// free buffers
	free(rbuf);
	free(sbuf);
}
