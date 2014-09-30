#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <string.h>

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
char name[100];
char foldername[100];
char fullname_iter[100];
char fullname_Wfields[100];
char fullname_Wsub[100];
char fullname_free_en[100];
char fullname_conv[100];
char fullname_snap[100];
char fullname_param[100];
char fullname_press[100];
char fullname_test[100];
char path[50];
char folder[50];
double write_snapshot;

/* global properties of the system*/
double kbT;					// kbT :-)
double rho_liq;			// density in the bulk (at the top of the box)
double lambda;				// relaxation parameter of the SC-scheme

/* properties of the run */
int	restart;					// restart flag
int	iteration_num;		// if restarted, the last existing iteration
int	iterations;				// desired iterations limit
double convergence;	// convergence

/* size and grid related properties */
VecR L;								// box length
VecI nRep;						// number of replicas of the initial box
VecI grid;						// number of grid points
VecR corr;						// dimensions of the corrugation
VecR cav;							// cavity dimensions (z-one is irrelevant)
VecR rho_fd;					// forward difference in density
VecR rho_bd;					// backward difference in density
VecR rho_cd;					// central difference in density
double dx,dy,dz;			// size of the grid cell
double dx2, dy2, dz2;	// size of the grid cell squared
double volume;				// volume of the computation domain
double dV;						// volume of the elementary cell
double nMol;					// number of water molecules for NVT-ensemble
int nDims;						// number of dimensions
int NVT;							// canonical ensemble

/* substrate' parameters */
double rCut;					// cut off distance for potential
double rrCut;					// squared cut off distance for potential
double sigma_sub;			// sigma of the substrate potential
double sigma_sub2;		// sigma of the substrate potential squared
double eps;						// amplitude of the potential

/* statistical properties of the ensemble */
double total_N;				// total number of molecules in the system
double Q;							// statistical sum
double check_N;

/* time variables */
time_t now, later;		// two marks of the time-measurement
double seconds;				// length of time interval between 2 measurements

/* integration variables */
double *koeff_x_semi;
double *koeff_y_semi;
double *koeff_z_open;
double *XYa,*YYa;
double *XYat,*YYat;

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
 c
 c       |<------------- L.x ------------->|
 c       |                                 |
 c   0__ |__1__2__.....__grid.x-1__grid.x__|__grid.x+1
 c       |                                 |
 c       |                                 |
 c       |
 c
 c       0 = grid.x, grid.x+1 = 1   _____ = L.x / grid.x
 c
 c
 c    Make sure since the integrations will be performed by
 c    Simpson method grid.x grid.y grid.z must be odd numbers
 c
 */

int main (int argc, char **argv){
	time(&now);
	
	PassConsoleParams (argc, argv);

	double W2A, tempA;
	double r_sub, ri2, ri6;
	double distance2;
	VecR dist2;
	int i1,j1,k1,kk;
	int is, js, ks, sn;
	int corrShift;
	int subDepth;
	VecR cm0;
	double rrToCm0;
	
	AllocArrays ();
	
	FILE  *Wtest, *Wrestart, *Wfields, *Wsub,
				*snapshot, *iterkeeper, *converge, *wPressure;
		
	/* Lengths are measured in nm */
	iterations = 6000000;
	kbT = 4.116404397;
	rCut = 5.;
	sigma_sub = .3;
	sigma_sub2 = SQR(sigma_sub);

	rrCut = SQR(rCut);
	strcpy (foldername, path);
	strcat (foldername, folder);
    
	InitParameters ();
	StoreParameters ();
	
	if (nDims == 3) {
		rrToCm0 = pow(.75 * nMol / (rho_liq * M_PI),1./3.);
	} else if (nDims == 2) {
		rrToCm0 = sqrt(nMol / (rho_liq * M_PI * L.y));
	}
	rrToCm0 += .02 * dx;
	V_SET(cm0, .5 * (L.x - cav.x), .5 * L.y, rrToCm0 + dz);
	rrToCm0 *= rrToCm0;
	dV = dx*dy*dz;

	printf ("rrToCm0 is %8.4f\n", rrToCm0);
	printf ("cm0 is %8.4f %8.4f %8.4f \n", cm0.x, cm0.y, cm0.z);
	
	DefineSimpson ();
	
	/* set initial values of field, gradients, etc. to zero. 
	 Create nodes of flat substrate */
	nSub = 0;
	corrShift = (int) (corr.z / dz) + 1;
	subDepth = (int) (rCut / dz);
	V_SET (nRep, (int) (2 * ((int) (rCut / L.x) + 1) + 1),
							 (int) (2 * ((int) (rCut / L.y) + 1) + 1),
							 0);
	// total number of box replicas in y-direction, inc. the original one

	printf ("total number of box replicas in x is %8d, in y is %8d\n", nRep.x, nRep.y);
	
	// set all the fields to zero
	for (i1 = 1; i1 < grid.x+1; i1++){
		for (j1 = 1; j1 < grid.y+1; j1++){
			for (k1 = 1; k1 < grid.z+1; k1++){
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
	
	if (restart == 0) {
		// set up the flat substrate
		for (i1 = 1; i1 < nRep.x * grid.x + 1; i1++){
			for (j1 = 1; j1 < nRep.y * grid.y + 1; j1++){
				// create "0"th substrate plane
				V_SET(subNode[nSub].rs, i1, j1, 0);
				++nSub;
//			printf ("nSub is %8d\n", nSub);
			}
		}
		printf ("set up the flat substrate\n");
	
		/* calculate z-dependence of the potential created by a flat substrate */
		for (k1 = 1; k1 < grid.z+1; k1++){
			for(sn = 0; sn < nSub; sn++){
				dist2.x = SQR(subNode[sn].rs.x - ((int) ((nRep.x - 1) / 2) + 0.5) * grid.x - 1);
//			dist2.x = MIN(dist2.x, SQR(subNode[sn].rs.x - 1 + grid.x));
//			dist2.x = MIN(dist2.x, SQR(subNode[sn].rs.x - 1 - grid.x));
			
				dist2.y = SQR(subNode[sn].rs.y - ((int) ((nRep.y - 1) / 2) + 0.5) * grid.y - 1);
//			dist2.y = MIN(dist2.y, SQR(subNode[sn].rs.y - 1 + grid.y));
//			dist2.y = MIN(dist2.y, SQR(subNode[sn].rs.y - 1 - grid.y));
			
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
			subPotZ[k1] *=	4. * eps;
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
		
		// create corrugation sites (inkl. replicas of the box)
		for (i1 = 1; i1 < nRep.x * grid.x + 1; i1++){
			for (j1 = 1; j1 < nRep.y * grid.y + 1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create substrate nodes taking into account corrugation */
					if ((((i1 - 1) % (int)((corr.x + cav.x) / dx)) * dx < corr.x) &&
						 (((j1 - 1) % (int)((corr.y + cav.y) / dy)) * dy < corr.y) &&
						 ((k1 - .5)*dz < corr.z) ) {
						V_SET(subNode[nSub].rs, i1, j1, k1);
						++nSub;
					} else {
					}
				}
			}
		}
		
		// create liquid
		for (i1 = 1; i1 < grid.x+1; i1++){
			for (j1 = 1; j1 < grid.y+1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create initial liquid field */
					if (NVT == 1) {
						if (nDims == 3 && (SQR((i1 - 0.5)*dx-cm0.x)+SQR((j1 - 0.5)*dy-cm0.y)
															 +SQR((k1 - 0.5 - 1.)*dz - corr.z - cm0.z) < rrToCm0)
													 && (filled_x_grids < nMol / (rho_liq * dV))) {
							wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
							total_N = total_N + 1.;
							++filled_x_grids;
						} else if	(nDims == 2 &&
											 (SQR((i1 - 0.5)*dx-cm0.x) +
												SQR((k1 - 0.5 - 1.)*dz - corr.z - cm0.z) < rrToCm0) &&
											 (filled_x_grids < nMol / (rho_liq * dV))) {
							wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
							total_N = total_N + 1.;
							++filled_x_grids;
						}
					} else {
						wa[i1][j1][k1] = 0.;
//						if ((0.5 * dz + (k1 - 1)*dz) > 2. * corr.z){
						if ((0.5 * dz + (k1 - 1)*dz) > 2. * corr.z){
							wa[i1][j1][k1] = V11*rho_liq + W111*SQR(rho_liq);
							total_N = total_N + 1.;
						}
					}
				}
			}
		}
		printf ("filled_x_grids %d\n", filled_x_grids);
		total_N = total_N * rho_liq * dx * dy * dz;
		printf ("total_N is %8.4f\n", total_N);
		
		printf ("Made the corrugations. Now the substrate has %d nodes\n", nSub);

		MAKE_FILENAME(fullname_Wsub,"substrate.dat");
		Wsub=fopen(fullname_Wsub,"a");
				
		// assign the potential provided by the corrugations
		for (i1 = 1; i1 < grid.x+1; i1++){
			for (j1 = 1; j1 < grid.y+1; j1++){
				for (k1 = 1; k1 < grid.z+1; k1++){
					/* create substrate potential */
					sub[i1][j1][k1] = subPotZ[k1]; // set 1d potential of a flat substrate

					if ((((i1 - 1) % (int)((corr.x + cav.x) / dx)) * dx < corr.x) &&
						 (((j1 - 1) % (int)((corr.y + cav.y) / dy)) * dy < corr.y) &&
						 ((k1 - .5) * dz < corr.z) ) {
						// if inside a corrugation, set large potential:
						sub[i1][j1][k1] =	150000.;
					} else {
						// if the node is accesible to fluid, calculate sub.potential there
						for(sn = 0; sn < nSub; sn++){
							// account for PBC in x
							dist2.x = SQR(subNode[sn].rs.x - (int) ((nRep.x - 1) / 2) * grid.x - i1);
//							dist2.x = MIN(dist2.x, SQR(subNode[sn].rs.x - i1 + grid.x));
//							dist2.x = MIN(dist2.x, SQR(subNode[sn].rs.x - i1 - grid.x));

							// account for PBC in y
							dist2.y = SQR(subNode[sn].rs.y - (int) ((nRep.y - 1) / 2) * grid.y - j1);
//							dist2.y = MIN(dist2.y, SQR(subNode[sn].rs.y - j1 + grid.y));
//							dist2.y = MIN(dist2.y, SQR(subNode[sn].rs.y - j1 - grid.y));
							
							dist2.z = SQR(subNode[sn].rs.z - k1); // no PBC in z!
		
							distance2 = dist2.x * dx2 + dist2.y * dy2 + dist2.z * dz2;
							if (distance2 < rrCut) {
								ri2 = sigma_sub2/distance2;
								ri6 = CUBE(ri2);
								sub[i1][j1][k1] +=	4. * eps * ri6 * (ri6 - 1.);
							}
						}
					}
					fprintf(Wsub,"%16.12f \n", sub[i1][j1][k1]);
				}
			}
			printf ("Finished subsrtrate creation for i1 = %d \n", i1);
			printf (" dist2.x is %d \n", SQR(- (int) ((nRep.x - 1) / 2) * grid.x - i1));
		}
		fclose(Wsub);
		printf ("Finished with initial fields creation! total_N is %10.8f\n",total_N);
		
	} else {
		// read iteration number and total number of particles for restart
		MAKE_FILENAME(fullname_iter,"iteration.dat");
		iterkeeper = fopen(fullname_iter,"r");
		fscanf(iterkeeper,"%d %lf\n",&iteration_num, &total_N);
		fclose(iterkeeper);
		
		// read filename of the fields for restart
		sprintf(name,"Wfields_%d.dat",iteration_num);
		MAKE_FILENAME(fullname_Wfields,name);
		Wfields=fopen(fullname_Wfields,"r");
		printf ("fullname_Wfields is %s\n", fullname_Wfields);
		
		// read substrate potential for restart
		MAKE_FILENAME(fullname_Wsub,"substrate.dat");
		Wsub = fopen(fullname_Wsub,"r");

		for(i1 = 1; i1 < grid.x+1; i1++){
			for(j1 = 1; j1 < grid.y+1; j1++){
				for(k1 = 1; k1 < grid.z+1; k1++){
					/* set free energy terms to zero*/
					SET_TO_ZERO (df_drho,i1,j1,k1);
					SET_TO_ZERO (grad,i1,j1,k1);
					
					/* read in the field and the substrate potential */
					fscanf(Wfields,"%lf \n",&wa[i1][j1][k1]);
					fscanf(Wsub,"%lf \n",&sub[i1][j1][k1]);
					
				}
			}
		}
		printf ("Finished with scanning substrate file!\n");
		printf ("total_N is %10.8f\n", total_N);

		fclose(Wfields);
	}
	
	for (kk = 0; kk < iterations; kk++){
		
//		if (kk % TWRITE != 0 || (kk == 0 && iteration_num != 0) ) {
		if (kk % TWRITE != 0) {
			write_snapshot = 0;
		} else {
			write_snapshot = 1;
		}
		
		if (kk != 0) {
			/* Set initial conditions and check for convergence */
			W2A = 0.;
			tempA = 0.;
			
			for(i1 = 1; i1 < grid.x+1; i1++){
				for(j1 = 1; j1 < grid.y+1; j1++){
					for(k1 = 1; k1 < grid.z+1; k1++){
						wa[i1][j1][k1] = (1.-lambda) * wa[i1][j1][k1] + lambda * waout[i1][j1][k1];
						W2A = W2A +  wa[i1][j1][k1] * wa[i1][j1][k1];
						tempA = tempA + SQR(wa[i1][j1][k1] - waout[i1][j1][k1]);
					}
				}
			}
			
			convergence = tempA/W2A;
			
			if((int)(kk*10) % TWRITE == 0){
				converge = fopen(fullname_conv,"a");
				fprintf(converge, "%d %g \n", kk + iteration_num, convergence);
				fclose(converge);
			}
		} else {
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
		
		/*
		 c
		 c In order to implement easier
		 c the periodic boundary conditions
		 c we have the
		 c   q(0,..,..) and the q(grid.x+1,..,..) cells
		 c   q(..,0,..) and the q(..,grid.y+1,..) cells
		 c   q(..,..,0) and the q(..,..,grid.z+1) cells
		 c
		 */
		
		/* account for periodic boundaries */
		for(i1 = 1; i1 < grid.x+1; i1++){
			// z-direction
			for(j1 = 1; j1 < grid.y+1; j1++){
				rho[i1][j1][0] = 0.; //rho[i1][j1][grid.z+1] = 0.;
				// an attempt to implement reflecting boundary
				rho[i1][j1][grid.z+1] = rho[i1][j1][grid.z-1];
			}
		}
		
		for(k1 = 1; k1 < grid.z+1; k1++){
			// x-direction
			for(i1 = 1; i1 < grid.x+1; i1++){
				rho[i1][0][k1] = rho[i1][grid.y][k1]; rho[i1][grid.y+1][k1] = rho[i1][1][k1];
			}
			// y-direction
			for(j1 = 1; j1 < grid.y+1; j1++){
				rho[0][j1][k1] = rho[grid.x][j1][k1]; rho[grid.x+1][j1][k1] = rho[1][j1][k1];
			}
		}

		/* FINISHED density calculation */
		
		PrintSnapField (0, kk, 0, 0, 0);
		
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

					PrintSnapField (1, kk, i1, j1, k1);
					
				}
			}
		}
		
		PrintSnapField (2, kk, 0, 0, 0);
	}  /* FINISHED SCF iteration loop */
	
	return(0);
} /* FINISHED main function */

void CalcFreeEn (int temp_iteration_num, int temp_kk, int temp_NVT) {
	V_ZERO (rho_fd);
	V_ZERO (rho_bd);
	V_ZERO (rho_cd);
	
	int i, j, k;
	
	double  term1,term2,term3,term4,term_sub;
	double 	sum_final;
	
	FILE *free_energy_out;
	
	sprintf(name,"free_energy_%d.dat",temp_iteration_num);
	MAKE_FILENAME(fullname_free_en,name);
	free_energy_out = fopen(fullname_free_en,"a");
	if (temp_kk-iteration_num == 0) {
		fprintf (free_energy_out, "#   step          Q     FreeEn   total_N\n");
	} else {
		// do nothing
	}
	
	if (temp_NVT == 1) {
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
	
	sum_final =  sum_final + term1;
	
	/* print into the file */
	fprintf(free_energy_out,"%d %10.8f %10.8f %10.8f\n", temp_kk, Q,
			sum_final * kbT, total_N);
	
	/* print onto the screen */
	printf("%d %10.8f %10.8f %g %10.8f %8.4f\n", temp_kk, sum_final * kbT, Q, convergence, total_N, seconds);
	
	fclose(free_energy_out);
}

void PrintSnapField (int _id, int _kk, int _i1, int _j1, int _k1) {
	FILE *snapshot, *Wrestart;
	int i1, j1, k1;
	
	if (_id == 0 && _kk % TWRITE == 0) {
		FILE *iterkeeper;
		extern double check_N;
		
		/* output into files */
		time(&later);
		seconds = difftime(later,now); now = later;
			
		/* calculate check_N and test it output onto the screen every step */
		check_N = 0.;
		for (i1 = 1; i1 < grid.x+1; i1++){
			XYat[i1] = 0.;
			for (j1 = 1; j1 < grid.y+1; j1++){
				YYat[j1] = 0.;
				for (k1 = 1; k1 < grid.z+1; k1++){
					YYat[j1] = YYat[j1] + koeff_z_open[k1]*rho[i1][j1][k1]*dz;
				}
				XYat[i1] = XYat[i1] + koeff_y_semi[j1]*YYat[j1]*dy;
			}
			check_N = check_N + koeff_x_semi[i1]*XYat[i1]*dx;
		}
		if (write_snapshot == 1) {
			// save current interation number
			MAKE_FILENAME(fullname_iter,"iteration.dat");
			iterkeeper = fopen(fullname_iter,"w");
			fprintf(iterkeeper,"%d %10.8f\n",_kk + iteration_num, total_N);
			fclose(iterkeeper);
				
			// save snapshots (here is a caption ONLY)
			sprintf(name,"snapshot_%d.dat",_kk + iteration_num);
			MAKE_FILENAME(fullname_snap,name);
			snapshot = fopen(fullname_snap,"w");
			fprintf (snapshot,"VARIABLES = \"X\", \"Y\", \"Z\",\"density\",\"df_drho\",\"grad\",\"subpot\",\"wa\",\n");
			fprintf (snapshot,"ZONE I=%d, J=%d, K=%d, F=POINT\n",grid.z, grid.y, grid.x);
				
			// make filename for the fields for future restart
			sprintf(name,"Wfields_%d.dat",_kk + iteration_num);
			MAKE_FILENAME(fullname_Wfields,name);
			Wrestart=fopen(fullname_Wfields,"w");
		}
	} else if (_id == 1 && write_snapshot != 0) {
		// save the snapshot (the rest) and fields for future restart
		fprintf(snapshot,"%g %g %g %g %g %g %g %g\n",
						(_i1 - 0.5)*dx, (_j1 - 0.5)*dy, (_k1 - 0.5)*dz,
						rho[_i1][_j1][_k1], df_drho[_i1][_j1][_k1], grad[_i1][_j1][_k1],
						sub[_i1][_j1][_k1], wa[_i1][_j1][_k1]);
			
		fprintf(Wrestart,"%16.12f \n", wa[_i1][_j1][_k1]);
	} else if (_id == 2 && write_snapshot == 1) {
		FILE *wPressure;
		// close snapshot and field files
		fclose(snapshot); fclose(Wrestart);
			
		// save the gas pressure
		wPressure = fopen(fullname_press,"a");
		fprintf (wPressure,	"step %d. The gas pressure is %6.10f \n",
						 _kk + iteration_num, kbT*10.*(rho[2][5][5] + (0.5 * V11 * SQR(rho[2][5][5])) +
																					(2. * W111 * CUBE(rho[2][5][5])/3.)) );
		fclose(wPressure);
			
		// calculate free energy
		CalcFreeEn(iteration_num, _kk + iteration_num, NVT);
	} else {
	}
}

/* calculation of the Q coefficient!!! */
void CalcPartSumQ (int _NVT) {
	int i1, j1, k1;
	
	Q = 0.;
	
	if (_NVT == 1) {
		for(i1 = 1; i1 < grid.x+1; i1++){
			XYa[i1] = 0.;
			for(j1 = 1; j1 < grid.y+1; j1++){
				YYa[j1] = 0.;
				for(k1 = 1; k1 < grid.z+1; k1++){
					YYa[j1] = YYa[j1] + koeff_z_open[k1]*exp(-wa[i1][j1][k1])*dz;
				}
				XYa[i1] = XYa[i1] + koeff_y_semi[j1]*YYa[j1]*dy;
			}
			Q = Q + koeff_x_semi[i1]*XYa[i1]*dx;
		}
	} else {
		total_N = 0.;
		
		for(i1 = 1; i1 < grid.x+1; i1++){
			XYa[i1]  = 0.;
			XYat[i1] = 0.;
			for(j1 = 1; j1 < grid.y+1; j1++){
				YYa[j1]  = 0.;
				YYat[j1] = 0.;
				for(k1 = 1; k1 < grid.z+1; k1++){
					YYa[j1]  = YYa[j1]  + koeff_z_open[k1]*exp(-wa[i1][j1][k1])*dz;
					YYat[j1] = YYat[j1] + koeff_z_open[k1]*rho[i1][j1][k1]*dz;
				}
				XYa[i1]  = XYa[i1]  + koeff_y_semi[j1]*YYa[j1]*dy;
				XYat[i1] = XYat[i1] + koeff_y_semi[j1]*YYat[j1]*dy;
			}
			Q = Q + koeff_x_semi[i1]*XYa[i1]*dx;
			total_N = total_N + koeff_x_semi[i1]*XYat[i1]*dx;
		}
	}
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
	koeff_x_semi[1] 	  = koeff_y_semi[1]    = 3./8.;
	koeff_x_semi[2] 	  = koeff_y_semi[2]    = 7./6.;
	koeff_x_semi[3] 	  = koeff_y_semi[3]    = 23./24.;
	
	koeff_x_semi[grid.x-2]	= koeff_y_semi[grid.y-2]	= 55./24.;
	koeff_x_semi[grid.x-1]	= koeff_y_semi[grid.y-1]	= -1./6.;
	koeff_x_semi[grid.x]		= koeff_y_semi[grid.y]		= 11./8.;
	
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
	ALLOC_MEM (subNode, (int)(nRep.x * grid.x * nRep.y * grid.y * ((int)(corr.z / dz) + 1)), substrateNode);
}
