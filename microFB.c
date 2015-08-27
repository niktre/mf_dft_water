#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "in_vdefs.h"

char name[120];
char foldername[120];
char path[80];
char folder[80];
char fullname_snap[120];
char fullname_stress[120];
char fullname_Wfields[120];
char fullname_Wsub[120];
char fullname_diff[120];

#define	V11			-1.086641506142464
#define	W111		0.023102120829070
#define	KAPPA		2 * 0.015739171

/* integration variables */
double *koeff_x_semi, *koeff_y_semi, *koeff_z_open;
double *XYa,*YYa;
double *XYat,*YYat;

/* relevant fields */
double ***wa;
double ***sub;
double ***waout;
double ***rho;

/* free energy terms (may be omitted in principal) */
double ***df_drho;
double ***grad;

//double ***p0;
double ***divSigma;
double ***sigmaXZ, ***sigmaYZ, ***sigmaZZ;
double ***fZ;

VecI grid;
double dx = 0.1;
double dy = 0.1;
double dz = 0.1;						// size of the grid cell

int main (int argc, char **argv){
	PassConsoleParams (argc, argv);
	
	AllocArrays ();
	
	ReadSubFunc ();

	//	PrintSubFile ();
	
	CalcForceZ ();
	
	CalcSigmas ();
	
	CalcDivSigma ();
	
	PrintDiff ();
}

/*******************************************************************************************/

void PassConsoleParams (int argc, char **argv) {
	int argz = 1;
	int grid_flag = 0;

	/* check for path */
	if (strcmp(argv[argz], "-path") == 0) {
		if (strcmp(argv[++argz], "pckr160") == 0) {
			strcpy(path,"/data/isilon/tretyakov/mf_lb/");
		} else if (strcmp(argv[argz], "mbpr") == 0) {
			strcpy(path,"/Users/niktre/simulation/60_mf_lb_sh/");
		} else {
			printf ("Error! Unknown system specification!\n");
		}
		++argz;
	} else {
		printf ("Error! No path is given!\n");
	}
	
	/* check for folder */
	if (strcmp(argv[argz], "-folder") == 0) {
		strcpy(folder,argv[++argz]);
		++argz;
	} else {
		printf ("Error! No folder is given!\n");
	}
	
	strcpy (foldername, path);
	strcat (foldername, folder);
	
	/* check for parameters */
	while ((argz < argc) && (argv[argz][0] == '-')) {
		if (strcmp(argv[argz], "-grid") == 0) {
			grid.x = atoi(argv[++argz]);
			grid.y = atoi(argv[++argz]);
			grid.z = atoi(argv[++argz]);
			grid_flag = 1;
		}
		argz++;
	}
}

/*******************************************************************************************/

void ReadSubFunc () {
	double _i, _j, _k;
	double _rho, _df_drho, _grad;
	double _sub, _wa;
	double _wa_temp, _sub_temp;
	char str [120];
	
	FILE *snap, *Wfields, *Wsub;
	
	sprintf(name,"snapshot_final.dat");
	MAKE_FILENAME(fullname_snap,name);
	snap=fopen(fullname_snap,"r");

	// skip first 10 words (header of tecplot)
	for (int i=0; i<10; i++) {
		fscanf (snap, "%s", str);
	}

	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				/* read in the field and the substrate potential */
				fscanf(snap,"%lf %lf %lf %lf %lf %lf %lf %lf \n",
							 &_i, &_j, &_k, &rho[i][j][k], &df_drho[i][j][k], &grad[i][j][k], &_sub_temp, &_wa_temp);
			}
		}
	}

	fclose(snap);
	printf("finished snapshot reading\n");

	// read filename of the fields for restart
	sprintf(name,"Wfields_final.dat");
	MAKE_FILENAME(fullname_Wfields,name);
	Wfields=fopen(fullname_Wfields,"r");

	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				fscanf(Wfields,"%lf \n",&_wa);
				wa[i][j][k] = _wa;
			}
		}
	}
	fclose(Wfields);
	printf("finished Wfields reading\n");
	
	// read substrate potential for restart
	MAKE_FILENAME(fullname_Wsub,"substrate.dat");
	Wsub = fopen(fullname_Wsub,"r");
	
	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				fscanf(Wsub,"%lf \n",&_sub);
				sub[i][j][k] = _sub;
			}
		}
	}
	fclose(Wsub);
	printf("finished substrate reading\n");
	
	// fill in array borders through periodic boundaries //
	// x-direction
	for(int k = 1; k < grid.z+1; k++){
		for(int j = 1; j < grid.y+1; j++){
			rho[0][j][k] = rho[grid.x][j][k]; rho[grid.x+1][j][k] = rho[1][j][k];
		}
	}
	
	// y-direction
	for(int i = 1; i < grid.x+1; i++){
		for(int k = 1; k < grid.z+1; k++){
			rho[i][0][k] = rho[i][grid.y][k]; rho[i][grid.y+1][k] = rho[i][1][k];
		}
	}
	
	// z-direction
	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			rho[i][j][0] = 0.;
			rho[i][j][grid.z+1] = rho[i][j][grid.z-1];
		}
	}
}

/*******************************************************************************************/

/* write down the sub potential into the output file */
void PrintSubFile () {
	FILE *snap, *Wsub, *Wfields;

	MAKE_FILENAME(fullname_snap,"snapshot2.dat");
	
	snap=fopen(fullname_snap,"a");
	for (int i = 1; i < grid.x+1; i++){
		for (int j = 1; j < grid.y+1; j++){
			for (int k = 1; k < grid.z+1; k++){
				fprintf(snap,"%g %g %g %g %g %g %g %g\n",
								(i - 0.5)*dx, (j - 0.5)*dy, (k - 0.5)*dz,
								rho[i][j][k], df_drho[i][j][k], grad[i][j][k],
								sub[i][j][k], wa[i][j][k]);
			}
		}
	}
	fclose(snap);
	
	MAKE_FILENAME(fullname_Wsub,"substrate2.dat");
	
	Wsub=fopen(fullname_Wsub,"a");
	for (int i = 1; i < grid.x+1; i++){
		for (int j = 1; j < grid.y+1; j++){
			for (int k = 1; k < grid.z+1; k++){
				fprintf(Wsub,"%16.12f \n", sub[i][j][k]);
			}
		}
	}
	fclose(Wsub);
	
	MAKE_FILENAME(fullname_Wfields,"Wfields2.dat");
	
	Wfields=fopen(fullname_Wfields,"a");
	for (int i = 1; i < grid.x+1; i++){
		for (int j = 1; j < grid.y+1; j++){
			for (int k = 1; k < grid.z+1; k++){
				fprintf(Wfields,"%16.12f \n", wa[i][j][k]);
			}
		}
	}
	fclose(Wfields);
}

/*******************************************************************************************/

void CalcForceZ () {
	double invDZ = .5 / dz;
	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z; k++){
				fZ[i][j][k] = - rho[i][j][k] * (sub[i][j][k+1] - sub[i][j][k-1]) * invDZ;
			}
		}
	}
	printf("calculated fZ\n");
}

/*******************************************************************************************/

void CalcSigmas () {
	double invDX = .5 / dx;
	double invDY = .5 / dy;
	double invDZ = .5 / dz;
	double p0;
	double grad2rho;
	double gradrho2;
	double WabSec;
	for(int i = 1; i < grid.x; i++){
		for(int j = 1; j < grid.y; j++){
			for(int k = 1; k < grid.z; k++){
				grad2rho = (rho[i+1][j][k] - 2. * rho[i][j][k] + rho[i-1][j][k]) * 4. * SQR(invDX) +
									 (rho[i][j+1][k] - 2. * rho[i][j][k] + rho[i][j-1][k]) * 4. * SQR(invDY) +
									 (rho[i][j][k+1] - 2. * rho[i][j][k] + rho[i][j][k-1]) * 4. * SQR(invDZ);
				gradrho2 = SQR((rho[i+1][j][k] - rho[i-1][j][k]) * invDX) +
									 SQR((rho[i][j+1][k] - rho[i][j-1][k]) * invDY) +
									 SQR((rho[i][j][k+1] - rho[i][j][k-1]) * invDZ);
				WabSec = KAPPA * (rho[i][j][k] * grad2rho + .5 * gradrho2);
				sigmaXZ[i][j][k] = KAPPA * (rho[i+1][j][k] - rho[i-1][j][k]) * invDX * (rho[i][j][k+1] - rho[i][j][k-1]) * invDZ - WabSec;

				sigmaYZ[i][j][k] = KAPPA * (rho[i][j+1][k] - rho[i][j-1][k]) * invDY * (rho[i][j][k+1] - rho[i][j][k-1]) * invDZ - WabSec;
				
				sigmaZZ[i][j][k] = KAPPA * (rho[i][j][k+1] - rho[i][j][k-1]) * invDZ * (rho[i][j][k+1] - rho[i][j][k-1]) * invDZ - WabSec;

				p0 = - 0.5 * V11 * SQR(rho[i][j][k]) - W111 * CUBE(rho[i][j][k]) / 3. +
						 rho[i][j][k] * (1 + wa[i][j][k] - sub[i][j][k]);
				sigmaZZ[i][j][k] += p0;
			}
		}
	}
	printf("calculated sigmaXZ, sigmaYZ and sigmaZZ\n");
}

/*******************************************************************************************/

void CalcDivSigma () {
	double invDX = .5 / dx;
	double invDY = .5 / dy;
	double invDZ = .5 / dz;
	for(int i = 2; i < grid.x-1; i++){
		for(int j = 2; j < grid.y-1; j++){
			for(int k = 2; k < grid.z-1; k++){
				divSigma[i][j][k] = (sigmaXZ[i+1][j][k] - sigmaXZ[i-1][j][k]) * invDX +
														(sigmaYZ[i][j+1][k] - sigmaYZ[i][j-1][k]) * invDY +
														(sigmaZZ[i][j][k+1] - sigmaZZ[i][j][k-1]) * invDZ;
			}
		}
	}
	printf("calculated divergence of sigma\n");
}

/*******************************************************************************************/

void PrintDiff () {
	double 	kbT = 41.16404397;

	FILE *diff;
	
	MAKE_FILENAME(fullname_diff,"microFB.dat");
	
	diff=fopen(fullname_diff,"a");
	fprintf(diff,"VARIABLES = \"X\", \"Y\", \"Z\",\"sigmaXZ\",\"sigmaYZ\",\"sigmaZZ\",\"divSigma\",\"denForce\",\"diff\",\n");
	fprintf(diff,"ZONE I=%d, J=%d, K=%d, F=POINT\n", grid.x - 3, grid.y - 3, grid.z - 3);
	for (int i = 2; i < grid.x-1; i++){
		for (int j = 2; j < grid.y-1; j++){
			for (int k = 2; k < grid.z-1; k++){
				fprintf(diff,"%5.2f %5.2f %5.2f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
								i*dx, j*dy, k*dz,
								kbT * sigmaXZ[i][j][k],kbT * sigmaYZ[i][j][k],kbT * sigmaZZ[i][j][k],
								kbT * divSigma[i][j][k], kbT * fZ[i][j][k], kbT * (divSigma[i][j][k] - fZ[i][j][k]) );
			}
		}
	}
	fclose(diff);
	printf("finished output. Thank you very much!\n");
}

/*******************************************************************************************/

void AllocArrays () {
	ALLOC_MEM (koeff_x_semi, grid.x+1, double);
	ALLOC_MEM (koeff_y_semi, grid.y+1, double);
	ALLOC_MEM (koeff_z_open, grid.z+1, double);
	
	ALLOC_MEM (XYa, grid.x+1, double);
	ALLOC_MEM (YYa, grid.y+1, double);
	ALLOC_MEM (XYat, grid.x+1, double);
	ALLOC_MEM (YYat, grid.y+1, double);
	
	ALLOC_MEM3 (wa, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sub, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (rho, grid.x+2, grid.y+2, grid.z+2, double);
	
	ALLOC_MEM3 (df_drho, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (grad, grid.x+1, grid.y+1, grid.z+1, double);
	
	ALLOC_MEM3 (divSigma, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (fZ, grid.x+1, grid.y+1, grid.z+1, double);
	
	ALLOC_MEM3 (sigmaXZ, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sigmaYZ, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sigmaZZ, grid.x+1, grid.y+1, grid.z+1, double);
}
