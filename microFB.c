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

double ***p0;
double ***Wab;
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
	
	PrintSubFile ();
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
			printf("path is %s\n", path);
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
	printf ("fullname_snap is %s\n", fullname_snap);

	// skip first 10 words (header of tecplot)
	for (int i=0; i<10; i++) {
		fscanf (snap, "%s", str);
		printf ("%s\n", str);
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

	// read filename of the fields for restart
	sprintf(name,"Wfields_final.dat");
	MAKE_FILENAME(fullname_Wfields,name);
	Wfields=fopen(fullname_Wfields,"r");
	printf ("fullname_Wfields is %s\n", fullname_Wfields);

	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				fscanf(Wfields,"%lf \n",&_wa);
				wa[i][j][k] = _wa;
			}
		}
	}
	fclose(Wfields);
	
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
	
	// "periodic" boundaries //
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
	fclose(Wsub);
	
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
	
	ALLOC_MEM3 (Wab, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (fZ, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (p0, grid.x+1, grid.y+1, grid.z+1, double);
	
	ALLOC_MEM3 (sigmaXZ, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sigmaYZ, grid.x+1, grid.y+1, grid.z+1, double);
	ALLOC_MEM3 (sigmaZZ, grid.x+1, grid.y+1, grid.z+1, double);
}