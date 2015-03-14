#ifndef V11
	#define V11 -1.086641506142464
#endif

#ifndef W111
	#define W111 0.023102120829070
#endif

void InitParameters () {
	extern VecR L, corr;
	extern VecI grid;
	extern double dx, dy, dz;
	extern double rCut, sigma_sub, eps;
	extern double nMol;
	extern int nDims;
	extern int NVT;
	extern double rho_liq, kbT;
	
	extern char fullname_press[];
	extern char foldername[];

	FILE *wPressure;
	
	printf("NVT-Ensemble %8d\n", NVT);
	printf("Lx, Ly and Lz are %8.4f %8.4f %8.4f\n", L.x, L.y, L.z);
	printf("grid.x, grid.y and grid.z are %8d %8d %8d\n", grid.x, grid.y, grid.z);
	printf("dx, dy and dz are %8.4f %8.4f %8.4f\n",dx, dy, dz);
	printf("rCut %8.4f  sigma_sub %8.4f  eps %8.4f\n", rCut, sigma_sub, eps);
	printf("The pressure is %6.10f \n",  kbT*10.*(rho_liq + (0.5*V11*SQR(rho_liq)) + (2.*W111*CUBE(rho_liq)/3.)));
	printf("corrugation	is %6.2f %6.2f %6.2f\n", corr.x, corr.y, corr.z);
	printf("cavity is %6.2f %6.2f %6.2f\n", L.x - corr.x, L.y - corr.y, L.z - corr.z);
	printf("rCut %8.4f  sigma_sub %8.4f  eps %8.4f\n", rCut, sigma_sub, eps);
	printf("nMol %8.4f \n", nMol);
	printf("nDims %d \n", nDims);
	
	MAKE_FILENAME(fullname_press, "pressure.dat");
	wPressure = fopen(fullname_press,"a");
	fprintf(wPressure,"The pressure is %6.10f \n",  kbT*10.*(rho_liq + (0.5*V11*SQR(rho_liq)) + (2.*W111*CUBE(rho_liq)/3.)));
	fclose(wPressure);
}

void StoreParameters () {
	extern VecR L, corr;
	extern VecI grid;
	extern double rCut, sigma_sub, eps;
	extern double nMol;
	extern int nDims;
	extern int NVT;
	extern double rho_liq, lambda;
	extern int iterations;
	
	extern char fullname_param[];
	extern char foldername[];
	
	FILE *inParam;
	
	MAKE_FILENAME(fullname_param, "in_param.dat");
	inParam = fopen(fullname_param,"a");
	fprintf(inParam, "Info about parameters used\n");
	fprintf(inParam, "NVT-Ensemble %8d\n", NVT);
	fprintf(inParam, "Grid(grid.x,grid.y,grid.z)	%4d %4d %4d\n", grid.x, grid.y, grid.z);
	fprintf(inParam, "Box (L.x,L.y,L.z)	%6.2f %6.2f %6.2f\n", L.x, L.y, L.z);
	fprintf(inParam, "rho_liq	%8.4f\n", rho_liq);
	fprintf(inParam, "iterations	%8d\n", iterations);
	fprintf(inParam, "lambda %8.6f\n", lambda);
	fprintf(inParam, "corrugation	%6.2f %6.2f %6.2f\n", corr.x, corr.y, corr.z);
	fprintf(inParam, "cavity is %6.2f %6.2f %6.2f\n", L.x - corr.x, L.y - corr.y, L.z - corr.z);
	fprintf(inParam, "rCut %8.4f  sigma_sub %8.4f  eps %8.4f\n", rCut, sigma_sub, eps);
	fprintf(inParam, "nMol %8.4f \n", nMol);
	fprintf(inParam, "nDims %d \n", nDims);
	fclose(inParam);
}

/* a function to handle parameters passed through console */
void PassConsoleParams (int argc, char **argv) {
	int argz = 1;
	extern char path[];
	extern char folder[];
	
	extern int restart;
	extern double lambda;
	
	int rho_liq_flag = 0;
	int restart_flag = 0;
	int lambda_flag = 0;
	int L_flag = 0;
	int grid_flag = 0;
	int corr_flag = 0;
	int eps_flag = 0;
	int nMol_flag = 0;
	int nDims_flag = 0;
	int NVT_flag = 0;
	
	extern double rho_liq;
	extern VecR L;
	extern VecI grid;
	extern VecR corr;
	extern double eps;
	extern double nMol;
	extern int nDims;
	extern int NVT;
	
	extern double dx, dy, dz, dx2, dy2, dz2, volume;
	
	/* check for path */
	if (strcmp(argv[argz], "-path") == 0) {
		if (strcmp(argv[++argz], "pckr160") == 0) {
			strcpy(path,"/data/isilon/tretyakov/mf_lb/");
		} else if (strcmp(argv[argz], "vaio") == 0) {
			strcpy(path,"/home/niktre/Documents/MPIP/LB_codes/mf_code/");
		} else if (strcmp(argv[argz], "rocks") == 0) {
			strcpy(path,"/home/tretyakov/mf_freen/");
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
	
	/* check for parameters */
	while ((argz < argc) && (argv[argz][0] == '-')) {
		if (strcmp(argv[argz], "-den") == 0) {
			rho_liq = atof(argv[++argz]);
			rho_liq_flag = 1;
		} else if (strcmp(argv[argz], "-restart") == 0) {
			restart = atof(argv[++argz]);
			restart_flag = 1;
		} else if (strcmp(argv[argz], "-lambda") == 0) {
			lambda = atof(argv[++argz]);
			lambda_flag = 1;
		} else if (strcmp(argv[argz], "-Lbox") == 0) {
			L.x = atof(argv[++argz]);
			L.y = atof(argv[++argz]);
			L.z = atof(argv[++argz]);
			L_flag = 1;
		} else if (strcmp(argv[argz], "-grid") == 0) {
			grid.x = atoi(argv[++argz]);
			grid.y = atoi(argv[++argz]);
			grid.z = atoi(argv[++argz]);
			grid_flag = 1;
		} else if (strcmp(argv[argz], "-corr") == 0) {
			corr.x = atof(argv[++argz]);
			corr.y = atof(argv[++argz]);
			corr.z = atof(argv[++argz]);
			corr_flag = 1;
		} else if (strcmp(argv[argz], "-eps") == 0) {
			eps = atof(argv[++argz]);
			eps_flag = 1;
		} else if (strcmp(argv[argz], "-nMol") == 0) {
			nMol = atof(argv[++argz]);
			nMol_flag = 1;
		} else if (strcmp(argv[argz], "-nDims") == 0) {
			nDims = atof(argv[++argz]);
			nDims_flag = 1;
		} else if (strcmp(argv[argz], "-NVT") == 0) {
			NVT = atof(argv[++argz]);
			NVT_flag = 1;
		}
		argz++;
	}
	
	/* check for set flags */
	if (rho_liq_flag == 0) {
		rho_liq = 33.33;
	}
	
	if (restart_flag == 0) {
		restart = 0;
	}
	
	if (lambda_flag == 0) {
		lambda = 0.001;
		printf("lambda_flag was 0: lambda is %8.4f\n", lambda);
	}
	
	if (L_flag == 0) {
		L.x = 6.4; L.y = 6.4; L.z = 6.4;
		printf("L_flag was 0: Lx, Ly and Lz are %8.4f %8.4f %8.4f\n", L.x, L.y, L.z);
	}
	
	if (grid_flag == 0) {
		grid.x = (int)(L.x * 10.);
		grid.y = (int)(L.y * 10.);
		grid.z = (int)(L.z * 10.);
		printf("grid_flag was 0: grid.x, grid.y and grid.z are %8d %8d %8d\n", grid.x, grid.y, grid.z);
	}
	
	if (corr_flag == 0) {
		corr.x = 0.8; corr.y = 0.8; corr.z = 0.8;
	}
	
	if (eps_flag == 0) {
		eps = 1.;
	}
	
	if (nMol_flag == 0) {
		nMol = 1092.;
	}
	
	if (nDims_flag == 0) {
		nDims = 3;
	}
	
	if (NVT_flag == 0) {
		NVT = 1;
	}
	
	dx = L.x/((double)grid.x);
	dy = L.y/((double)grid.y);
	dz = L.z/((double)grid.z);
	
	dx2 = SQR(dx);
	dy2 = SQR(dy);
	dz2 = SQR(dz);
	
	volume = L.x * L.y * L.z;
}
