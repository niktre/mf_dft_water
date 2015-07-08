#ifndef V11
	#define V11 -1.086641506142464
#endif

#ifndef W111
	#define W111 0.023102120829070
#endif

#ifndef TWRITE
	#define TWRITE 50000
#endif

void InitParameters () {
	extern VecR L, corr;
	extern VecI grid;
	extern double dx, dy, dz;
	extern double rCut, sigma_sub, eps;
	extern double nMol;
	extern int nDims;
	extern int NVT;
	extern int filled_init;
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
	printf("filled_init %d \n", filled_init);
	
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
	extern int filled_init;
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
	fprintf(inParam, "filled_init %d \n", filled_init);
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
	int filled_init_flag = 0;
	
	extern double rho_liq;
	extern VecR L;
	extern VecI grid;
	extern VecR corr;
	extern double eps;
	extern double nMol;
	extern int nDims;
	extern int NVT;
	extern int filled_init;
	
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
		} else if (strcmp(argv[argz], "-filled_init") == 0) {
			filled_init = atof(argv[++argz]);
			filled_init_flag = 1;
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
	
	if (filled_init_flag == 0) {
		filled_init = 0;
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

/* write down the sub potential into the output file */
void PrintSubFile () {
	extern int myRank, MPIsize;
	extern char foldername[], fullname_Wsub[], fullname_profiling[];
	extern VecI grid;
	extern double ***sub, seconds;
	extern time_t begin;
	time_t endInit;

	FILE *Wsub;
  
	MPI_Status status;
	int *sbuf=NULL, *rbuf=NULL;

	sbuf = (int*) malloc(sizeof(int));
	rbuf = (int*) malloc(sizeof(int));

	sbuf[0] = 0;
  
	MAKE_FILENAME(fullname_Wsub,"substrate2.dat");
  
	if (myRank == 0) {
		Wsub=fopen(fullname_Wsub,"a");
		for (int i = 1; i < grid.x+1; i++){
			for (int j = 1; j < grid.y+1; j++){
				for (int k = 1; k < grid.z+1; k++){
					fprintf(Wsub,"%16.12f \n", sub[i][j][k]);
				}
			}
		}
		fclose(Wsub);
		
		if (MPIsize > 1) {
			MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 101, MPI_COMM_WORLD);
		} else {
			// do not send anything, there is only 1 CPU
		}
		
	} else {
		MPI_Recv (rbuf, 1, MPI_INT, myRank - 1, 101, MPI_COMM_WORLD, &status);
		Wsub=fopen(fullname_Wsub,"a");
		for (int i = 1; i < grid.x+1; i++){
			for (int j = 1; j < grid.y+1; j++){
				for (int k = 1; k < grid.z+1; k++){
					fprintf(Wsub,"%16.12f \n", sub[i][j][k]);
				}
			}
		}
		fclose(Wsub);
		if (myRank != MPIsize - 1) {
			MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 101, MPI_COMM_WORLD);
		} else {
		}
	}

	if (myRank == MPIsize - 1) {
		/* output into profiling statistics */
		time(&endInit);
		seconds = difftime(endInit,begin);
		
		FILE *timeprof;
		MAKE_FILENAME(fullname_profiling, "timeprof.dat");
		timeprof = fopen(fullname_profiling,"a");
		fprintf (timeprof, "initialisation and substrate are done in %g seconds\n", seconds);
		fclose(timeprof);
		
		printf ("initialisation and substrate are done in %g seconds\n", seconds);
	}
	
}
void PrintSnapField (int _id, int _kk) {
	extern int write_snapshot, myRank, MPIsize, iteration_num, globStartIdx, subGridX;
	extern char foldername[], fullname_iter[], fullname_profiling[], name[];
	extern char fullname_snap[], fullname_Wfields[];
	extern double total_N, seconds, dz, dx, dy;
	extern VecI grid;
	extern time_t begin, end;
	extern double ***rho, ***sub, ***wa, ***df_drho, ***grad;
	FILE *snapshot, *Wrestart;
	
	MPI_Status status;
	int *sbuf=NULL, *rbuf=NULL;
	
	sbuf = (int*) malloc(sizeof(int));
	rbuf = (int*) malloc(sizeof(int));
	
	sbuf[0] = 0;
	
	if (_id == 0 && _kk % TWRITE == 0) {
		if (write_snapshot == 1 && myRank == 0) {
			// save current interation number
			FILE *iterkeeper;
			MAKE_FILENAME(fullname_iter,"iteration.dat");
			iterkeeper = fopen(fullname_iter,"w");
			fprintf(iterkeeper,"%d %10.8f\n",_kk + iteration_num, total_N);
			fclose(iterkeeper);
			
			/* output into profiling statistics */
			time(&end);
			seconds = difftime(end,begin);
			begin = end;
			
			FILE *timeprof;
			MAKE_FILENAME(fullname_profiling, "timeprof.dat");
			timeprof = fopen(fullname_profiling,"a");
			fprintf (timeprof, "%d %g\n", _kk + iteration_num, seconds);
			fclose(timeprof);
			
		} else {
		}
	} else if (_id == 1 && write_snapshot == 1) {
		FILE *snapshot, *Wrestart;
		
		// set the name of the files for all CPUs and prepare the HEADERS by CPU0
		// create snapshot header and open the file
		sprintf(name,"snapshot_%d.dat",_kk + iteration_num);
		MAKE_FILENAME(fullname_snap,name);
		if (myRank == 0) {
			snapshot = fopen(fullname_snap,"w");
			fprintf (snapshot,"VARIABLES = \"X\", \"Y\", \"Z\",\"density\",\"df_drho\",\"grad\",\"subpot\",\"wa\",\n");
			fprintf (snapshot,"ZONE I=%d, J=%d, K=%d, F=POINT\n",grid.z, grid.y, subGridX);
			fclose(snapshot);
		}
		// make filename for the fields and open the file
		sprintf(name,"Wfields_%d.dat",_kk + iteration_num);
		MAKE_FILENAME(fullname_Wfields,name);
		
		if (myRank == 0) {
			// saves the snapshot (the rest) and fields for future restart
			snapshot = fopen(fullname_snap,"a");
			Wrestart = fopen(fullname_Wfields,"a");
			for (int i = 1; i < grid.x+1; i++){
				for (int j = 1; j < grid.y+1; j++){
					for (int k = 1; k < grid.z+1; k++){
						fprintf(snapshot,"%g %g %g %g %g %g %g %g\n",
					  (globStartIdx + i - 0.5)*dx, (j - 0.5)*dy, (k - 0.5)*dz,
					  rho[i][j][k], df_drho[i][j][k], grad[i][j][k],
					  sub[i][j][k], wa[i][j][k]);
						fprintf(Wrestart,"%16.12f \n", wa[i][j][k]);
					}
				}
			}
			fclose(snapshot); fclose(Wrestart);
			
			if (MPIsize > 1) {
				MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 103, MPI_COMM_WORLD);
			} else {
				// do not send anything, there is only 1 CPU
			}
			
		} else {
			MPI_Recv (rbuf, 1, MPI_INT, myRank - 1, 103, MPI_COMM_WORLD, &status);
			// saves the snapshot (the rest) and fields for future restart
			snapshot = fopen(fullname_snap,"a");
			Wrestart = fopen(fullname_Wfields,"a");
			for (int i = 1; i < grid.x+1; i++){
				for (int j = 1; j < grid.y+1; j++){
					for (int k = 1; k < grid.z+1; k++){
						fprintf(snapshot,"%g %g %g %g %g %g %g %g\n",
								  (globStartIdx + i - 0.5)*dx, (j - 0.5)*dy, (k - 0.5)*dz,
								  rho[i][j][k], df_drho[i][j][k], grad[i][j][k],
								  sub[i][j][k], wa[i][j][k]);
						fprintf(Wrestart,"%16.12f \n", wa[i][j][k]);
					}
				}
			}
			fclose(snapshot); fclose(Wrestart);
			if (myRank != MPIsize - 1) {
				MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 103, MPI_COMM_WORLD);
			} else {
			}
		}
	} else {
	}
}

void ReadSubFile () {
	extern int myRank, MPIsize;
	
	MPI_Status status;
	int *sbuf=NULL, *rbuf=NULL;
	
	sbuf = (int*) malloc(sizeof(int));
	rbuf = (int*) malloc(sizeof(int));
	
	sbuf[0] = 0;
 
	if (myRank == 0) {
		ReadSubFunc();
		
		if (MPIsize > 1) {
			MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 102, MPI_COMM_WORLD);
		} else {
			// do not send anything, there is only 1 CPU
		}
	} else {
		MPI_Recv (rbuf, 1, MPI_INT, myRank - 1, 102, MPI_COMM_WORLD, &status);

		ReadSubFunc();

		if (myRank != MPIsize - 1) {
			MPI_Send (sbuf, 1, MPI_INT, myRank + 1, 102, MPI_COMM_WORLD);
		} else {
		}
	}
}

void ReadSubFunc () {
	extern char name[], foldername[], fullname_iter[], fullname_Wfields[], fullname_Wsub[];
	extern int iteration_num, globStartIdx;
	extern double total_N, dz, rho_liq, ***df_drho, ***grad, ***wa, ***sub;
	
	extern VecI grid;
	extern VecR corr;

	FILE *iterkeeper, *Wfields, *Wsub;

	// read iteration number and total number of particles for restart
	MAKE_FILENAME(fullname_iter,"iteration.dat");
	iterkeeper = fopen(fullname_iter,"r");
	fscanf(iterkeeper,"%d %lf\n",&iteration_num, &total_N);
	fclose(iterkeeper);
	
	// read filename of the fields for restart
	if (iteration_num != 0) {
		sprintf(name,"Wfields_%d.dat",iteration_num);
		MAKE_FILENAME(fullname_Wfields,name);
		Wfields=fopen(fullname_Wfields,"r");
		printf ("fullname_Wfields is %s\n", fullname_Wfields);
	}
	
	// read substrate potential for restart
	MAKE_FILENAME(fullname_Wsub,"substrate.dat");
	Wsub = fopen(fullname_Wsub,"r");
	
	// skip some lines because they should be read by another processor
	int skipLines = 0;
	double tempBuf = 0.;
	while (skipLines != globStartIdx*grid.y*grid.z) {
		if (iteration_num != 0) {
			fscanf(Wfields,"%lf \n",&tempBuf);
		}
		fscanf(Wsub,"%lf \n",&tempBuf);
		skipLines++;
	}
	
	for(int i = 1; i < grid.x+1; i++){
		for(int j = 1; j < grid.y+1; j++){
			for(int k = 1; k < grid.z+1; k++){
				/* set free energy terms to zero*/
				SET_TO_ZERO (df_drho,i,j,k);
				SET_TO_ZERO (grad,i,j,k);
				
				/* read in the field and the substrate potential */
				if (iteration_num == 0) {
					// create fiels from scratch
					wa[i][j][k] = 0.;
					if ((0.5 * dz + (k - 1) * dz) > 2. * corr.z){
						wa[i][j][k] = V11*rho_liq + W111*SQR(rho_liq);
						total_N = total_N + 1.;
					}
				} else {
					fscanf(Wfields,"%lf \n",&wa[i][j][k]);
				}
				fscanf(Wsub,"%lf \n",&sub[i][j][k]);
			}
		}
	}
	fclose(Wsub);
	printf ("Finished with scanning substrate file!\n");
	printf ("total_N is %10.8f\n", total_N);
	
	if (iteration_num != 0) {
		fclose(Wfields);
	}
}
	
	
void BcastConsoleParams () {
  extern int myRank;
  
  extern char path[];
  extern char folder[];
  
  extern int nDims;
  extern int NVT;
  extern int filled_init;
  extern int restart;
  
  extern double rho_liq;
  extern double eps;
  extern double nMol;
  extern double lambda;
  extern double dx, dy, dz, dx2, dy2, dz2, volume;
  
  extern VecI grid;
  extern VecR L;
  extern VecR corr;
  
  int countInt = 4;
  int countDouble = 11;
  int countVecI = 3;
  int countVecR = 6;

//  char *pbuf=NULL, *fbuf=NULL;
  int *ibuf=NULL, *vecIbuf=NULL;
  double *dbuf=NULL, *vecRbuf=NULL;

    	if (myRank == 0) {
//  pbuf = malloc((strlen(path)+1)*sizeof(char));
//  fbuf = malloc((strlen(folder)+1)*sizeof(char));
  ibuf = malloc(countInt*sizeof(int));
  dbuf = malloc(countDouble*sizeof(double));
  vecIbuf = malloc(countVecI*sizeof(int));
  vecRbuf = malloc(countVecR*sizeof(double));
  
/*	  for (int i = 0; i < strlen(path)+1; i++) {
		 pbuf[i]=path[i];
	  }
	  for (int i = 0; i < strlen(folder)+1; i++) {
		 fbuf[i]=folder[i];
	  }
*/
		  ibuf[0] = nDims;
		  ibuf[1] = NVT;
		  ibuf[2] = filled_init;
		  ibuf[3] = restart;
//	  &fbuf = folder;
//	  printf ("length of path and folder are %d and %d\n", strlen(path), strlen(folder));
//  /* Find out process rank */
//  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	}
//  MPI_Bcast(pbuf, strlen(path)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
//  MPI_Bcast(fbuf, strlen(folder)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(ibuf, countInt, MPI_INT, 0, MPI_COMM_WORLD);
  
  
  if (myRank != 0) {
	 nDims = ibuf[0]; NVT = ibuf[1]; filled_init = ibuf[2]; restart = ibuf[3];
//  path = &pbuf;
//	 folder = &fbuf;
  }
	 
//  printf("Hello world, I am process %d and the path is %s folder is %s \n", myRank, path, folder);
}