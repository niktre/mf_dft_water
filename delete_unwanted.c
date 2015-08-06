#include <stdio.h>
#include <string.h>

char name[120];
char foldername[120];
char path[120];
char folder[120];
char fullname_iter[120];

#define MAKE_FILENAME(full, file)						\
	{strcpy (full, foldername);							\
	 strcat (full,file);}

#define TWRITE 50000

int main (int argc, char **argv){
	double total_N;
	int iteration_num;
	FILE *iterkeeper, *source, *target;
	char ch;
	int argz = 1;

	/* check for path */
	if (strcmp(argv[argz], "-path") == 0) {
		if (strcmp(argv[++argz], "isilon") == 0) {
			strcpy(path,"/data/isilon/tretyakov/mf_lb/");
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

	strcpy (foldername, path);
	strcat (foldername, folder);

	// read iteration number and total number of particles for restart
	MAKE_FILENAME(fullname_iter,"iteration.dat");
	iterkeeper = fopen(fullname_iter,"r");
	fscanf(iterkeeper,"%d %lf\n",&iteration_num, &total_N);
	fclose(iterkeeper);

	char source_file[120], target_file[120];
	/* copy the newest snapshot and Wfields file */
	sprintf(name,"Wfields_%d.dat",iteration_num);
	MAKE_FILENAME(source_file,name);
	sprintf(name,"Wfields_final.dat");
	MAKE_FILENAME(target_file,name);

	source = fopen(source_file, "r");
	target = fopen(target_file, "w");
	while( ( ch = fgetc(source) ) != EOF )
      fputc(ch, target);

	fclose(source);
	fclose(target);

	sprintf(name,"snapshot_%d.dat",iteration_num);
	MAKE_FILENAME(source_file,name);
	sprintf(name,"snapshot_final.dat");
	MAKE_FILENAME(target_file,name);

	source = fopen(source_file, "r");
	target = fopen(target_file, "w");
	while( ( ch = fgetc(source) ) != EOF )
      fputc(ch, target);

	fclose(source);
	fclose(target);

	/* remove old snapshots and Wfields */
	int remove_status, fileNumToDelete;
	char fileToDel[120];
	fileNumToDelete = iteration_num - 2 * TWRITE;
	while (fileNumToDelete > 0) {
		// make filename to delete the file
		sprintf(name,"Wfields_%d.dat",fileNumToDelete);
		MAKE_FILENAME(fileToDel,name);
		remove_status = remove(fileToDel);
					
		sprintf(name,"snapshot_%d.dat",fileNumToDelete);
		MAKE_FILENAME(fileToDel,name);
		remove_status = remove(fileToDel);
		
		// reduce file_number
		fileNumToDelete -= TWRITE;
	}
	return 0;
}
