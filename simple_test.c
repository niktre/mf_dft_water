//
//  simple_test.c
//  mf_lb
//
//  Created by niktre on 20.03.14.
//  Copyright (c) 2014 niktre. All rights reserved.
//

#include <stdio.h>

int main () {
	int i, j, k;
	int grid_x = 10;
	int grid_y = 10;
	int grid_z = 10;
	
	double XYa[grid_x+1];
	double YYa[grid_y+1];
	double koeff_x_close[grid_x+1];
	double koeff_y_close[grid_y+1];
	double koeff_z_open[grid_z+1];
	
	double total_N;
	
	koeff_x_close[1] 	  = koeff_y_close[1]    = 25./24.;
	koeff_x_close[2] 	  = koeff_y_close[2]    = 25./24.;
	koeff_x_close[3] 	  = koeff_y_close[3]    = 3./4.;
	
/*	koeff_x_close[grid_x-2]	= koeff_y_close[grid_y-2]	= 23./24.;
	koeff_x_close[grid_x-1]	= koeff_y_close[grid_y-1]	= 7./6.;
	koeff_x_close[grid_x]	= koeff_y_close[grid_y]		= 3./8.;
*/
	koeff_x_close[grid_x-2]	= koeff_y_close[grid_y-2]	= 3./4.;
	koeff_x_close[grid_x-1]	= koeff_y_close[grid_y-1]	= 25./24.;
	koeff_x_close[grid_x]	= koeff_y_close[grid_y]		= 25./24.;
	
	koeff_z_open[1] = koeff_z_open[grid_z]		= 11./ 24.;
	koeff_z_open[2] = koeff_z_open[grid_z-1]	= 7./ 8.;
	koeff_z_open[3] = koeff_z_open[grid_z-2]	= 5./ 3.;

/*	koeff_z_open[1] = koeff_z_open[grid_z]		= 13./ 12.;
	koeff_z_open[2] = koeff_z_open[grid_z-1]	= 7./ 8.;
	koeff_z_open[3] = koeff_z_open[grid_z-2]	= 5./ 3.;
*/
	for(i = 4; i <= grid_x - 3; i++){
		koeff_x_close[i] = 1.;
	}
	
	for(i = 4; i <= grid_y - 3; i++){
		koeff_y_close[i] = 1.;
	}
	
	for(i = 4; i <= grid_x - 3; i++){
		if (i % 2 == 0) koeff_z_open[i] = 4./3.;
		else koeff_z_open[i] = 2./3.;
	}
/*	for(i = 4; i <= grid_z - 3; i++){
		koeff_z_open[i] = 1.;
	}
*/
	total_N = 0.;
	
	for (i = 1; i < grid_x+1; i++){
		XYa[i] = 0.;
//		for (j = 1; j < grid_y+1; j++){
//			YYa[j] = 0.;
//			for(k = 1; k < grid_z+1; k++){
//				YYa[j] = YYa[j] + koeff_y_close[j];
//			}
//			XYa[i] = XYa[i] + koeff_y_close[j]*YYa[j];
//		}
		total_N = total_N +  koeff_x_close[i];
	}
//	total_N = YYa
	printf ("total_N is %8.4f\n", total_N);
}