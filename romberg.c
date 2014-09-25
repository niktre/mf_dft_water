//	Romberg integration from the "Numerical Recipies for C"
//	with some adaptations and account for our special file input
//  romberg.c
//  mf_lb
//  ls
//  Created by niktre on 09.03.14.
//  Copyright (c) 2014 niktre. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x))
#define EPS 1.0e-6		// fractional accuracy desired, as determined by the extrapolation error estimate;
#define JMAX 20			// limits the total number of steps
#define JMAXP (JMAX+1) 
#define K 5				// the number of points used in the extrapolation
#define PIO2 1.5707963
#define NR_END 1
#define FREE_ARG char*

float func (float);
float qromb (float (float), float, float);
float trapzd(float (float), float, float, int);
void polint(float *, float *, int, float, float *, float *);
float fint(float);
void free_vector(float (*), long, long);
float *vector(long, long);

/* test function */
float func(float x) {
	return x*x*(x*x-2.)*sin(x);
}

int main(void) {
	float a=0.0, b = PIO2, s, t;
	float left_val, right_val;
	
	left_val = fint(a);
	right_val = fint(b);
	
	printf("Integral of func computed with QROMB\n\n");
	printf("Actual value of the integral is %12.6f\n", right_val-left_val);
	s = qromb(func,a,b);
	printf("Result from routine QROMB is %11.6f\n",s);
	return 0;
}

float qromb(float (*func)(float), float a, float b) {
//Returns the integral of the function func from a to b. Integration is performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.

//	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
//	float trapzd(float (*func)(float), float a, float b, int n);
	
	float ss,dss;
	float s[JMAXP],h[JMAXP+1];
	// These store the successive trapezoidal approximations and their relative stepsizes.
	int j;
	
	h[1]=1.0;
	for (j=1; j<=JMAX; j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];	// This is a key step: The factor is 0.25 even though the stepsize is decreased by only 0.5. This makes the extrapolation a polynomial in h^2 as allowed by equation (4.2.1), not just a polynomial in h.
	}
	printf ("Too many steps in routine qromb");
	return 0.0;			// Never get here.
}

float trapzd(float (*func)(float), float a, float b, int n) {
/*This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
￼as a pointer to the function to be integrated between limits a and b, also input. When called with
n=1, the routine returns the crudest estimate of int^b_a f (x)dx. Subsequent calls with n=2,3,... a
(in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
*/
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1; j<n-1; j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;				// This is the spacing of the points to be added.
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);	// This replaces s by its refined value.
		return s;
	}
}

void polint(float xa[], float ya[], int n, float x, float *y, float *dy) {
// Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) = yai,i = 1,...,n, then the returned value y = P(x).

	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;
	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1; i<=n; i++) { 	// Here we find the index ns of the closest table entry,
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1; m<n; m++) {
		for (i=1; i<=n-m; i++) {
			ho=xa[i]-x;		// and initialize the tableau of c’s and d’s.
							// This is the initial approximation to y.
							// For each column of the tableau,
							// we loop over the current c’s and d’s and update them.
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
				printf("Error in routine polint");
			//This error can occur only if two input xa’s are (to within roundoff) identical.
			den=w/den;
			d[i]=hp*den;	// Here the c’s and d’s are updated.
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		// After each column in the tableau is completed, we decide which correction, c or d, we want to add to our accumulating value of y, i.e., which path to take through the tableau—forking up or down. We do this in such a way as to take the most “straight line” route through the tableau to its apex, updating ns accordingly to keep track of where we are. This route keeps the partial approximations centered (insofar as possible) on the target x. The last dy added is thus the error indication.
			}
	free_vector(d,1,n); free_vector(c,1,n);
}

/* Integral of the test function func, i.e., test result */
float fint(float x) {
	return 4.0*x*(x*x-7.)*sin(x)-(pow(x,4.)-14.*x*x+28.)*cos(x);
}

/* allocate a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh) {
	float *v;
	v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) printf ("allocation error in vector()");
	return v-nl+NR_END;
}

/* free a vector allocated with vector() */
void free_vector(float *v, long nl, long nh) {
	free((FREE_ARG) (v+nl-NR_END));
}