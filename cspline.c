/*******************************************************************************
 * This work is about cubic spline interpolation algorithms. In order to make 
   life easier, I am going to code these algorithms step by step. However, you 
   must read the references that I mention below .
 * For general understanding please refer to the Wikia page
   http://en.wikipedia.org/wiki/Spline_interpolation
 * To understand the program, please refer to
   http://www.physics.arizona.edu/~restrepo/475A/Notes/sourcea-/node35.html
 * And you should not ignore the book:
	Numerical Recipes in c: The Art of Scientific Computing by
	William H. Press,
	Brian P. Flannery,
	Saul a. Teukolsky, and
	William T. Vetterling .
	Copyright 1988 (and 1992 for the 2nd edition)
	http://www.nr.com/
 * These codes are written for studying purpose only. However, you could feel
   free to modify and improve as your need.
 * Author: Thanh Nguyen
 * Email: thanhng1985@gmail.com 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cspline.h"
/*----------------------------------------------------------------------------*/
#define MALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)malloc((num) * sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

static 		FILE *fp= NULL;

/*----------------------------------------------------------------------------*/
/* 
 * Natural cubic spline 
 * params:
	INPUTS:
	- x[]: together with y[] creates a tabulated data
	- y[] = f(x[]): f(x) can be linear or nonlinear
	- nx = length of x[]
	- ny = length of y[]
	- n: number of points will be taken to calculate, n <= min(nx, ny)
	OUTPUTS:	
	- ypp[]: second derivatives of the interpolating function, or fpp(x[])
 */
void ncspline(double x[], double y[], int nx, int ny, int n, double ypp[])
{

	/* m = min (nx, ny)*/
	int nxy= nx > ny ? ny : nx;
	/* n = min (n, m)*/
	n = n > nxy? nxy : n;
	int i;
	//double h[n-1];
	double *h; 
	MALLOC(h, double, n - 1);
	for (i = 0; i < n - 1; i++) 
	{
		h[i] = (x[i+1] - x[i]);
	}
	/* after the for loop, i = n - 1 */
	//double r[n-1];
	double *r;
	MALLOC(r, double, n - 1);

	for (i = 0; i < n - 1; i++)
	{
		r[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
	}
	
	//double A[n-2], B[n-2], C[n-2], D[n-2];
	double *a, *b, *c, *d;
	MALLOC(a, double, n - 2);
	MALLOC(b, double, n - 2);
	MALLOC(c, double, n - 2);
	MALLOC(d, double, n - 2);
	
	for (i = 0; i < n - 2; i++)
	{ 
		a[i] = h[i];
		b[i] = 2.0 * (h[i] + h[i+1]);
		c[i] = h[i+1];
		d[i] = 6.0 * (r[i+1] - r[i]);
	}
	/* after the for loop, i = n - 2 */
	ypp[i+1] = 0.0; // ypp[n-1] = 0.0, natural cubic spline
	d[i-1] = d[i-1] - c[i-1] * ypp[i+1];
	c[i-1] = 0.0;
	
	i = 0;
	ypp[i] = 0.0; // natural cubic spline
	d[i] = 6.0 * (r[i+1] - r[i])- a[i] * ypp[i];
	a[i] = 0.0;
	
	/*
	 * tridialogal matrix algorimth
	 * we can call tridiag(d, n-2, a, b, c) instead, then define ypp(0) and 
	 * ypp(n-1), because d = ypp(1) ... ypp(n-2)
	 */
	//double dp[n-2], dp[n-2];
	double *cp, *dp;
	MALLOC(cp, double, n - 2);
	MALLOC(dp, double, n - 2);
	
	if (b[0] == 0)
	{
		fprintf(stderr,"-E- %s line %d: error 1 in tridiagonal\n",
            __FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	
	cp[0] = c[0] / b[0];
	dp[0] = d[0] / b[0];
	double m = 0.0;
	for (i = 1; i < n - 3; i++)
	{
		m = (b[i] - a[i] * cp[i-1]);
		if (m == 0)
		{
			fprintf(stderr,"-E- %s line %d: error 2 in tridiagonal\n",
            __FILE__,__LINE__);
			exit(EXIT_FAILURE);
		}
		cp[i] = c[i] / m;
		dp[i] = (d[i] - a[i] * dp[i-1]) / m;
	}
	//m = (b[n-3] - a[n-3] * cp[n-4]);
	if (m == 0)
	{
		fprintf(stderr,"-E- %s line %d: error 3 in tridiagonal\n",
            __FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	//dp[n-3] = (d[n-3] - a[n-3] * dp[n-4]) / m;
	// ypp[n-2] = dp[n-3];
	dp[i] = (d[i] - a[i] * dp[i-1]) / m;
	ypp[i+1] = dp[i];
	
	for (i = n - 2; i-- > 1;)
	{
		ypp[i] = dp[i-1] - cp[i-1] * ypp[i+1];
	}

	free(h);
	free(r);
	free(a);
	free(b);
	free(c);
	free(d);
	free(cp);
	free(dp);
}
/*----------------------------------------------------------------------------*/
/* 
 * Clamped cubic spline 
 * params:
	INPUTS:
	- x[]: together with y[] creates a tabulated data
	- y[] = f(x[]): f(x) can be linear or nonlinear
	- nx = length of x[]
	- ny = length of y[]
	- n: number of points will be taken to calculate, n <= min(nx, ny)
	- yp0 = the first derivate at x[0]
	- ypn_1 = the first derivate at x[n-1]
	OUTPUTS:	
	- ypp[]: second derivatives of the interpolating function, or fpp(x[])
 */
void clcspline(double x[], double y[], int nx, int ny, int n, double yp0, double ypn_1, double ypp[])
{
	/* m = min (nx, ny)*/
	int m = nx > ny ? ny : nx;
	/* n = min (n, m)*/
	n = n > m ? m : n;
	
	double *h; 
	MALLOC(h, double, n - 1);
	int i;
	
	for (i = 0; i < n - 1; i++) 
	{
		h[i] = (x[i+1] - x[i]);
	}
	/* after the for loop, i = n - 1, but h stop at h[n-2] */
	double *r; 
	MALLOC(r, double, n - 1);
	
	for (i = 0; i < n - 1; i++)
	{
		r[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
	}
	
	double *a, *b, *c, *d; 
	MALLOC(a, double, n);
	MALLOC(b, double, n);
	MALLOC(c, double, n);
	MALLOC(d, double, n);
	
	a[0] = 0.0;
	b[0] = 2.0 * h[0]; 
	c[0] = h[0];
	d[0] = 6.0 * (r[0] - yp0);
	for (i = 1; i < n - 1; i++)
	{ 
		a[i] = h[i-1];
		b[i] = 2.0 * (h[i-1] + h[i]);
		c[i] = h[i];
		d[i] = 6.0 * (r[i] - r[i-1]);
	}
	a[n-1] = h[n-2];
	b[n-1] = 2 * h[n-2];
	c[n-1] = 0.0;
	d[n-1] = 6.0 * (ypn_1 - r[n-2]);

	tridiag(d, n, a, b, c);
	
	for (i = 0; i < n; i++)
	{
		ypp[i] = d[i];
	}
	
	free(a);
	free(b);
	free(c);
	free(c);
	free(h);
	free(r);
}
/*----------------------------------------------------------------------------*/
/*
 This code is based on the cubic spline interpolation code presented in:
 Numerical Recipes in C: The Art of Scientific Computing
 by
 William H. Press,
 Brian P. Flannery,
 Saul A. Teukolsky, and
 William T. Vetterling .
 Copyright 1988 (and 1992 for the 2nd edition)
 */
void csplint(double x[], double y[], double ypp[], int n, double xi, double *yi)
{

	int			klo, khi, k;
	double		h, b, a;
	static int	pklo = 0,pkhi = 1;

	if(x[pklo] <= xi && x[pkhi] > xi)
	{
		klo = pklo;
		khi = pkhi;
	}
	else
	{
		klo = 0;
		khi = n - 1;
		while(khi - klo > 1)
		{
			k = (khi + klo) >> 1;
			if(x[k] > xi) 
				khi = k;
			else
				klo = k;
		}
	}

	h = x[khi] - x[klo];
	if(h == 0)
	{
		fprintf(stderr,"-E- %s line %d: Bad x input to function csplint()\n",
            __FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	a = (x[khi] - xi)/h;
	b = (xi - x[klo])/h;
	*yi = a*y[klo] + b*y[khi] +
		((a*a*a - a)*ypp[klo] + (b*b*b - b)*ypp[khi])*(h*h)/6.0;
}
/*----------------------------------------------------------------------------*/
/* 	   http://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
 */

void tridiag(double x[], int N, double a[], double b[], double c[]) 
{
    int i;
 
    /* Allocate scratch space. */
    double* cp;
	MALLOC(cp, double, N);
	
	if (b[0] == 0)
	{
		fprintf(stderr,"-E- %s line %d: error 1 in tridiagonal\n",
            __FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	
    cp[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
 
    /* loop from 1 to N - 1 inclusive */
    for (i = 1; i < N; i++) 
	{
        double m = (b[i] - a[i] * cp[i - 1]);
		if (m == 0)
		{
			fprintf(stderr,"-E- %s line %d: error 2 in tridiagonal\n",
				__FILE__,__LINE__);
			exit(EXIT_FAILURE);
		}
        cp[i] = c[i] / m;
        x[i] = (x[i] - a[i] * x[i - 1]) / m;
    }
 
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (i = N - 1; i-- > 0; )
        x[i] = x[i] - cp[i] * x[i + 1];
 
    /* free scratch space */
    free(cp);
}
/*----------------------------------------------------------------------------*/
void getxytable(double x[], double y[], int *n)
{
	*n = 0;
	fp = fopen("xytable.txt", "r");
	if (fp == NULL) {
        perror("Failed to open file \"xytable.txt\"");
        exit(EXIT_FAILURE);
    }
	else
	{
		int size = 32, pos;
		int c;
		char *buffer = (char *)malloc(size);
		char *tbuffer = (char *)malloc(size/2);
		do
		{ /* read all lines in file */
			pos = 0;
			do
			{ 
				/* read one line */
				c = fgetc(fp);
				if (c != EOF) /* neccessary? */
					buffer[pos++] = c;

				if (pos >= size - 1) 
				{ 
					/* increase buffer length - leave room for 0 */ 
					size *=2;
					buffer = (char *)realloc(buffer, size);
					tbuffer = (char *)realloc(tbuffer, size/2);
				}
				
			} while (c != EOF && c != '\n');
			
			buffer[pos] = '\0';
			tbuffer = strtok(buffer, ", ;");
			if (tbuffer != NULL)
			{
				x[*n] = atof(tbuffer);
				
				tbuffer = strtok((char *)(buffer - tbuffer), ", ;");
				if (tbuffer != NULL)
				{
					y[*n] = atof(tbuffer);
					tbuffer = strtok (NULL, ", ;");
				}
				tbuffer = strtok (NULL, ", ;");
			}
			
			(*n)++;
		} while (c != EOF);
		
		fclose(fp);
		free(buffer);
		free(tbuffer);	  
    }
}

/*----------------------------------------------------------------------------*/
