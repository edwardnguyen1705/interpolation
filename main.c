#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include "cspline.h"
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
	/* clear the command prompt screen */
	system("cls");
	int n;
	/* 
	 * look at the data file to find out number of lines 
	 * delete lines without data
	 */
	double x1[128], y1[128]; // 128 >= n
	
	getxytable(y1, x1, &n);
	double yp0 = 0.0, ypn_1 = 0.0;
	double ypp[n];
	
	printf("Please choose cubic spline type\n");
	printf("n: for natural cubic spline\n");
	printf("c: for clampled cubic spline\n");
	printf("Default: natural cubic spline\n");
	char c = getch();
	switch(c)
	{
		case 'n':
		case 'N':
			printf("\nNatural cubic spline case\n");
			ncspline(x1, y1, n, n, n, ypp);
			break; // case 'n'
		case 'c':
		case 'C':
			printf("\nClamped cubic spline case\n");
			clcspline(x1, y1, n, n, n, yp0, ypn_1, ypp);
			break; // case 'c'
		default:
			printf("\nNatural cubic spline case\n");
			ncspline(x1, y1, n, n, n, ypp);
			break;
	}
	
	double x[n], y[n];
	int i;
	for (i = 0; i < n; i++)
	{
		x[i] = x1[i];
		y[i] = y1[i];
	}
	
	double xi, fxi = 0.0, fxi_1 = y[0], fpxi = 0.0, fpxi_1 = 0.0, fppxi = 0.0;
	double A = (x[n-1] - x[0]);
	i = 0;
	int m = 100;
	
	printf("\n%9s\t %9s\t %9s\t %9s\n\n", "xi", "fxi", "fpxi", "fppxi");
	
	while (i <= m)
	{
		double f = i / ((double) m);
		/* xi is a s-curve */
		xi = (3.0 - 2.0 * f) * f * f * A + x[0];
		csplint(x, y, ypp, n, xi, &fxi);
		/* calculate the first derivative of fx */
		fpxi = (fxi - fxi_1) / (1.0 / ((double)m));
		/* calculate the second derivative of fx */
		fppxi = (fpxi - fpxi_1) / (1.0 / ((double)m));
		printf("%9.3f\t %9.3f\t %9.3f\t %9.3f\n", xi, fxi, fpxi, fppxi);
		//printf("%9.3f\n", fpxi);
		i++;
		fxi_1 = fxi;
		fpxi_1 = fpxi;
	}
	
	return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
