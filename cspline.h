#ifndef CSPLINE_H
#define CSPLINE_H
#define PI		3.1416
/*----------------------------------------------------------------------------*/
#ifdef __cplusplus      // avoid c++ name mangling
  extern "C" {
#endif
void ncspline(double x[], double y[], int nx, int ny, int n, double ypp[]);
void clcspline(double x[], double y[], int nx, int ny, int n, double yp0, double yn_1, double ypp[]);
void tridiag(double x[], int N, double a[], double b[], double c[]);
void getxytable(double x[], double y[], int *n);
void csplint(double x[], double y[], double ypp[], int n, double xi, double *yi);
#ifdef __cplusplus
  }
#endif
/*----------------------------------------------------------------------------*/
#endif //CSPLINE_H
