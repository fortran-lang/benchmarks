#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define M 200
#define analytical_tol 5.0e-3

/*
 *
 * Code by ARC, @arunningcroc
 *
 */
void swap(double (*a)[M], double (*b)[M], int n)
{
  double swp;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      a[i][j] = b[i][j];
    }
  }
}
double mmax(double (*phi)[M], double (*phip)[M], int n)
{
  double max = 0.0;
  double diff = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      diff = fabs(phi[i][j]-phip[i][j]);
      if (diff > max)
        max = diff;
    }
  }
  return max;
}
void normalize(double (*phi)[M]) {
  double max = -50000;
  for(int i=0; i<M; i++) {
    for(int j=0; j<M; j++) {
      if(phi[i][j] > max) max = phi[i][j];
    }
  }
  for(int i=0; i<M; i++) {
    for(int j=0; j<M; j++) {
      phi[i][j] /= max;
    }
  }
  
}
double get_average_error(double (*sol)[M], double (*analytical)[M]) {
  double sum = 0.0;
  normalize(sol);
  normalize(analytical);
  for(int i=0; i<M; i++) {
    for(int j=0; j<M; j++) {
      sum += fabs(sol[i][j]-analytical[i][j]);
    }
  }
  return sum/(M*M);
}
double rho(double x, double y)
{
  return -(6.0*x*y*(1.0-y)-2.0*pow(x,3.0));
}
double boundary(double y)
{
  return y*(1.0-y);
}
double analytical_sol(double x,double y)
{
  return y*(1.0-y)*pow(x,3.0);
}
int run(double toler, double a)
{
  double a2;

  double (*phi)[M];
  double (*phip)[M]; 
  double (*rhoa)[M]; 
  double (*analytical)[M];
  phi = malloc(sizeof(double[M][M]));
  phip = malloc(sizeof(double[M][M]));
  rhoa = malloc(sizeof(double[M][M]));
  analytical = malloc(sizeof(double[M][M]));
  for(int i=0; i<M; i++) {
    for (int j=0; j<M; j++) {
      analytical[i][j] = analytical_sol(i*a,j*a);
    }
  }
  int iter = 0;

  memset(phip, 0, sizeof(phip[0][0])*M*M);
  memset(phi, 0, sizeof(phi[0][0])*M*M);
  for (int i=1; i<M-1; i++)  {
    for (int j=1; j<M-1; j++) {
      rhoa[i][j] =  rho(i*a,j*a);
    }
  }
  for (int j=0; j < M; j++) {
    phip[M-1][j] = j*a*(1-j*a);
    phi[M-1][j] = j*a*(1-j*a);
  }
  double delta = 1.0;
  a2 = pow(a,2.0);
  while (delta > toler) {
    iter += 1;
    
    for (int i=1; i < M-1; i++) {
      for (int j=1; j < M-1; j++) {
        phip[i][j] = (phi[i+1][j] + phi[i-1][j] +
                      phi[i][j+1] + phi[i][j-1])/4.0 +
                      a2*rhoa[i][j]/4.0;
      }
    }
    delta = mmax(phi, phip, M);
    swap(phi, phip, M);

  }
  if(get_average_error(phi, analytical) > analytical_tol) {
    printf("Failed to reach analytical solution, error %f\n", get_average_error(phi, analytical) );
  }
  return iter;

}

int main(int argc, char *argv[])
{
  int iter;
  clock_t start = clock();
  iter = run(1e-8, 1.0/(M-1));
  clock_t end = clock();
  double total = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Execution time %f, iters: %d \n",total, iter);
}
