#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define M 300
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
double rho(double x, double y)
{
  double s1 = 0.6;
  double e1 = 0.8;
	double s2 = 0.2;
	double e2 = 0.4;

	if (x > s1 && x < e1 && y > s1 && y < e1) {
		return 1.0;
	} else if (x > s2 && x < e2 && y > s2 && y < e2 ) {
		return -1.0;
	} else {
		return 0.0;
	}
}
void run(double toler, double a)
{
  double epsilon0 = 8.85e-12;
  double a2;

  double (*phi)[M];
  double (*phip)[M]; 
  double (*rhoa)[M]; 
  phi = malloc(sizeof(double[M][M]));
  phip = malloc(sizeof(double[M][M]));
  rhoa = malloc(sizeof(double[M][M]));

  int iter = 0;

  memset(phip, 0, sizeof(phip[0][0])*M*M);
  memset(phi, 0, sizeof(phi[0][0])*M*M);
  for (int i=1; i<M-1; i++)  {
    for (int j=1; j<M-1; j++) {
      rhoa[i][j] =  rho(i*a,j*a);
    }
  }
  double delta = 1.0;
  a2 = pow(a,2.0);
  while (delta > toler) {
    iter += 1;
    
    for (int i=1; i < M-1; i++) {
      for (int j=1; j < M-1; j++) {
        phip[i][j] = (phi[i+1][j] + phi[i-1][j] +
                      phi[i][j+1] + phi[i][j-1])/4.0 +
                      a2/(4.0 * epsilon0)*rhoa[i][j];
      }
    }
    delta = mmax(phi, phip, M);
    swap(phi, phip, M);

  }
  printf("iters %d", iter);

}

int main(int argc, char *argv[])
{
  clock_t start = clock();
  run(1e-6, 0.01);
  clock_t end = clock();
  double total = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Execution time: %f\n",total);
}
