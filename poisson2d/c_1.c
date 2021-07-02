#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "c_common.h"

// COMMENT
// I know, this is extremely bad practice.
// But to accommodate argument sizes for the arguments
// I felt this was ok.
// It makes it much more versatile when testing
// END COMMENT
void swap(int M, double (*a)[M], double (*b)[M])
{
  double swp;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      swp = a[i][j];
      a[i][j] = b[i][j];
      b[i][j] = swp;
    }
  }
}

// See comment above
double mmax(int M, double (*phi)[M], double (*phip)[M]) {
  double max = 0.0;
  double diff = 0.0;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      diff = fabs(phi[i][j]-phip[i][j]);
      if (diff > max)
        max = diff;
    }
  }
  return max;
}

double rho(double x, double y) {
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

clock_t run(int M, int N_ITER) {

  // COMMENT:
  // I wouldn't even imagine any novice in C would
  // do something like this. It really makes no sense :)
  // END OF COMMENT
  double (*phi)[M];
  double (*phip)[M]; 
  double (*rhoa)[M]; 
  phi = malloc(sizeof(double[M][M]));
  phip = malloc(sizeof(double[M][M]));
  rhoa = malloc(sizeof(double[M][M]));

  clock_t tic = clock();

  // initialize matrices
  memset(phip, 0, sizeof(phip[0][0])*M*M);
  memset(phi, 0, sizeof(phi[0][0])*M*M);
  
  double delta = 1.0;
  double a2 = pow(SPACE,2.0);
  
  int iter = 0;
  while (iter < N_ITER ) {
    iter += 1;
    
    for (int i=1; i < M-1; i++) {
      for (int j=1; j < M-1; j++) {
        phip[i][j] = (phi[i+1][j] + phi[i-1][j] +
                      phi[i][j+1] + phi[i][j-1])/4.0 +
	  a2/(4.0 * EPSILON0)*rho(i*SPACE,j*SPACE);
      }
    }
    delta = mmax(M, phi, phip);
    swap(M, phi, phip);
    
  }

  clock_t toc = clock();

  free(phi);
  free(phip);
  free(rhoa);

  write_timing(M, iter, delta, toc - tic);
}


int main(int argc, char *argv[]) {

  int N, N_ITER;

  // Retrieve arguments
  parse_args(argc, argv, &N, &N_ITER);
  run(N, N_ITER);
}
