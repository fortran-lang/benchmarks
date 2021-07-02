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
void swap(int M, double a[M][M], double b[M][M]) {
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
double mmax(int M, double a[M][M], double b[M][M]) {
  double max = 0.0;
  double diff = 0.0;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      diff = fabs(a[i][j]-b[i][j]);
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

void run(int M, int N_ITER) {

  double phi[M][M];
  double phip[M][M]; 
  double rhoa[M][M]; 

  clock_t tic = clock();
  int iter = 0;

  // initialize matrices
  memset(&phip[0][0], 0, sizeof(double)*M*M);
  memset(&phi[0][0], 0, sizeof(double)*M*M);
  for (int i=1; i<M-1; i++)  {
    for (int j=1; j<M-1; j++) {
      rhoa[i][j] =  rho(i*SPACE,j*SPACE);
    }
  }

  double delta;
  double a2 = pow(SPACE,2.0);
  
  while (iter < N_ITER ) {
    iter += 1;
    
    for (int i=1; i < M-1; i++) {
      for (int j=1; j < M-1; j++) {
        phip[i][j] = (phi[i+1][j] + phi[i-1][j] +
                      phi[i][j+1] + phi[i][j-1])/4.0 +
	  a2/(4.0 * EPSILON0)*rhoa[i][j];
      }
    }
    delta = mmax(M, phi, phip);
    swap(M, phi, phip);
    
  }

  clock_t toc = clock();

  write_timing(M, iter, delta, toc - tic);
}


int main(int argc, char *argv[]) {

  int N, N_ITER;

  // Retrieve arguments
  parse_args(argc, argv, &N, &N_ITER);

  run(N, N_ITER);
}
