#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "c_common.h"

double ** malloc_2d(int m, int n) {
  if (m <= 0 || n <= 0)
    return NULL;

  // Allocate pointers for columns
  double **array2D = malloc(m * sizeof(double *));
  if (array2D == NULL) {
    return NULL;
  }

  // Allocate actual data
  array2D[0] = malloc(m * n * sizeof(double));
  if (array2D[0] == NULL) {
    free(array2D);
    return NULL;
  }

  // assign pointers to make it 2D
  for(int i = 1; i < m; i++) {
    array2D[i] = array2D[0] + i * n ;
  }

  return array2D;
}

void free_2d(double **array2D) {
  free(array2D[0]);
  free(array2D);
}



double rho(const double x, const double y) {
  const double s1 = 0.6;
  const double e1 = 0.8;
  const double s2 = 0.2;
  const double e2 = 0.4;

  if (s1 < x && x < e1 && s1 < y && y < e1) {
    return 1.0;
  } else if ( s2 < x && x < e2 && s2 < y && y < e2 ) {
    return -1.0;
  } else {
    return 0.0;
  }
}

double iterate(int M, double *restrict*restrict phi, double *restrict*restrict phinew, double *restrict*restrict rhoa) {
  double delta = 0, err;
  for (int i=1; i < M-1; i++) {
    for (int j=1; j < M-1; j++) {
      phinew[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] + rhoa[i][j])*0.25;
      err = fabs(phinew[i][j] - phi[i][j]);
      if ( err > delta ) delta = err;
    }
  }
  return delta;
}

void init_rho(int M, double *restrict*restrict rhoa, const double epsilon, const double a) {
  const double a2 = a * a / epsilon;
  for (int i=1; i<M-1; i++)  {
    for (int j=1; j<M-1; j++) {
      rhoa[i][j] = rho(i*a,j*a) * a2;
    }
  }
}


void run(int M, int N_ITER) {

  double **phi;
  double **phip;
  double **rhoa;
  double **tmp;

  // A real world program will definitely use malloc
  phi = malloc_2d(M, M);
  phip = malloc_2d(M, M);
  rhoa = malloc_2d(M, M);

  clock_t tic = clock();

  memset(&phi[0][0], 0, sizeof(double)*M*M);
  memset(&phip[0][0], 0, sizeof(double)*M*M);

  // In C one tries to avoid using pow because
  // it assumes floating point powers (integers are faster)
  // So better do it directly
  init_rho(M, rhoa, EPSILON0, SPACE);
  
  int iter = 0;
  double delta;
  while ( iter < N_ITER ) {
    iter += 1;

    delta = iterate(M, phi, phip, rhoa);

    // swap pointers (no copies)
    tmp = phi;
    phi = phip;
    phip = tmp;
  }

  clock_t toc = clock();

  free_2d(phi);
  free_2d(phip);
  free_2d(rhoa);

  write_timing(M, iter, delta, toc - tic);
}

int main(int argc, char *argv[]) {
  int N, N_ITER;

  // Retrieve arguments
  parse_args(argc, argv, &N, &N_ITER);
  run(N, N_ITER);
}
