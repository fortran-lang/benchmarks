#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void parse_args(int argc, char **argv, int *N, int *N_iter) {
  if ( argc < 3 ) {
    printf("Not enough arguments [N N_ITER]\n");
  }

  *N = atoi(argv[1]);
  *N_iter = atoi(argv[2]);
}

void write_timing(int N, int iter, double delta, clock_t timing) {
  printf("%16d %16d %22.15e %22.15e\n", N, iter, delta, ((double)timing)/CLOCKS_PER_SEC/iter);
}
