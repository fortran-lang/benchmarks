## Poisson solver in various languages

This solver uses a simple Jacobi iteration to solve a Poisson equation. For details, see, for example,
"Computational Physics" by Mark Newman, Chap 9. Compiler commands were:

```gfortran sourcefile.f95 -O3 -o program```

```gcc sourcefile.c -O3 -o program```

For Cython module (import it to run):

```python setup.py build_ext --inplace```


There are various implementations using various strategies.
All codes are named with a 1-1 correspondance, i.e. `c_i_*` has the same strategy as `fortran_i_*` files.
The suffixed `_*` is only used if there are multiple alternatives that utilizes the same strategy between
the different languages. For example the `c_4_1d` and `c_4_2d` are both using the same strategy for
solving the Poisson equation, while the former uses a 1D array for the data while the latter uses a 2D
array with some slight overhead of another pointer level.

The two files `c_common.c` and `fortran_common.f90` are baseline programs used for printing
out data and for parsing the command line arguments.

All executables accept 2 arguments:

    # N = 100
    # N_ITER = 2000
    ./c_1 100 2000
    100             2000  5.658391137993336e+04  1.448613750000000e-03

which will use a 2D array of size `100x100` and run for 2000 iterations.
The output line consists of 4 values, 1) N, 2) N_ITER, 3) max-difference, 4) time per iteration.

Some thoughts on the code at https://runningcrocodile.fi/articles/pythonperformance.html. Also, there was
discussion about this problem at Discourse: https://fortran-lang.discourse.group/t/performance-c-vs-fortran/1461 
with important contributions from numerous users (Mason, Beliavsky, septc, implicitall, nncarlson, han190 and pmk).

The entire code base has been re-written to allow command-line arguments allowing to sample different
matrix sizes with the same executable (user zerothi).

To run a preset benchmark suite, simply execute `bash run.sh` which will run the currently
implemented benchmarks with dimensions up to 600.


## Old timings

The grid is a 2-dimensional array of size MxM. The timings for different values of M are, in seconds:

| Language            | M=100           | M=200                 | M=300                  |
|---------------------|-----------------|-----------------------|------------------------|
| Python (pure)       | 276             | n/a                   | n/a                    |
| Cython              | 1.02            | 32.8                  | 229                    |
| Fortran (naive)     | 0.34            | 13.2                  | 69.7                   |
| Fortran (optimized) | 0.18            | 6.25                  | 31.4                   |
| C (naive)           | 0.42*           | 7.25                  | 33.7                   |
| C (optimized)       | 0.37*           | 6.80                  | 32.8                   |

* For all of these results, the amount of iterations performed by the respective codes was approximately
the same, with the exception of the 100x100 grid C codes, which did nearly a double amount of iterations
compared to the rest, for some reason. The timings are on AMD Ryzen 5 3600 @3.6GHz, using WSL2 on Windows 10.
