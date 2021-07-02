## Poisson solver in various languages

This solver uses a simple Jacobi iteration to solve a Poisson equation. For details, see, for example,
"Computational Physics" by Mark Newman, Chap 9. Compiler commands were:

```gfortran sourcefile.f95 -Ofast -o program```

```gcc sourcefile.c -Ofast -o program```

For Cython module (import it to run):

```python setup.py build_ext --inplace```

The grid is a 2-dimensional array of size MxM. The timings for different values of M are, in seconds:

| Language            | M=100           | M=200                 | M=300                  |
|---------------------|-----------------|-----------------------|------------------------|
| Python (pure)       | 198             | n/a                   | n/a                    |
| Cython              |	0.46            | 5.77                  | 43.9                   |
| Fortran (trad1)     | 0.23            | 3.59                  | 18.3                   |
| Fortran (trad2)     | 0.19            | 4.29                  | 19.6                   |
| C (naive)           | 0.18            | 2.46                  | 11.6                   |
| C (optimized)       | 0.13            | 1.83                  | 8.26                   |
| Fortran (modern)    | 0.19            | n/a                   | 8.08                   |
| Python (vectorized) | 0.92            | 10.8                  | 65.7                   |

* For all of these results, the amount of iterations performed was exactly the same.
 The timings are on AMD Ryzen 5 3600 @3.6GHz, using WSL2 on Windows 10.

Some thoughts on the code at https://runningcrocodile.fi/articles/pythonperformance.html . Also, there was
discussion about this problem at Discourse: https://fortran-lang.discourse.group/t/performance-c-vs-fortran/1461 
with important contributions from numerous users (Mason, Beliavsky, septc, implicitall, nncarlson, han190 and pmk)

Codes written by a running crocodile (traditional1.f90, traditional2.f90, naive.c, optimized.c, poisson.py, poissonvectorized.py, poisson.pyx)
and Damian Rouson (modern.f90) with the aforementioned help from Discourse.
