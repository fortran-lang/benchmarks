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
| Python (pure)       | 276             | n/a                   | n/a                    |
| Cython              | 1.02            | 32.8                  | 229                    |
| Fortran (naive)     | 0.34            | 13.2                  | 69.7                   |
| Fortran (optimized) | 0.18            | 6.25                  | 31.4                   |
| C (naive)           | 0.42*           | 7.25                  | 33.7                   |
| C (optimized)       | 0.37*           | 6.80                  | 32.8                   |

* For all of these results, the amount of iterations performed by the respective codes was approximately
the same, with the exception of the 100x100 grid C codes, which did nearly a double amount of iterations
compared to the rest, for some reason. The timings are on AMD Ryzen 5 3600 @3.6GHz, using WSL2 on Windows 10.

Some thoughts on the code at https://runningcrocodile.fi/articles/pythonperformance.html . Also, there was
discussion about this problem at Discourse: https://fortran-lang.discourse.group/t/performance-c-vs-fortran/1461 
with important contributions from numerous users (Mason, Beliavsky, septc, implicitall, nncarlson, han190 and pmk)
