Python-fast-cubic-spline
========================

Implementation of 1D and 2D fast spline interpolation algorithm ([Habermann
and Kindermann 2007] [1]) in Python.

Calculation of spline coefficients are in NumPy, and actual interpolation
routine is coded in Cython. This is advantageous since if your main routine is
coded in Cython, once you have coefficients, you can call interpolation
functions without any Python overhead.

Documentation
-------------

Compile cython module:
```sh
$ python setup.py build_ext --inplace
```

Run the main module for an example:
```sh
$ python fast_cubic_spline.py
```


References
-------------

Habermann, C., & Kindermann, F. (2007). Multidimensional Spline Interpolation:
Theory and Applications. Computational Economics, 30(2), 153–169.

[1]: http://www.springerlink.com/index/10.1007/s10614-007-9092-4
"Habermann, C., & Kindermann, F. (2007). Multidimensional Spline
Interpolation: Theory and Applications. Computational Economics, 30(2),
153–169."
