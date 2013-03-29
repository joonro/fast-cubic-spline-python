************
Introduction
************

:Date: Mar 28, 2013
:Version: 0.1.0
:Authors: Joon Ro, joonhyoung.ro[at]gmail.com
:Web site: https://github.com/joonro/fast-cubic-spline-python
:Copyright: This document has been placed in the public domain.
:License: Fast-Cubic-Spline-Python is released under the GPLv3.


Purpose
=======

Fast-Cubic-Spline-Python provides an implementation of 1D and 2D fast spline
interpolation algorithm `Habermann and Kindermann 2007 <http://www.springerlink.com/index/10.1007/s10614-007-9092-4 "Habermann, C., & Kindermann, F. (2007). Multidimensional Spline Interpolation: Theory and Applications. Computational Economics, 30(2), 153–169.">`_) in Python.
Higher dimensional interpolation is possible with this code, but only 1D and
2D examples are provided.

Calculation of spline coefficients are in NumPy, and actual interpolation
routine is coded in Cython. This is advantageous since if your main routine is
coded in Cython, once you have coefficients, you can call interpolation
functions without any Python overhead.

*****
Usage
*****

Run the main module for an example:

.. code-block:: sh

    $ python fast_cubic_spline.py


************
Installation
************

Dependencies
============

* Python
* Cython (http://cython.org)

Compiling Cython Module
=======================

.. code-block:: sh

    $ python setup.py build_ext --inplace

References
==========

Habermann, C., & Kindermann, F. (2007). Multidimensional Spline Interpolation:
Theory and Applications. Computational Economics, 30(2), 153–169.

