# Fast-Cubic-Spline-Python provides an implementation of 1D and 2D fast spline
# interpolation algorithm (Habermann and Kindermann 2007) in Python.
# Copyright (C) 2012, 2013 Joon H. Ro

# This file is part of Fast-Cubic-Spline-Python.

# Fast-Cubic-Spline-Python is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Fast-Cubic-Spline-Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#cython: boundscheck = False
#cython: wraparound = False
#cython: cdivision = True

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport fabs, fmin
from cython.parallel import prange

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.

'''
Cubic spline interpolation using Habermann and Kindermann (2007)'s algorithm 
'''
#----------------------------------------------------------------------
cdef double Pi(double t) nogil:
    cdef:
        double abs_t = fabs(t)
    if abs_t <= 1:
        return(4 - 6 * abs_t**2 + 3 * abs_t **3)
    elif abs_t <= 2:
        return((2 - abs_t)**3)
    else:
        return(0)

#----------------------------------------------------------------------        
cdef double u(double x, int k, double a, double h) nogil:
    return(Pi((x - a)/h - (k - 2)))

#----------------------------------------------------------------------
def interpolate(double x,
                double a, double b,
                double[:] c,
                ):
    '''
    Return interpolated function value at x
    
    Parameters
    ----------
    x : float
        The value where the function will be approximated at
    a : double
        Lower bound of the grid
    b : double
        Upper bound of the grid
    c : ndarray
        Coefficients of spline
    
    Returns
    -------
    out : float
        Approximated function value at x
    '''
    return(_interpolate(x, a, b, c))

cdef double _interpolate(double x,
                         double a, double b,
                         double[:] c,
                         ) nogil:

    cdef:
        int n = c.shape[0] - 3
        double h = (b - a)/n
        int l = <int>((x - a)//h) + 1
        int m = <int>(fmin(l + 3, n + 3))
        int i1
        double s = 0

    for i1 in xrange(l, m + 1):
        s += c[i1 - 1] * u(x, i1, a, h)

    return(s)

#----------------------------------------------------------------------
def interpolate_2d(double x, double y,
                   double a1, double b1,
                   double a2, double b2,
                   double[:, :] c,
                   ):
    '''
    Return interpolated function value at x
    
    Parameters
    ----------
    x, y : float
        The values where the function will be approximated at
    a1, b1 : double
        Lower and upper bounds of the grid for x
    a2, b2 : double
        Lower and upper bounds of the grid for y
    c : ndarray
        Coefficients of spline
    
    Returns
    -------
    out : float
        Approximated function value at (x, y)
    '''
    return(_interpolate_2d(x, y, a1, b1, a2, b2, c))
    
cdef double _interpolate_2d(double x, double y,
                            double a1, double b1,
                            double a2, double b2,
                            double[:, :] c,
                            ) nogil:
    cdef:
        int n1 = c.shape[0] - 3
        int n2 = c.shape[1] - 3
        double h1 = (b1 - a1)/n1
        double h2 = (b2 - a2)/n2
        int l1 = <int>((x - a1)//h1) + 1
        int l2 = <int>((y - a2)//h2) + 1
        int m1 = <int>(fmin(l1 + 3, n1 + 3))
        int m2 = <int>(fmin(l2 + 3, n2 + 3))
        int i1, i2
        double s = 0
        double u_x, u_y

    for i1 in xrange(l1, m1 + 1):
        u_x = u(x, i1, a1, h1)
        for i2 in xrange(l2, m2 + 1):
            u_y = u(y, i2, a2, h2)
            s += c[i1 - 1, i2 - 1] * u_x * u_y

    return(s)
