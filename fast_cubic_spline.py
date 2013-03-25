from __future__ import division

from numpy import arange, array, linspace, ones, zeros
from scipy.linalg import solve_banded

import fast_cubic_spline as _spline

def cal_coefs(a, b, y, c=None, alpha=0, beta=0):
    '''
    Return spline coefficients 

    Parameters
    ----------
    a : float
        lower bound of the grid.
    b : float
        upper bound of the grid.
    y : ndarray
        actual function value at grid points.
    c : (y.shape[0] + 2, ) ndarray, optional
        ndarry to be written
    alpha : float
        Second-order derivative at a. Default is 0.
    beta : float
        Second-order derivative at b. Default is 0.

    Returns
    -------
    out : ndarray
        Array of coefficients.
    '''
    n = y.shape[0] - 1
    h = (b - a)/n

    if c is None:
        c = zeros((n + 3, ))
        ifreturn = True
    else:
        assert(c.shape[0] == n + 3)
        ifreturn = False

    c[1] = 1/6 * (y[0] - (alpha * h**2)/6)
    c[n + 1] = 1/6 * (y[n] - (beta * h**2)/6)

    # ab matrix here is just compressed banded matrix
    ab = ones((3, n - 1))
    ab[0, 0] = 0
    ab[1, :] = 4
    ab[-1, -1] = 0

    B = y[1:-1].copy()
    B[0] -= c[1]
    B[-1] -=  c[n + 1]

    c[2:-2] = solve_banded((1, 1), ab, B)

    c[0] = alpha * h**2/6 + 2 * c[1] - c[2]
    c[-1] = beta * h**2/6 + 2 * c[-2] - c[-3]

    if ifreturn:
        return(c)

# aliases
interpolate = _spline.interpolate
interpolate_2d = _spline.interpolate_2d

if __name__ == '__main__':

    from matplotlib import pyplot as plt

    f = lambda x: x**2

    a = -1
    b = 1
    n = 49  # n + 1 grid points (0,..., n)

    h = (b - a)/n
    grid = arange(n + 1) * h + a

    y = f(grid)
    alpha = 0
    beta = 0

    c = cal_coefs(a, b, y)

    grid_hat = linspace(a, b, 100)
    fhat = array([interpolate(x, a, b, c) for x in grid_hat])

    plt.plot(grid_hat, fhat)
    plt.plot(grid_hat, f(grid_hat))
    plt.show()

    # 2D interpolation
    f2d = lambda x, z: x**2 + 2*x + 1 + z ** 0.5 + 3 * z
    
    a1, a2 = 0, 0
    b1, b2 = 1, 1
    n1, n2 = 49, 39  # n + 1 grid points (0,..., n)

    h1, h2 = (b1 - a1)/n1, (b2 - a2)/n2
    grid_x = arange(n1 + 1) * h1 + a1
    grid_z = arange(n2 + 1) * h2 + a2

    y = zeros((n1 + 1, n2 + 1))

    for i, x in enumerate(grid_x):
        for j, z in enumerate(grid_z):
            y[i, j] = f2d(x, z)

    alpha = 0
    beta = 0

    c_tmp = zeros((n1 + 3, n2 + 1))
    cal_coefs(a1, b1, y, c_tmp)

    c = zeros((n1 + 3, n2 + 3))
    # NOTE: here you have to pass c_tmp.T and c.T
    cal_coefs(a2, b2, c_tmp.T, c.T)

    fhat = zeros((n1 + 1, n2 + 1))
    for i, x in enumerate(grid_x):
        for j, z in enumerate(grid_z):
            fhat[i, j] = interpolate_2d(x, z, a1, b1, a2, b2, c)

    real_val = zeros((n1 + 1, n2 + 1))
    for i, x in enumerate(grid_x):
        for j, z in enumerate(grid_z):
            real_val[i, j] = f2d(x, z)

    from mayavi import mlab

    def draw_3d(grid_x, grid_y, fval, title='pi'):
        mlab.figure()
        mlab.surf(grid_x, grid_y, fval)#, warp_scale="auto")
        mlab.axes(xlabel='x', ylabel='z', zlabel=title)
        mlab.orientation_axes(xlabel='x', ylabel='z', zlabel=title)
        mlab.title(title)

    draw_3d(grid_x, grid_z, fhat, title='interpolated')
    draw_3d(grid_x, grid_z, real_val, title='real')
    mlab.show()
