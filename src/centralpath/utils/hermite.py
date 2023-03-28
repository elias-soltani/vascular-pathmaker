"""
Utility functions to fit a cubic hermite curve to data.
"""

import numpy as np
from scipy.optimize import minimize


def cubic_hermite_function(xi):
    """
    one-dimensional cubic Hermite basis functions.
    x = f1*x1 + f2*x2 + f3*dx1/dxi + f4*dx2/dxi
    f1: psi^0_1, f2: psi^0_2, f3: psi^1_1, f4: psi^1_2
    :param xi: parameter
    :return: two numpy arrays including the functions and derivatives
    """
    f1 = 1 - 3 * xi**2 + 2 * xi**3
    f2 = xi**2 * (3 - 2 * xi)
    f3 = xi * (xi - 1)**2
    f4 = xi**2 * (xi - 1)

    f1p = -6 * xi + 6 * xi**2
    f2p = 2 * xi * (3 - 2 * xi) - 2 * xi**2
    f3p = (xi - 1)**2 + 2 * xi * (xi - 1)
    f4p = 2 * xi * (xi - 1) + xi**2

    f = np.array([f1, f2, f3, f4])
    fp = np.array([f1p, f2p, f3p, f4p])
    return f, fp


def interpolate_cubic_hermite(parameters, xi):
    """
    interpolate a field at xi. dx1 is dx/dxi at xi = 0 and dx2 is dx/dxi at xi=1
    :param parameters: a list of x1, x2, dx1, dx2 in form of numpy array.
    :param xi: xi parameter along the curve. ndarray, shape(n,)
    :return: x, dx coordinates evaluated at xi. ndarray, shape(n,3)
    """
    x1, x2, dx1, dx2 = parameters
    f, fp = cubic_hermite_function(xi)
    f1, f2, f3, f4 = f
    f1p, f2p, f3p, f4p = fp
    x = f1.reshape(-1, 1)*x1.reshape(1, 3) + f2.reshape(-1, 1)*x2.reshape(1, 3) + \
        f3.reshape(-1, 1)*dx1.reshape(1, 3) + f4.reshape(-1, 1)*dx2.reshape(1, 3)
    dx = f1p.reshape(-1, 1)*x1.reshape(1, 3) + f2p.reshape(-1, 1)*x2.reshape(1, 3) + \
         f3p.reshape(-1, 1)*dx1.reshape(1, 3) + f4p.reshape(-1, 1)*dx2.reshape(1, 3)

    return x, dx


def cost_function(params, xi, p1, p2, p):
    """
    calculate loss function. Unknown parameters are dx1 and dx2.
    :param params: params=[dx1, dy1, dz1, dx2, dy2, dz2]
    :param xi: xi parameter. ndarray, shape (n,)
    :param p1: point 1 coordinates
    :param p2: point 2 coordinates
    :param p: given data points. ndarray, shape (n,3) where n is the number of points.
    :return: cost function or objective function.
    """

    dx1 = params[:3]
    dx2 = params[3:]
    xpred, dxpred = interpolate_cubic_hermite(np.array([p1, p2, dx1, dx2]), xi)

    n = xi.size  # the number of data points.

    J = np.sum((xpred - p)**2)/(2*n)

    return J


def jacobian(params, xi, p1, p2, p):
    """
    calculate jacobian function. Unknown parameters are dx1 and dx2.
    :param params: params=[dx1, dy1, dz1, dx2, dy2, dz2]
    :param xi: xi parameter. ndarray, shape (n,)
    :param p1: point 1 coordinates
    :param p2: point 2 coordinates
    :param p: given data points. ndarray, shape (n,3) where n is the number of points.
    :return: jacobian. ndarray, shape (6,)
    """
    jac = np.empty((params.size))
    f, fp = cubic_hermite_function(xi)
    f1, f2, f3, f4 = f

    n = xi.size

    dx1 = params[:3]
    dx2 = params[3:]
    xpred, dxpred = interpolate_cubic_hermite(np.array([p1, p2, dx1, dx2]), xi)

    jac[0] = np.sum((xpred[:, 0] - p[:, 0])*f3)/n
    jac[1] = np.sum((xpred[:, 1] - p[:, 1])*f3)/n
    jac[2] = np.sum((xpred[:, 2] - p[:, 2])*f3)/n
    jac[3] = np.sum((xpred[:, 0] - p[:, 0])*f4)/n
    jac[4] = np.sum((xpred[:, 1] - p[:, 1])*f4)/n
    jac[5] = np.sum((xpred[:, 2] - p[:, 2])*f4)/n

    return jac


def compute_derivatives(params, args, disp=True):
    """
    compute derivatives at the end points of the curve.
    :param params: params=[dx1, dy1, dz1, dx2, dy2, dz2]
    :param args: (xi, p1, p2, p). where xi is ndarray, shape (n,). p1 and p2 two endpoints, shape(3,). p represents
        the points at different xi. shape (n,3)
    :param disp: If true, display information about optimization.
    :return: d1 at p1 and p2. or dx/dxi at p1 and dx/dxi at p2.
    """
    res = minimize(cost_function, params, method='BFGS', jac=jacobian, args=args,
                   options={'disp': disp})

    d1_1 = res.x[:3]
    d1_2 = res.x[3:]

    return d1_1, d1_2


def compute_curve_length(points):
    """
    Calculate the curve length
    :param points: ndarray, shape(n,3) where n is number of points
    :return: curve length.
    """
    return np.sum(np.linalg.norm(points[1:, :] - points[:-1, :], axis=1))


def compute_xi_from_point(points):
    """
    Calculate xi for a point with index as point_idx and all the curve points.
    :param points: curve points. ndarray, shape(n,3) where n is number of points
    :return: xi.
    """
    curve_length = compute_curve_length(points)
    if curve_length == 0:
        raise 'the curve length should not be zero'

    n, _ = points.shape
    xi = np.zeros(n)
    for idx in range(n):
        if idx == 0:
            xi[idx] = 0
        else:
            xi[idx] = np.sum(np.linalg.norm(points[1:idx+1, :] - points[:idx, :], axis=1))/curve_length

    return xi


def initialize_d1(points):
    """
    Calculate d1_1 and d1_2 approximately.
    :param points: curve points. ndarray, shape(n,3) where n is number of points
    :return: d1_1 and d1_2
    """
    d1_1 = points[1, :] - points[0, :] / np.linalg.norm(points[1, :] - points[0, :])
    d1_2 = points[-1, :] - points[-2, :] / np.linalg.norm(points[-1, :] - points[-2, :])

    return d1_1, d1_2


def fit_compute_derivatives(curve_points, disp=True, plot=False):
    """

    :param curve_points: ndarray, shape(n,3) where n is the number of points
    :param disp: If true it outputs optimization information.
    :param plot: If true it plots the fitted curve and data
    :return: parameters for cubic Hermite spline interpolation. x1, x2, dx1 and dx2.
    """
    xi_delta = 0.01
    # get xi for each point on the curve
    xi = compute_xi_from_point(curve_points)
    # first guess for the d1_1 and d1_2
    dxp1, dxp2 = initialize_d1(curve_points)
    params = np.hstack((dxp1, dxp2))
    # if we have only two points then it's a line
    if len(curve_points) == 2:
        d1_1 = curve_points[-1] - curve_points[0]
        d1_2 = curve_points[-1] - curve_points[0]
    else:
        d1_1, d1_2 = compute_derivatives(params, args=(xi, curve_points[0], curve_points[-1], curve_points), disp=disp)

    xi_test = np.arange(0, 1.0, xi_delta)
    if plot:
        x, dx = interpolate_cubic_hermite(np.array([curve_points[0], curve_points[-1], d1_1, d1_2]), xi_test)
        plot_curve(curve_points, x)

    return curve_points[0], curve_points[-1], d1_1, d1_2


def plot_curve(points, x):
    """
    Plot fitted cubic Hermite spline with data
    :param points: data points
    :param x: interpolated x values along the curve.
    """
    import matplotlib.pyplot as plt

    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(points[:, 0], points[:, 1], points[:, 2], 'r.', label='data points')
    ax.plot(x[:, 0], x[:, 1], x[:, 2], label='cubic Hermite spline', linewidth=3, linestyle='-')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('fitting 1D cubic Hermite spline to data points')
    ax.set_box_aspect((np.ptp(x[:, 0]), np.ptp(x[:, 1]), np.ptp(x[:, 2])))
    plt.legend()
    plt.show()


