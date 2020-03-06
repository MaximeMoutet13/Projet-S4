__author__ = "maxime"
__file__ = "Implicite.py"
__date__ = "10/02/20"

import matplotlib.pyplot as plt
import numpy as np

A = 201.4
B = 1.45
D = 0.310
C = 14 / B
Q = 340  # W.m-2


def S(x):
    """
    Mean annual distribution of radiation at each latitude
    """
    S2 = -0.482
    return 1 + S2 * (3 * x ** 2 - 1) / 2


def a(x, xs):
    """
    Absorption function, equal to 1-alb(x) where alb(x) is the albedo
    :param xs: sine of the latitude of ice-sheet edge (positive)
    """
    b0 = 0.38
    a0, a2 = 0.697, 0.0779
    if abs(x) < xs:
        return a0 + a2 * (3 * x ** 2 - 1) / 2
    else:
        return b0


def second_member(x, xs):
    """
    Gives the second member of the EDP
    """
    return Q * S(x) * a(x, xs)


x0, xf = 0, 1
dx = 0.01
Nx = round((xf - x0) / dx)
x = np.linspace(x0, xf, Nx)

t0, tf = 0, 1
dt = 0.01
Nt = round((tf - t0) / dt)

xs = 0.95


def schema(x0, xf, dx, t0, tf, dt, T0, F, func):
    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)
    xs = 0.95
    I = A + B * T0(x)
    plt.plot(x, I, color=plt.get_cmap('copper')(0 / Nt))

    for n in range(Nt):
        print("________________________")
        print(n * 100 / Nt, "%")
        M = F(x)
        q = np.linalg.inv(M)
        z = delta() * I + Vec(x, xs, I, func)
        R = np.dot(q, z)
        print(R)
        # print(R)
        I = R
        plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt))


def alpha(x):
    return x * D / dx - (1 - x ** 2) * D / dx ** 2


def beta(x):
    return C / dt + 2 * (1 - x ** 2) / dx ** 2


def gamma(x):
    return - x * D / dx - (1 - x ** 2) * D / dx ** 2


def delta():
    return C / dt - 1


def Mat(x):
    A = np.zeros((Nx, Nx))
    for i in range(Nx):
        A[i, i] = beta(x[i])
    for i in range(1, Nx - 1):
        A[i, i - 1] = gamma(x[i])
        A[i, i + 1] = alpha(x[i])

    A[0, :] = 0
    A[-1, :] = 0
    A[0, 0], A[-1, -1] = delta(), delta()
    return A


def Vec(x, xs, I, func):
    vec = np.zeros(Nx)
    for i in range(Nx):
        vec[i] = func(x[i], xs)
    vec[0], vec[-1] = 0, 0
    return vec


def homogene(x, xs):
    return 0


def T0(x):
    """
    Sea level temperature as a function of latitude for the climate in 1975
    """
    return -50 * x ** 2 + 28


schema(x0, xf, dx, t0, tf, dt, T0, Mat, second_member)

plt.plot(np.linspace(x0, xf, Nx), 186.9 * np.ones(Nx), "r")
plt.xlabel(u'$x$', fontsize=26)
plt.ylabel(u'$I$', fontsize=26, rotation=0)
plt.title(u'Schéma explicite')
plt.legend()
plt.show()
