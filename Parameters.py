__author__ = "maxime"
__file__ = "Parameters.py"
__date__ = "24/01/20"

import numpy as np

# I = A + B T
A = 201.4  # W.m-2
B = 1.45  # W.m-2.°C-1

########################################################################################################################
"""
Left member of the EDP
"""

Rt = 6371e3  # Earth radius in meters
C = 14 / B  # Heat capacity of the relevant layers of athmosphere divided by B
D = 3.5e6 * B * Rt ** 2  # Phenomenological thermal diffusion coefficient which absorbed other things

########################################################################################################################

""""
Right Member of the EDP
"""

# Solar constant divided by 4
Q = 340  # W.m-2

S2 = -0.482


def S(x):
    """
    Mean annual distribution of radiation at each latitude
    """
    return 1 + S2 * (3 * x ** 2 - 1) / 2


b0 = 0.38
a0, a2 = 0.697, 0.0779


def a(x, xs):
    """
    Absorption function, equal to 1-alb(x) where alb(x) is the albedo
    :param xs: sine of the latitude of ice-sheet edge
    """
    if x > xs:
        return b0
    else:
        return a0 + a2 * (3 * x ** 2 - 1) / 2


def second_member(x, xs):
    """
    Gives the second member of the EDP
    """
    return Q * S(x) * a(x, xs)


########################################################################################################################

"""
Conditions
"""

x0 = -1
xf = 1
dx = 0.01
Nx = round((xf - x0) / dx)

X = x0 + dx * np.arange(Nx + 1)

# Initial conditions


def T0(x):
    """
    Sea level temperature as a function of latitude for the climate in 1975
    """
    return -50 * x ** 2 + 28


# Initial condition
It0 = A + B * T0(X)

# Boundary conditions
Ix0, Ix1 = -44.2, -44.2
