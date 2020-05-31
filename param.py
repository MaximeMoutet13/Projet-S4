import numpy as np
import matplotlib.pyplot as plt

# I = A + B T
A_cste = 201.4
B_cste = 1.45

D = 0.310
tau = 14 / 1.45
Q_cste = 340

def S(x):
    """
    Mean annual distribution of radiation at each latitude
    """
    S2 = -0.482
    return 1 + S2 * (3 * x ** 2 - 1) / 2


def a(x, xs):
    """
    Absorption function, equal to 1-alb(x) where alb(x) is the albedo
    :param xs: sine of the latitude of ice-sheet edge
    """
    a1 = 0.38
    a0 = 0.69
    if abs(x) < xs:
        return a0
    else:
        return a1


def second_member(x, xs, q):
    """
    Gives the second member of the EDP
    """
    return q * S(x) * a(x, xs)


def homogene(x, xs, q):
    return 0


def T0(x):
    """
    Sea level temperature as a function of latitude for the climate in 1975
    """
    return -50 * x ** 2 + 28

Ix0, Ix1 = -44.2, -44.2
