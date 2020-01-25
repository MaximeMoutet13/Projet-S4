__author__ = "maxime"
__file__ = "Numerical_schema.py"
__date__ = "24/01/20"


import matplotlib.pyplot as plt
from Parameters import *
import os

os.system("start" + "/home/maxime/Documents/S4/Projet_Modélisation_Climat/ImagesClimat")
files = os.listdir("/home/maxime/Documents/S4/Projet_Modélisation_Climat/ImagesClimat")
for i in range(0, len(files)):
    os.remove("/home/maxime/Documents/S4/Projet_Modélisation_Climat/ImagesClimat" + '/' + files[i])

xs = 0.95

t0 = 0
tf = 1
dt = 0.01
Nt = round((tf - t0) / dt)

alpha = (D * dt) / (C * dx ** 2)


# Coefficients that will be needed in the matrix
def a(x):
    return -2 * alpha * dx * x + alpha * (1 - x ** 2)


def b(x):
    return 1 + 2 * alpha * x * dx - 2 * (1 - x ** 2) * alpha - dt / C


def c(x):
    return alpha * (1 - x ** 2)


def schema(x0, xf, dx, t0, tf, dt, I0, F, xs, affichage=False):
    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)
    I = I0
    x = x0 + dx * np.arange(Nx + 1)
    if affichage:
        # plt.figure()
        # plt.xlim(x0, xf)
        plt.plot(x, I)
        # plt.show()
        # plt.savefig("/home/maxime/Documents/S4/Modélisation_Climat/ImagesClimat/im")
    for n in range(Nt - 1):
        I = F(I, x, xs)

        if affichage:
            # plt.figure()
            # plt.xlim(x0, xf)
            plt.plot(x, I)
            # plt.show()
            # plt.savefig("/home/maxime/Documents/S4/Modélisation_Climat/ImagesClimat/im" + str(n))
    plt.show()


def F(I, x, xs):
    """
    Function for the time schema
    :param I: vector I at the time n
    :param x: vector which contains Nx positions
    :return: vector I at the time n+1
    """

    A = np.eye(Nx + 1)
    for i in range(Nx + 1):
        A[i, i] = b(x[i])
    for i in range(Nx):
        A[i, i + 1] = a(x[i])
    for i in range(Nx):
        A[i + 1, i] = c(x[i])
    A[0, 0], A[0, 1], A[-1, -1], A[-1, -2] = 0, 0, 0, 0

    B = np.zeros(Nx + 1)
    for i in range(Nx + 1):
        B[i] = dt * second_member(x[i], xs) / C
    B[0] += Ix0
    B[-1] += Ix1
    return A.dot(I) + B


schema(x0, xf, dx, t0, tf, dt, It0, F, xs, affichage=True)
