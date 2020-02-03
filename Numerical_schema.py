import numpy as np
import matplotlib.pyplot as plt
from Parameters import second_member, D, C, A, B, T0
import time

x0, xf = -1, 1
dx = 0.01

t0, tf = 0, 10
dt = 0.001

xs = 0.95

alpha = (D * dt) / (C * dx ** 2)


def schema(x0, xf, dx, t0, tf, dt, F):
    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)

    I = A + B * T0(x)
    RI = np.zeros(Nx)

    for n in range(Nt):
        for j in range(1, Nx - 1):
            RI[j] = I[j + 1] * (-2 * alpha * dx * x[j] + alpha * (1 - x[j] ** 2)) + I[j] * (1 + 2 * alpha * x[j] * dx
                                                                                            - 2 * (1 - x[
                        j] ** 2) * alpha - dt / C) + I[j - 1] * alpha * (1 - x[j] ** 2) + dt * F(x[j], xs) / C
        for j in range(1, Nx - 1):
            I[j] = RI[j]

        if n % 100 == 0:
            plotlabel = "t = %1.2f" % (n * dt)
            plt.plot(x, I, label=plotlabel, color=plt.get_cmap('copper')(float(n) / Nt))
    plt.show()


def homogene(x, xs):
    return 0


debut = time.time()
schema(x0, xf, dx, t0, tf, dt, homogene)
fin = time.time()
print(fin - debut)
