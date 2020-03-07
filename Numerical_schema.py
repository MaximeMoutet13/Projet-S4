import matplotlib.pyplot as plt
from Parameters import *
import time

x0, xf = -1, 1
dx = 0.01
Nx = round((xf - x0) / dx)
x = np.linspace(x0, xf, Nx)

t0, tf = 0, 10
dt = 0.001
Nt = round((tf - t0) / dt)

xs = 0.95


def schema(x0, xf, dx, t0, tf, dt, T0, F, ice, symetrique=True):
    """
    :param x0:
    :param xf:
    :param dx:
    :param t0:
    :param tf:
    :param dt:
    :param T0: Initial T as a function of x
    :param F: Second member of the EDP as a function of x, xs
    :param ice: sine of the latitude of ice sheet edge as a function of T0
    :return:
    """

    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)

    I = A + B * T0(x)
    RI = np.zeros(Nx)
    alpha = (D * dt) / (C * dx ** 2)
    plt.plot(x, I, color=plt.get_cmap('copper')(float(0) / Nt), label="temps initial")
    xs = ice(I, x)

    for n in range(Nt):
        print(n * 100 / Nt, "%")
        for j in range(1, Nx - 1):
            if symetrique:
                RI[j] = I[j] + dt * (F(x[j], xs) - x[j] * D * (I[j + 1] - I[j - 1]) / dx + (1 - x[j] ** 2) * D * (
                        I[j + 1] + I[j - 1] - 2 * I[j]) / dx ** 2 - I[j]) / C
            else:
                RI[j] = I[j + 1] * (-2 * alpha * dx * x[j] + alpha * (1 - x[j] ** 2)) + I[j] * (1 + 2 * alpha * x[j] * dx
                                                                                            - 2 * (1 - x[
                        j] ** 2) * alpha - dt / C) + I[j - 1] * alpha * (1 - x[j] ** 2) + dt * F(x[j], xs) / C

        for j in range(1, Nx - 1):
            I[j] = RI[j]
        xs = ice(I, x)
        # plt.plot(xs, Is, "+b")

        if n % 100 == 0:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt))

        if n == Nt - 1:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt), label="temps final")
    return xs


def homogene(x, xs):
    return 0


def T0(x):
    """
    Sea level temperature as a function of latitude for the climate in 1975
    """
    return -50 * x ** 2 + 28


def latitudeS(I0, x):
    if I0[0] >= 186.9:
        return 1
    for i, I in enumerate(I0[1:]):
        if I >= 186.9 - 0.001:
            return abs((x[i - 1] + x[i]) / 2)
    return 0


debut = time.time()
print(schema(x0, xf, dx, t0, tf, dt, T0, second_member, latitudeS, symetrique=False))
fin = time.time()

plt.plot(np.linspace(x0, xf, Nx), 186.9 * np.ones(Nx), "b", label="ligne de glace")
plt.xlabel(u'$x$', fontsize=20)
plt.ylabel(u'$I (W.m^{-2})$', fontsize=20, rotation=90)
plt.title(u'Schéma explicite asymétrique')
plt.legend()
# plt.savefig("Exp_sol_asym.png")
plt.show()
