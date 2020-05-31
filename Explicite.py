from parametres import *
import time
import os


x0, xf = -1, 0
dx = 0.01
Nx = round((xf - x0) / dx)
x = np.linspace(x0, xf, Nx)

t0, tf = 0, 10
dt = 0.0001
Nt = round((tf - t0) / dt)

xs = 0.95


def schema(x0, xf, dx, t0, tf, dt, T0, F, ice):


    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)

    I = A_cste + B_cste * T0(x)
    RI = np.zeros(Nx)
    plt.plot(x, I, color=plt.get_cmap('copper')(float(0) / Nt), label="$t_0=0$")
    plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(0 / Nt))
    xs = ice(I, x)
    for n in range(Nt):
        print("Ploting...", int(n * 100 / Nt), "%")
        for j in range(Nx):
            if j==0:
                RI[j] = I[j] + dt * (
                            F(x[j + 1], xs, Q_cste) - x[j + 1] * D * (I[j + 2] - I[j]) / dx + (1 - x[j + 1] ** 2) * D * (
                            I[j + 2] + I[j] - 2 * I[j + 1]) / dx ** 2 - I[j + 1]) / tau
            elif j == Nx - 1:
                RI[j] = I[j - 1] + dt * (
                        F(x[j - 1], xs, Q_cste) - x[j - 1] * D * (I[j] - I[j - 2]) / dx + (
                            1 - x[j - 1] ** 2) * D * (
                                I[j] + I[j - 2] - 2 * I[j - 1]) / dx ** 2 - I[j]) / tau
            else:
                RI[j] = I[j] + dt * (
                            F(x[j], xs, Q_cste) - x[j] * D * (I[j + 1] - I[j - 1]) / dx + (1 - x[j] ** 2) * D * (
                            I[j + 1] + I[j - 1] - 2 * I[j]) / dx ** 2 - I[j]) / tau
        I = RI
        xs = ice(I, x)

        if n % 100 == 0:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt))
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))

        if n == Nt - 1:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt), label="$t_f=10$")
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))
    return xs


def latitudeS(I0, x):
    if I0[0] >= 186.9:
        return 1
    for i, I in enumerate(I0[1:]):
        if I >= 186.9 - 0.001:
            return abs((x[i - 1] + x[i]) / 2)
    return 0


v = schema(x0, xf, dx, t0, tf, dt, T0, second_member, latitudeS)
print(v)

plt.plot(np.linspace(x0, xf + 1, Nx), 186.9 * np.ones(Nx), "b", label="ligne de glace")
plt.xlabel(u'$x$', fontsize=20)
plt.ylabel(u'$I (W.m^{-2})$', fontsize=20, rotation=90)
plt.title(u'Schéma explicite avec second membre')
plt.legend()
# plt.savefig("im35.png")
plt.show()
