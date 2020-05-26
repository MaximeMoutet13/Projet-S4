from param import*


x0, xf = -1, 0
dx = 0.01
Nx = round((xf - x0) / dx)
x = np.linspace(x0, xf, Nx)

t0, tf = 0, 10
dt = 0.01
Nt = round((tf - t0) / dt)

xs = 0.95


def schema(x0, xf, dx, t0, tf, dt, T0, F, func):
    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)
    xs = 0.95
    I = A_cste + B_cste * T0(x)
    plt.plot(x, I, color=plt.get_cmap('copper')(0 / Nt), label="temps initial")
    plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(0 / Nt))
    for n in range(Nt):
        # print("________________________")
        print(n * 100 / Nt, "%")
        M = F(x)
        # q = np.linalg.inv(M)
        z = delta() * I + Vec(x, xs, I, func)
        R = np.linalg.solve(M, z)
        # print(R)
        # print(R)
        I = R
        if n % 100 == 0:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt))
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))
        if n == Nt - 1:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt), label="temps final")
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))
def alpha(x):
    return ((x * D) / dx) - (((1 - x ** 2) * D) / (dx ** 2))


def beta(x):
    return (tau / dt) + 2 * (((1 - x ** 2) * D) / (dx ** 2))  + 1


def gamma(x):
    return - (x * D / dx) -((1 - x ** 2) * D) / (dx ** 2)


def delta():
    return tau / dt


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


schema(x0, xf, dx, t0, tf, dt, T0, Mat, second_member)

plt.plot(np.linspace(x0, xf + 1, Nx), 186.9 * np.ones(Nx), "b", label="ligne de glace")
plt.xlabel(u'$x$', fontsize=20)
plt.ylabel(u'$I (X.m^{-2})$', fontsize=20, rotation=90)
plt.title(u'Schéma implicite')
plt.legend()
# plt.savefig("Imp.png")
plt.show()
