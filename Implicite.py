from parametres import*


x0, xf = -1, 0
dx = 0.01
Nx = round((xf - x0) / dx)
x = np.linspace(x0, xf, Nx)
t0, tf = 0, 10
dt = 0.01
Nt = round((tf - t0) / dt)


def schema(x0, xf, dx, t0, tf, dt, T0, F, func, q):
    Nx = round((xf - x0) / dx)
    Nt = round((tf - t0) / dt)

    x = np.linspace(x0, xf, Nx)
    I = A_cste + B_cste * T0(x)
    xs = latitudeS(I, x)

    M = F(x)

    plt.plot(x, I, color=plt.get_cmap('copper')(0 / Nt), label="$t_0=0$ ")
    plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(0 / Nt))
    for n in range(Nt):
        # print("________________________")
        print("Ploting...", n * 100 / Nt, "%")
        # q = np.linalg.inv(M)
        z = delta() * I + Vec(x, xs, func, q)
        R = np.linalg.solve(M, z)
        I = R
        I[0] = I[1]
        I[-1] = I[-2]
        xs = latitudeS(I, x)
        # print("/////")
        # print(xs)
        # print("/////")
        if n % 100 == 0:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt))
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))
        if n == Nt - 1:
            plt.plot(x, I, color=plt.get_cmap('copper')(float(n) / Nt), label="$t_f=10$")
            plt.plot(x + 1, I[::-1], color=plt.get_cmap('copper')(float(n) / Nt))
    return(xs)

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

    A[0, :], A[-1, :] = 0, 0
    A[0, 0], A[-1, -1] = beta(x[0]), delta() + 1

    return A

def latitudeS(I0, x):
    if I0[0] >= 186.9:
        return 1
    for i, I in enumerate(I0[1:]):
        if I >= 186.9 - 0.001:
            return abs((x[i - 1] + x[i]) / 2)
    return 0

def Vec(x, xs, func, q):
    vec = np.zeros(Nx)
    for i in range(Nx):
        vec[i] = func(x[i], xs, q)
    return vec

print(schema(x0, xf, dx, t0, tf, dt, T0, Mat, second_member, Q_cste))

plt.plot(np.linspace(x0, xf + 1, Nx), 186.9 * np.ones(Nx), "b", label="ligne de glace")
plt.xlabel(u'$x$', fontsize=20)
plt.ylabel(u'$I (W.m^{-2})$', fontsize=20, rotation=90)
plt.title(u'Euler implicite avec second membre')
plt.legend()
plt.savefig("im99.png")
plt.show()
