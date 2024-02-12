#ruff: noqa: E402 F405 F403
from math import *
from integration import *
from root import *
from differentiation import *

def q2_2():
    def fa(x: float) -> float:
        return cos(x)
    
    def fb(x: float) -> float:
        return 1/x**2
    
    def fc(x: float) -> float:
        return x**2 + x + 1
    
    def fd(x: float) -> float:
        return cos(pi/2 * x**2)

    for (i, f, a, b) in [('a', fa, 0.0, pi/2), ('b', fb, 1.0, 3.0), ('c', fc, 2.0, 4.0), ('d', fd, 0.0, 6.9)]:
        print(f"---------------\n{i}:")
        for i in range(6):
            h = 10**(-i)
            print(f"dx = {h}")
            s = int_euler(f, a, b, h)
            print(f"\tEuler     : {s}")

            s = int_trapezoid(f, a, b, h)
            print(f"\tTrapezoid : {s}")

            s = int_simpson(f, a, b, h)
            print(f"\tSimpson   : {s}")


def q2_4():
    def f(x, y):
        d2 = 5**2 - x**2 - y**2
        return 0.0 if d2 <= 0.0 else sqrt(d2)

    dx = 1e-3
    dy = 1e-3
    Nx = int(5/dx)
    Ny = int(5/dy)
    fn = [[f(x*dx, y*dy) for y in range(-Ny, Ny)] for x in range(-Nx, Nx)]
    print(euler2(fn, dx, dy))


def q2_6():
    def f_FD(mu, E):
        return 1 / (exp(40*(E - mu)) + 1)

    def f(mu):
        dE = 1e-3
        s = trapezoid([f_FD(mu, E*dE) for E in range(int(2/dE))], dE)
        return s - 1

    mu = newton(f, 0.2, dx=1e-3)
    print(mu)


def q2_7():
    dx = 1e-3
    df = df1([sin(x * dx) for x in range(int(pi/dx) + 1)], dx)
    errs = []
    for (x, y) in enumerate(df):
        errs.append((cos(x*dx) - y)**2)

    err = sum(errs)
    print(err)


def q2_8():
    dx = 1e-3
    df = df2([sin(x * dx) for x in range(int(pi/dx) + 1)], dx)
    errs = []
    for (x, y) in enumerate(df):
        errs.append((sin(x*dx) + y)**2)

    err = sum(errs)
    print(err)
