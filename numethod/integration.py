from typing import Callable, List

def int_euler(f: Callable[[float], float], a: float, b: float, dx: float =1e-6) -> float:
    """
    A simple rectangular slice approximation of an integral.
    
    f: Function to integrate
    a, b: limits of integration
    dx: step size. Defaults to 1e-6
    """

    return dx * sum(f(a + dx*j) for j in range(int((b - a) / dx) + 1))


def euler(f: List[float], dx: float) -> float:
    """
    A simple rectangular slice approximation of an integral.
    
    f: List of values each separated by interval of dx
    dx: step size
    
    Limits of integration are `f[0]` and `f[0] + dx * len(f)`
    """

    return dx * sum(f[0:-1])


def int_trapezoid(f: Callable[[float], float], a: float, b: float, dx: float =1e-6) -> float:
    """
    A simple trapezoid slice approximation of an integral.
    
    f: Function to integrate
    a, b: limits of integration
    dx: step size. Defaults to 1e-6
    """

    return int_euler(f, a, b, dx) - (f(a) + f(b)) / 2 *dx


def trapezoid(f: List[float], dx: float) -> float:
    """
    A simple trapezoid slice approximation of an integral.
    
    f: List of values each separated by interval of dx
    dx: step size
    
    Limits of integration are `f[0]` and `f[0] + dx * len(f)`
    """

    return euler(f, dx) - (f[0] + f[-1]) / 2 * dx


def int_simpson(f: Callable[[float], float], a: float, b: float, dx: float =1e-6) -> float:
    """
    Composite Simpson's 1/3 rule for calculating integrals
    
    f: Function to integrate
    a, b: limits of integration
    dx: step size. Defaults to 1e-6
    """

    s = f(a) + f(b)
    n = (b - a) / dx

    for i in range(1, int(n/2) + 1):
        s += 4 * f(a + (2*i - 1)*dx)
    
    for i in range(1, int(n/2)):
        s += 2 * f(a + 2*i*dx)
    
    if int(n) % 2 == 1:
        s += 0.25 * (5*f(b) + 8*f(b - dx) - f(b - 2*dx))
    
    return s * 1/3 *dx


def simpson(f: List[float], dx: float) -> float:
    """
    Composite Simpson's 1/3 rule for calculating integrals
    
    f: List of values each separated by interval of dx
    dx: step size
    
    Limits of integration are `f[0]` and `f[0] + dx * len(f)`
    """

    s = f[0] + f[-1]
    n = len(f)

    for i in range(1, n-1, 2):
        s += f[i - 1] + 4 * f[i] + f[i + 1]

    if n % 2 == 1:
        s += 0.25 * (5*f[-1] + 8*f[-2] - f[-3])
    
    return s * 1/3 * dx


def euler2(f: List[List[float]], dx: float, dy: float) -> float:
    """
    A simple cuboid slice approximation of a double integral.
    
    f[x][y]: List of values each separated by interval of (dx, dy)
    dx: step size in x direction
    dy: step size in y direction
    
    Limits of integration are `f[0]` and `(f[0][0] + dx * len(f), f[0][1] + dy * len(f))`
    """

    return dx * dy * sum(sum(f[i][0:-1]) for i in range(len(f)))
