from typing import Callable
from differentiation import diff_1

def bisection(f: Callable[[float], float], a: float, b: float, tol: float =1e-6) -> float:
    """
    Uses the bisection method to find a value x between a and b for which 
    f(x) = 0, to within the tolerance given.
    
    Raises error if both a and b are on the same side of the root
    """

    if f(a) * f(b) > 0:
        raise

    dx = abs(a - b)
    x = (a + b) / 2
    while dx > tol:
        if f(x) * f(a) > 0:
            a = x
        else:
            b = x
        dx = abs(a - b)
        x = (a + b) / 2

    return x


def newton(f: Callable[[float], float], guess: float, tol: float =1e-6, dx: float =1e-6) -> float:
    """
    Uses Newton's method to find a value x near "guess" for which 
    f(x) = 0, to within the tolerance given. The derivative is calculated by the slope of line with end
    points as x+dx and x-dx
    """

    h = dx
    x = guess
    dx = 2*tol
    while abs(dx) > tol:
        df = diff_1(f, x, h)
        dx = f(x)/df
        x = x - dx
    return x
