from typing import List, Callable

def diff_1(f: Callable[[float], float], x: float, dx: float =1e-6) -> float:
    """
    Calculates first derivative using 2 point formula
    """
    return (f(x + dx) - f(x-dx)) / 2 / dx


def diff_2(f: Callable[[float], float], x: float, dx: float =1e-6) -> float:
    """
    Calculates second derivative using 3 point formula
    """
    return (f(x+dx) - 2*f(x) + f(x-dx)) / dx**2


def diff_1_4(f: Callable[[float], float], x: float, dx: float =1e-6) -> float:
    """
    Calculates first derivative using 4 point formula
    """
    return (f(x - 2*dx) - 8*f(x-dx) + 8*f(x+dx) - f(x+2*dx))/12/dx


def diff_2_5(f: Callable[[float], float], x: float, dx: float =1e-6) -> float:
    """
    Calculates second derivative using 5 point formula
    """
    return (-f(x - 2*dx) + 16*f(x-dx) - 30*f(x) + 16*f(x+dx) - f(x+2*dx))/12/dx**2


def df1(f: List[float], dx: float) -> List[float]:
    """
    Calculates first derivative using 2 point formula
    """
    d = [0.0]
    for x in range(1, len(f) - 1):
        d.append((f[x + 1] - f[x - 1]) / 2 / dx)
    d.append(d[-1])
    d[0] = d[1]
    return d


def df2(f: List[float], dx: float) -> List[float]:
    """
    Calculates second derivative using 3 point formula
    """
    d = [0.0]
    for x in range(1, len(f) - 1):
        d.append((f[x+1] - 2*f[x] + f[x-1]) / dx**2)
    d.append(d[-1])
    d[0] = d[1]
    return d


def df1_4(f: List[float], dx: float) -> List[float]:
    """
    Calculates first derivative using 4 point formula
    """
    d = [0.0]
    for x in range(1, len(f) - 1):
        d.append((f[x-2] - 8*f[x-1] + 8*f[x+1] - f[x+2])/12/dx)
    d.append(d[-1])
    d[0] = d[1]
    return d


def df2_5(f: List[float], dx: float) -> List[float]:
    """
    Calculates second derivative using 5 point formula
    """
    d = [0.0]
    for x in range(1, len(f) - 1):
        d.append((-f[x-2] + 16*f[x-1] - 30*f[x] + 16*f[x+1] - f[x+2])/12/dx**2)
    d.append(d[-1])
    d[0] = d[1]
    return d

