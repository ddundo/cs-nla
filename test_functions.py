import numpy as np
from utilities import System
import matplotlib.pyplot as plt


def f1(x, y):
    return 13 * np.pi**2 * np.sin(2*np.pi*x) * np.sin(3*np.pi*y)


def f2(x, y):
    return -(x-1)**3 * (42*x**2 - 24*x + 2) * y*(y-1) - 2*x**2 * (x-1)**5


def sol1(x, y):
    return np.sin(2*np.pi*x) * np.sin(3*np.pi*y)


def sol2(x, y):
    return (x-1)**5 * x**2 * y*(y-1)


m = 100
n = 100

system = System(m, n, f2, sol2)

u, niter = system.solve_cg(1.e-8)

print(niter)
print(system.error())
