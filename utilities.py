import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RectBivariateSpline


def evaluate(f, x, y):
    feval = f(x[:, None], y[None, :])

    return feval

def get_a(m, n, dx, dy):
    C = sp.eye(m - 1) / dx ** 2

    B_diags = [np.full(m - 1, 2), -np.ones(m - 2), -np.ones(m - 2)]
    B = sp.csr_matrix(sp.diags(B_diags, [0, -1, 1])) / dy**2 + 2*C

    A = sp.bmat([[B if i == j else -C if abs(i-j) == 1 else None for i in range(n-1)]
                 for j in range(n-1)], format='csr')

    return A


def get_b(f, dx, dy):
    _x = np.arange(dx, 1, dx)
    _y = np.arange(dy, 1, dy)

    b = f(_x[:, None], _y[None, :]).reshape(-1)

    return b


def cg(A, b, tol, maxiter=None):
    if maxiter is None:
        maxiter = len(b)

    u = np.zeros(len(b))
    r = b - A.dot(u)
    p = r
    rnorm = np.linalg.norm(r)**2

    k = 0
    while np.sqrt(rnorm) > tol and k <= maxiter:

        Ap = A.dot(p)
        alpha = rnorm / (p.T @ A.dot(p))

        u += alpha * p
        r -= alpha * Ap
        rnorm_new = np.linalg.norm(r)**2

        beta = rnorm_new / rnorm
        p = r + beta*p

        rnorm = rnorm_new
        k += 1

    return u, k


class System(object):
    def __init__(self, m, n, f, exact=None):
        self.m = m
        self.n = n
        self.f = f

        self.xmin = 0.
        self.xmax = 1.
        self.ymin = 0.
        self.ymax = 1.

        self.A = get_a(self.m, self.n, self.dx, self.dy)
        self.b = get_b(self.f, self.dx, self.dy)

        if exact is not None:
            self.exact = exact

    @property
    def dx(self):
        return (self.xmax - self.xmin) / self.n

    @property
    def dy(self):
        return (self.ymax - self.ymin) / self.m

    def u_exact(self, x=None, y=None):
        if x is None:
            x = np.linspace(self.xmin, self.xmax, self.n + 1)
            y = np.linspace(self.ymin, self.ymax, self.m + 1)

        u_exact = evaluate(self.exact, x, y)

        return u_exact

    def solve(self):
        _u = spsolve(self.A, self.b)
        _u = _u.reshape(self.m - 1, self.n - 1)
        u = np.zeros((self.m+1, self.n+1))
        u[1:-1, 1:-1] = _u

        return u

    def solve_cg(self, tol, maxiter=None):
        _u, niter = cg(self.A, self.b, tol, maxiter)
        _u = _u.reshape(self.m - 1, self.n - 1)
        u = np.zeros((self.m + 1, self.n + 1))
        u[1:-1, 1:-1] = _u

        return u, niter

    def error(self, k=3):
        u = self.solve()
        x = np.linspace(self.xmin, self.xmax, self.n + 1)
        y = np.linspace(self.ymin, self.ymax, self.m + 1)

        ip = RectBivariateSpline(x, y, u)

        xi = np.linspace(self.xmin, self.xmax, (self.n+1) * k)
        yi = np.linspace(self.ymin, self.ymax, (self.m+1) * k)

        u_exact = self.u_exact(xi, yi)
        u_ip = ip(xi, yi)

        e = np.sum(abs(u_exact - u_ip) * self.dx * self.dy)

        return e
