import numpy as np

class Integrator:
    def __init__(self,V, xL=-5, xR=5, n=501):
        """
        Initializes the Numerov solver for the harmonic oscillator.

        The Schrödinger equation (in adimensional units) is:
            -ψ''(x) + x^2 ψ(x) = E ψ(x)
        which we write as:
            ψ''(x) = (x^2 - E) ψ(x).

        Parameters:
            xL: float, left boundary of the spatial grid.
            xR: float, right boundary of the spatial grid.
            n: int, number of grid points.
        """
        self.xL = xL
        self.xR = xR
        self.V = V
        self.n = n
        self.x, self.h = np.linspace(xL, xR, n, retstep=True, dtype=float)

    def numerov(self, E):
        """
        Computes the wavefunction using Numerov's method for a given energy E.

        Parameters:
            E: float, the trial energy.

        Returns:
            psi: 1D NumPy array, the computed wavefunction on the grid.

        The equation solved is:
            ψ''(x) = -[E-x^2]ψ(x)
        """
        psi = np.zeros(self.n, dtype=float)
        # Set boundary conditions.
        # We assume that at the far left (xL) we are in the classically forbidden region.
        # For simplicity, we set psi(xL)=0 and a small nonzero value at the next point.
        psi[0] = 0.0
        psi[1] = 1e-6

        # g(x)
        G = E-self.V(self.x)

        F = 1 + (self.h**2 / 12.0) * G


        # Numerov recurrence: for i = 1 to n-2:
        for i in range(1, self.n - 1):
            psi[i+1] = ((12-10*F[i])*psi[i] - F[i-1]*psi[i-1]) / F[i+1]
        return psi


