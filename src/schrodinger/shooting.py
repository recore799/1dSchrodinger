import numpy as np
from scipy.optimize import bisect

class Shooting:
    def __init__(self, integrator):
        """
        Initializes the shooting solver.

        Parameters:
            integrator: An instance of the Integrator class.
        """
        self.integrator = integrator

    def count_nodes(self, psi, tol):
        nodes = 0
        # Initialize with the first nonzero sign
        s_prev = np.sign(psi[2]) if abs(psi[2]) > tol else 0
        for val in psi[3:]:
            s = np.sign(val) if abs(val) > tol else 0
            # If we haven't got a previous sign (due to being near zero), update it
            if s_prev == 0 and s != 0:
                s_prev = s
                continue
            # Only count if we have a valid previous sign and it differs from the current
            if s != 0 and s != s_prev:
                nodes += 1
                s_prev = s
        return nodes - 1

    def shoot(self, E, target_nodes):
        psi = self.integrator.numerov(E)
        nodes = self.count_nodes(psi, tol=1e-6)
        diff = nodes - target_nodes
        print(f"shoot(E={E:.6f}): nodes = {nodes}, target_nodes = {target_nodes}, diff = {diff}")

        # Optionally, print a few representative psi values for debugging:
        indices = np.linspace(0, len(psi)-1, 10, dtype=int)
        psi_sample = psi[indices]
        print(f"  Sample psi values at indices {indices.tolist()}: {psi_sample}")

        return diff




    def find_eigenvalue(self, target_nodes, E_min, E_max, tol=1e-6):
        """
        Finds the eigenvalue using the shooting method.
        """
        return bisect(self.shoot, E_min, E_max, args=(target_nodes,), xtol=tol)

    def solve_state(self, target_nodes, E_min, E_max, tol=1e-6):
        """
        Solves for the eigenvalue and eigenfunction using the shooting method.
        """
        E = self.find_eigenvalue(target_nodes, E_min, E_max, tol)
        psi = self.integrator.numerov(E)
        norm = np.sqrt(np.trapezoid(psi**2, self.integrator.x))
        return E, psi / norm
