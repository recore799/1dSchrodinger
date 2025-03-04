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

    def boundary(self, E):
        psi = self.integrator.numerov(E)
        # Return the value of the wavefuncion at the boundary
        return psi[-1]

    def bracket_eigenvalue(self, target_nodes, E_min, E_max, num=1000, tol=1e-6):
        """
        Scans energies between E_min and E_max to find an interval where
        the number of nodes changes from target_nodes to something else.
        """
        energies = np.linspace(E_min, E_max, num)
        bracket = None
        prev_nodes = None
        for E in energies:
            psi = self.integrator.numerov(E)
            nodes = self.count_nodes(psi, tol)
            if prev_nodes is None:
                prev_nodes = nodes
                continue
            # If we've reached the target (or just crossed into it) and then left it,
            # we assume the eigenvalue lies in this small interval.
            if prev_nodes == target_nodes and nodes != target_nodes:
                bracket = (E_prev, E)
                break
            E_prev = E
            prev_nodes = nodes
        if bracket is None:
            raise ValueError(f"Could not bracket eigenvalue for target_nodes = {target_nodes}")
        print(f"Bracketing interval found: {bracket}")
        return bracket

    def eigenvalue(self, target_nodes, E_min, E_max, tol=1e-6):
        """
        Finds the eigenvalue for a given target number of nodes.
        First, it brackets the eigenvalue using the node count, using the boundary condition (psi(x_R)=0)."""

        # Bracket the eigenvalue using node counting
        E_bracket = self.bracket_eigenvalue(target_nodes, E_min, E_max, tol=tol)
        # Refine using boundary function
        E_eigen = bisect(self.boundary, E_bracket[0], E_bracket[1], xtol=tol)
        return E_eigen

        
    def solve_state(self, target_nodes, E_min, E_max, tol=1e-6):
        """
        Solves for the eigenvalue and eigenfunction using the shooting method.
        """
        E = self.eigenvalue(target_nodes, E_min, E_max, tol)
        psi = self.integrator.numerov(E)
        norm = np.sqrt(np.trapezoid(psi**2, self.integrator.x))
        return E, psi / norm
