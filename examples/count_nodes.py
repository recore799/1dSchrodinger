import numpy as np
from schrodinger import Integrator, Shooting, harmonic_oscillator

def main():
    # Initialize the solver for the first excited state (E = 1.5)
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=1001)
    solver = Shooting(integrator)

    energies = [2*n+1 for n in range(6)]
    i = 0
    for E in energies:
        psi = integrator.numerov(E)
        # Normalize and count nodes
        norm = np.sqrt(np.trapezoid(psi**2, integrator.x))
        psi_norm = psi / norm
        nodes = solver.count_nodes(psi_norm, tol=1e-3)
        print(f"Expected nodes: {i} | Actual nodes: {nodes}")
        i += 1

if __name__ == "__main__":
    main()
