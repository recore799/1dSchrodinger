import numpy as np
import matplotlib.pyplot as plt
from schrodinger import Integrator, Shooting, harmonic_oscillator

def main():
    # Initialize the integrator and shooting solver
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=501)
    solver = Shooting(integrator)

    # Find the ground state (0 nodes)
    E0, psi0 = solver.solve_state(target_nodes=1, E_min=2, E_max=4)
    nodes0 = solver.count_nodes(psi0, tol=1e-6)

    print(f"Ground state energy: {E0:.6f}")
    print(f"Ground state node count: {nodes0}")

    # Plot the ground state wavefunction
    plt.figure(figsize=(8, 4))
    plt.plot(integrator.x, psi0, label='Ground State (0 nodes)')
    plt.title('Ground State Wavefunction')
    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
