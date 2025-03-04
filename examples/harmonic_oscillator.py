import numpy as np
import matplotlib.pyplot as plt

from schrodinger import Integrator, harmonic_oscillator

def main():
    solver = Integrator(V=harmonic_oscillator)
    energies = [2*n+1 for n in range(5)]

    plt.figure(figsize=(10,6))
    for E in energies:
        psi = solver.numerov(E)
        norm = np.sqrt(np.trapezoid(psi**2, solver.x))
        plt.plot(solver.x, psi/norm, label=f"E = {E:.1f}")

    plt.xlabel("x")
    plt.ylabel("Normalized ψ(x)")
    plt.title("Harmonic Oscillator Eigenfunctions")
    plt.legend()
    # Save the plot instead of showing it
    plt.savefig("harmonic_oscillator.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
