import numpy as np
import matplotlib.pyplot as plt
from schrodinger import Integrator, harmonic_oscillator

def main():
    # Configure style
    plt.style.use("seaborn-v0_8-darkgrid")
    plt.rcParams.update({"font.size": 12, "figure.autolayout": True})

    solver = Integrator(V=harmonic_oscillator)
    energies = [2*n+1 for n in range(5)]  # n=0 (ground state) to n=4

    # Initialize figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Color palette for wavefunctions
    colors = plt.cm.viridis(np.linspace(0, 1, len(energies)))

    # Plot eigenfunctions
    for i, E in enumerate(energies):
        psi = solver.numerov(E)
        norm = np.sqrt(np.trapezoid(psi**2, solver.x))
        ax.plot(solver.x, psi/norm + E,  # Offset vertically by energy
                lw=2.5, color=colors[i],
                label=fr"$n = {i}$; $E_{i} = {E:.1f}$")

        # Fill under curve with transparency
        ax.fill_between(solver.x, E, psi/norm + E,
                        color=colors[i], alpha=0.15)

    # Formatting
    ax.set_xlabel("Position ($x$)", fontsize=14)
    ax.set_ylabel("Energy → Wavefunction Amplitude", fontsize=14)
    ax.set_title("Harmonic Oscillator Eigenfunctions", fontsize=16, pad=20)

    # Add energy levels
    for E in energies:
        ax.axhline(E, color='gray', linestyle='--', lw=1, alpha=0.7)

    ax.legend(loc='upper right', frameon=True, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)

    # Save figure for README
    plt.savefig("examples/harmonic_oscillator_eigenfunctions.png",
                dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    main()
