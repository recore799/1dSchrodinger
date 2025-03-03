import numpy as np
import matplotlib.pyplot as plt
from schrodinger import Integrator, harmonic_oscillator

def test_ground_state():
    solver = Integrator(V=harmonic_oscillator)
    E = 1
    psi = solver.numerov(E)
    norm = np.sqrt(np.trapezoid(psi**2, solver.x))
    psi /= norm

    # Debug plot
    plt.plot(solver.x, psi, label=f"E = {E:.1f}")
    plt.xlabel("x")
    plt.ylabel("ψ(x)")
    plt.title("Wavefunction (Debug)")
    plt.legend()
    plt.grid(True)
    plt.show()

    assert np.abs(psi[-1]) < 1e-2 # Test boundary condition



def test_node_count():
    solver = Integrator(V=harmonic_oscillator)
    psi = solver.numerov(1.5)
    nodes = np.sum(psi[:-1] * psi[1:] < 0)
    assert nodes == 1
