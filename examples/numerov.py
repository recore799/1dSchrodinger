import numpy as np
import matplotlib.pyplot as plt
from schrodinger import Integrator, Shooting, harmonic_oscillator

# Assume 'solver' is already created with the harmonic oscillator potential


solver = Integrator(V=harmonic_oscillator)
# Define your robust node counting method as a function (or use the method you implemented)
def count_nodes(psi, tol=1e-6):
    nodes = 0
    # Find the first significant value for the initial sign.
    s_prev = np.sign(psi[0]) if abs(psi[0]) > tol else 0
    for val in psi[1:]:
        s = np.sign(val) if abs(val) > tol else 0
        if s_prev == 0 and s != 0:
            s_prev = s
            continue
        if s != 0 and s != s_prev:
            nodes += 1
            s_prev = s
    return nodes

# Choose a range of energies to test (for the first excited state you might be interested in energies between 2 and 4)
energies = np.linspace(2.0, 4.0, 21)
node_counts = []

print("Energy scan for node counts:")
for E in energies:
    psi = solver.numerov(E)
    nodes = count_nodes(psi, tol=1e-6)
    node_counts.append(nodes)
    print(f"Energy: {E:.6f} -> Nodes: {nodes}")



energies1 = [2*n+1 for n in range(5)]

plt.figure(figsize=(10,6))
for E in energies:
    psi = solver.numerov(E)
    norm = np.sqrt(np.trapezoid(psi**2, solver.x))
    plt.plot(solver.x, psi/norm, label=f"E = {E:.1f}")

plt.xlabel("x")
plt.ylabel("Normalized ψ(x)")
plt.title("Harmonic Oscillator Eigenfunctions")
plt.legend()
plt.show()



# Optional: Plot the node count vs energy to visualize the behavior.
plt.figure(figsize=(8, 4))
plt.plot(energies, node_counts, 'o-', label='Node Count')
plt.xlabel("Energy")
plt.ylabel("Number of Nodes")
plt.title("Node Count vs Energy")
plt.grid(True)
plt.legend()
plt.show()
