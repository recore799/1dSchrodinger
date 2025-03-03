import numpy as np
import matplotlib.pyplot as plt
from .integrator import Integrator
from .shooting import Shooting
from .potentials import harmonic_oscillator

integrator = Integrator(V=harmonic_oscillator)

solver = Shooting(integrator)

# Ground state

E0, psi0 = solver.solve_state(target_nodes=0, E_min=0.8, E_max=1.2)
nodes0 = np.sum(psi0[:-1] * psi0[1:] < 0)

print("Ground state energy:", E0)
print("Ground state node count:", nodes0)

# Plot the ground state wavefunction
plt.figure(figsize=(8, 4))
plt.plot(integrator.x, psi0, label='Ground State (0 nodes)')
plt.title('Ground State Wavefunction')
plt.xlabel('x')
plt.ylabel(r'$\psi(x)$')
plt.legend()
plt.grid(True)
plt.show()
