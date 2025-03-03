import numpy as np
from schrodinger import Integrator, Shooting, harmonic_oscillator

def test_shooting_ground_state():
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=10, n=501)
    solver = Shooting(integrator)
    E, psi = solver.solve_state(target_nodes=0, E_min=0, E_max=2)
    assert np.isclose(E, 0.5, atol=1e-4) # Ground state energy
    assert np.sum(psi[:-1] * psi[1:] < 0) == 0

def test_shooting_first_excited_state():
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=501)
    solver = Shooting(integrator)
    E, psi = solver.solve_state(target_nodes=1, E_min=1, E_max=3)
    assert np.isclose(E, 1.5, atol=1e-4) #First excited state energy
    assert np.sum(psi[:-1] * psi[1:] <0) == 1 # 1 node
