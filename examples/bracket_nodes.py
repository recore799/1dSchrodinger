import numpy as np
from tabulate import tabulate
from schrodinger import Integrator, Shooting, harmonic_oscillator

def test_bracketing():
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=1001)
    solver = Shooting(integrator)

    # Test cases: (target_nodes, expected_E)
    test_cases = [
        (0, 1.0),
        (1, 3.0),
        (2, 5.0),
        (3, 7.0),
        (4, 9.0)
    ]

    results = []
    for target_nodes, expected_E in test_cases:
        # Set energy range around expected eigenvalue
        E_min = expected_E - 1.5
        E_max = expected_E + 1.5

        try:
            bracket = solver.bracket_eigenvalue(
                target_nodes=target_nodes,
                E_min=E_min,
                E_max=E_max,
                num=1000,
                tol=1e-6
            )
            # Verify expected energy is within the bracket
            in_bracket = bracket[0] <= expected_E <= bracket[1]
            results.append([
                target_nodes,
                expected_E,
                f"[{bracket[0]:.4f}, {bracket[1]:.4f}]",
                "✅" if in_bracket else "❌"
            ])
        except ValueError as e:
            results.append([
                target_nodes,
                expected_E,
                "FAILED",
                "❌"
            ])

    # Print results
    headers = ["Target Nodes", "Expected E", "Bracket Found", "Status"]
    print(tabulate(results, headers=headers, tablefmt="github", floatfmt=".4f"))

def test_bracketing2():
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=1001)
    solver = Shooting(integrator)

    E_min = 1
    E_max = 2

    bracket = solver.bracket_eigenvalue(target_nodes=0, E_min=E_min, E_max=E_max, num=1000, tol=1e-6)

    print(bracket)


if __name__ == "__main__":
    test_bracketing2()
