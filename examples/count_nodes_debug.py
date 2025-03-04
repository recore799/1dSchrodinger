import numpy as np
from tabulate import tabulate
from schrodinger import Integrator, ShootingDebug, harmonic_oscillator

def main():
    integrator = Integrator(V=harmonic_oscillator, xL=-5, xR=5, n=1001)
    solver = ShootingDebug(integrator)

    # Known energies and expected nodes
    test_cases = [
        {"energy": 1.0, "expected_nodes": 0},
        {"energy": 3.0, "expected_nodes": 1},
        {"energy": 5.0, "expected_nodes": 2},
        {"energy": 7.0, "expected_nodes": 3},
        {"energy": 9.0, "expected_nodes": 4},
        {"energy": 11.0, "expected_nodes": 5},
    ]

    results = []
    for case in test_cases:
        psi = integrator.numerov(case["energy"])
        norm = np.sqrt(np.trapezoid(psi**2, integrator.x))
        psi_norm = psi / norm
        nodes, node_locations = solver.count_nodes(psi_norm, tol=1e-3)

        # Format node locations as a string
        node_locations_str = ", ".join(f"{x:.4f}" for x in node_locations)

        results.append([
            case["energy"],
            case["expected_nodes"],
            nodes,
            node_locations_str,
            "✅" if nodes == case["expected_nodes"] else "❌"
        ])

    # Print table
    headers = ["Energy", "Expected Nodes", "Actual Nodes", "Node Locations", "Status"]
    print(tabulate(results, headers=headers, tablefmt="github"))

if __name__ == "__main__":
    main()
