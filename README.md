# 1D Schrödinger Equation Solver

A Python package for solving the 1D Schrödinger equation numerically using the **Numerov method**. This package provides tools to compute wavefunctions and energy eigenvalues for arbitrary potentials, with built-in support for the harmonic oscillator and hydrogen atom potentials.

## Features

- **Numerov Integration**: Solve the 1D Schrödinger equation for any potential.
- **Harmonic Oscillator**: Predefined potential for the quantum harmonic oscillator.
- **Hydrogen Atom**: Predefined potential for the radial Schrödinger equation of the hydrogen atom.
- **Wavefunction Visualization**: Plot wavefunctions and probability densities.
- **Extensible**: Easily add new potentials and solvers.

## Installation

### Prerequisites

- Python 3.8 or higher
- `numpy`, `scipy`, and `matplotlib` (installed automatically if using `pip`)

### Install from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/recore799/1dSchrodinger.git
   cd 1dSchrodinger
   ```
   
2. Install package in development mode:
   ```bash
   pip install -e .
   ```
   
