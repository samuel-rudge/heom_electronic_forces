# HEOM-Based Electronic Friction Framework

This repository implements a Hierarchical Equations of Motion (HEOM) framework for computing non-equilibrium electronic forces, currents, and friction kernels in vibronic quantum transport systems.

The nuclei are treated as classical degrees of freedom, while the electronic subsystem is treated as a fully quantum mechanical open system coupled to fermionic leads. The framework enables construction of Langevin-type nuclear dynamics with electronic backreaction derived from first principles.

---

# Physical Model

The code computes electronic feedback on nuclear motion in vibronic systems consisting of:

- Classical vibrational coordinate(s) x
- Quantum electronic subsystem (few-site system)
- Fermionic leads at finite bias and temperature

It evaluates:

- Adiabatic electronic mean force
- Electronic current (adiabatic and non-adiabatic contributions)
- Electronic friction tensor (Markovian approximation)
- Force–force correlation functions
- Optional non-Markovian friction integrands

These quantities can be used to construct generalized Langevin equations for nuclear dynamics.

---

# Model Flexibility

The framework supports general vibronic transport models with:

- Arbitrary number of electronic sites
- Classical and optional quantum vibrational modes
- Linear electron–vibration coupling
- User-defined electronic structure and interactions
- Fermionic leads with finite bias, finite temperature, and Lorentzian or wide-band spectral densities
- HEOM hierarchy truncation control

All parameters are defined in input_parameters.py and system.py. The first contains numerical and physical parameters, while the second defines the Hamiltonian structure and coupling operators.

---

# Quick Start

The program is executed from the terminal in the main directory containing all files.

For steady-state observables, run: python3 friction_heom_main.py ss

This computes steady-state quantities at the nuclear coordinate grid, including adiabatic forces and electronic currents. The outputs are written as .dat files in the same directory, including adiabatic_force.dat, adiabatic_force_mol.dat, adiabatic_force_molleads.dat, current_ad.dat, and current_na.dat. Additional intermediate files such as .p and .npy files may also be generated. A file simulation_info.dat is written containing information about runtime, grid size, and memory usage.

For Markovian friction and force–force correlations, run: python3 friction_heom_main.py markovian

This computes the Markovian electronic friction tensor and the force–force correlation function. The outputs include friction.dat, friction_mol.dat, corrfunc.dat, and corrfunc_mol.dat.

If print_integrand_yn is set to True in input_parameters.py, then additional non-Markovian quantities are computed at a single value of the nuclear coordinate. These are written to friction_integrand_heom.dat and corrfunc_integrand_heom.dat.

---

# Example Model (Holstein Junction)

A typical setup corresponds to a single-site Holstein model coupled to two electronic leads at finite bias.

The molecular Hamiltonian is of the form:

H_mol = (epsilon + lambda * x) d^\dag d + Omega/2 (x^2 + p^2)

with typical parameter values:

epsilon = 50 meV  
lambda = 10 meV  
Omega = 30 meV  
Gamma = 50 meV  
bias voltage = 0.1 V  
temperature = 300 K  
lead bandwidth = 10 eV (Lorentzian spectral density)

All energies are expressed in electron volts, and the model uses dimensionless units for the vibrational coordinate.

---

# Units

Energies are given in eV. Temperature is given in Kelvin and converted internally. The vibrational coordinate is dimensionless. The code assumes natural units with ħ = 1 internally.

---

# Numerical Method

The implementation is based on the Hierarchical Equations of Motion (HEOM) formalism with:

- Padé or barycentric decomposition of Fermi functions
- Krylov subspace GMRES solvers
- Adaptive time stepping using Dormand–Prince schemes
- Hierarchy truncation controlled by Nmax
- Lorentzian or wide-band lead spectral functions

---

# Output Files

Steady-state mode (ss):

- adiabatic_force.dat: total adiabatic mean force
- adiabatic_force_mol.dat: molecular contribution to force
- adiabatic_force_molleads.dat: molecule–lead coupling contribution
- current_ad.dat: adiabatic current
- current_na.dat: non-adiabatic current
- simulation_info.dat: simulation metadata including grid size and memory estimates

Markovian mode (markovian):

- friction.dat: Markovian electronic friction tensor
- friction_mol.dat: molecular contribution to friction
- corrfunc.dat: force–force correlation function
- corrfunc_mol.dat: molecular contribution to correlation function

Non-Markovian diagnostics (optional when print_integrand_yn is True):

- friction_integrand_heom.dat
- corrfunc_integrand_heom.dat

---

# System Definition

All model definitions and numerical parameters are controlled in input_parameters.py, while system.py defines the Hamiltonian structure, coupling terms, and vibrational form.

---

# Parallelization Strategy

The steady-state mode is optimized using efficient linear algebra routines (including GMRES and MKL acceleration). The Markovian mode parallelizes over the nuclear coordinate grid to improve performance.

---

# Notes

All output files are written to the working directory and will be overwritten upon rerunning the code. The Markovian mode requires a defined grid of nuclear coordinates. For single-point non-Markovian evaluation, the nuclear grid must be reduced to a single point by setting x_max_total = x_min_total + dx_grid.
