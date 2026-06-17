# HEOM-Based Electronic Friction Framework

This repository implements a Hierarchical Equations of Motion (HEOM) framework for computing electronic forces acting on classical nuclear (vibrational) degrees of freedom in vibronic systems.

The nuclei are treated classically, while the electronic subsystem is treated as a fully quantum mechanical open system. This allows the calculation of non-equilibrium electronic effects on nuclear motion in molecules coupled to metallic leads.

---

# Physical Overview

The code computes the following quantities as functions of nuclear (vibrational) coordinates:

- **Adiabatic electronic mean force**
- **Electronic friction tensor**
  - Markovian contribution
  - Non-Markovian contribution
- **Force–force correlation functions** of the stochastic electronic force
- **Electronic current (adiabatic and non-adiabatic components)**

These outputs can be used as input for classical **Langevin dynamics simulations** of nuclear motion.

---

# Model flexibility

The framework supports general vibronic models with:

- Arbitrary number of electronic sites
- Arbitrary number of vibrational degrees of freedom
- User-defined site energies and couplings
- Electron–vibration coupling strengths and functional forms
- Molecule–lead coupling strengths
- Lead parameters and spectral properties

All model parameters are defined in:
