# Caltech_constellation_design

Constellation Design - a quick overview

Dr. Luigi Mascolo, Sr. Astrodynamics Eng., Orbits R&D - Nov 25, 2025

---

## Overview

This repository provides a concise introduction to key concepts used in modern constellation design.  
It includes:

- Sun-synchronous inclination computation  
- Nodal regression under J2  
- Repeating groundtrack (RGT) conditions  
- Practical search for RGT orbits  
- Groundtrack generation and visualization  
- Tools useful for quick-look analysis 

The code is intentionally lightweight so it can be explored, modified, and extended easily.

## Disclaimer

The code in this repository provides a simplified, fast, yet effective approximation for generating periodic and repeating groundtracks.

Real mission design requires a far more complex dynamical environment, including multiple perturbations and operational constraints. As a result, for some orbital configurations you may test with this code, the computed repeating groundtrack may not close perfectly, or may fail to close altogether.

This outcome is expected given the model’s intentional simplifications and should not be interpreted as a numerical or physical error.

## Contact

For questions, suggestions, or discussions: luigi.mascolo [at] planet.com


