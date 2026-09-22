---
date: "2023-05-01T00:00:00+08:00"
title: "Aerospace Engineering Education — AeroFuse"
draft: false
index: true
weight: 2
description: A computational framework to teach aerospace engineering with multidisciplinary design optimization.
image: aerostruct-wing-tail.svg
categories:
    - Aerospace
    - Education
---

<!-- **A computational framework to teach aerospace engineering with multidisciplinary design optimization.** -->

This work is an interactive, computational platform for aircraft design designed to be used as an instructional tool: [AeroFuse](https://github.com/GodotMisogi/AeroFuse.jl).
* Used at The Hong Kong University of Science and Technology, and Imperial College London.
* Implemented geometry, aerodynamics (vortex-lattice), structures (beam-element), propulsion (blade-element momentum theory), and flight dynamics (rigid-body integrators) for coupled analyses.
* Enabled adjoint-based solutions of inverse and optimization problems between disciplines with automatic differentiation.

![eXtended Design Structure Matrix of features in AeroFuse](AeroFuse.svg)

AeroFuse couples the disciplinary solvers into a single differentiable model, so complete aircraft configurations can be analyzed and optimized from one interface for conceptual and preliminary design.

<!-- ![Spanwise lift distribution: initial vs. optimized wing](spanwise-cl.svg) -->

![Coupled aerostructural optimization of an aircraft with a T-tail configuration](aerostruct-wing-tail.svg)

The platform is open source and documented for classroom use:

* Repository: [github.com/GodotMisogi/AeroFuse.jl](https://github.com/GodotMisogi/AeroFuse.jl)
* Documentation: [godotmisogi.github.io/AeroFuse.jl/stable](https://godotmisogi.github.io/AeroFuse.jl/stable/)

### Publications

* Arjit Seth, Stephane Redonnet, and Rhea P. Liem. "MADE: A Multidisciplinary Computational Framework for Aerospace Engineering Education". *IEEE Transactions on Education* 66.6 (2023), pp. 622–631. [doi:10.1109/TE.2023.3281825](https://doi.org/10.1109/TE.2023.3281825)
* Arjit Seth and Rhea Liem. *AeroMDAO — A Multidisciplinary Aircraft Design Platform for Education*. Teaching and Learning Symposium, Center for Education Innovation, HKUST, June 2022. [Talk](https://youtu.be/_H5ig2tr7S4)
