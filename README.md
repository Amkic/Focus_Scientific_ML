Example: Reduced Order Modeling for a Nonlinear PDE (Burgers Equation)

This repository provides a working implementation of a reduced-order modeling (ROM) pipeline applied to a nonlinear partial differential equation (the 1D Burgers equation). The objective is to demonstrate how high-fidelity numerical simulations can be accelerated using projection-based model reduction techniques, while preserving good accuracy.

The code is intentionally written in a clear and modular way to highlight the full workflow used in industrial and scientific applications.

What this example demonstrates

This implementation shows how to:

- build a high-fidelity finite-volume solver (FOM) for a nonlinear PDE,
- generate solution snapshots over time,
- construct a reduced basis using Proper Orthogonal Decomposition (POD) or a Greedy algorithm,
- accelerate nonlinear evaluations using DEIM (Discrete Empirical Interpolation Method),
- project the governing equations onto a low-dimensional subspace,
- compare full-order and reduced-order solutions.

The reduced model is then solved using the same time integration scheme as the full model.


The example considers the one-dimensional viscous Burgers equation:

$$\partial_t \rho + \partial_x f(\rho) = \frac{1}{\mathrm{Re}}\, \partial_{xx} \rho,$$

with either:

linear advection: 
$$f(u) = u$$

nonlinear Burgers flux:

$$f(u) = \frac{1}{2}u^2$$.

The Reynolds number controls the diffusion intensity and allows testing different regimes.

Numerical method (Full Order Model)

The full-order solver uses:

- finite volume discretization on a uniform grid
- ghost cells for boundary conditions
- MUSCL-type reconstruction with slope limiting
- explicit second-order Runge–Kutta (Heun) time integration
- numerical diffusion for robustness

This represents a realistic CFD-style solver similar in spirit to industrial simulation codes.

Reduced Order Modeling strategy
POD (Proper Orthogonal Decomposition)

Snapshots of the solution are collected during the full simulation:

$$S = u(t_1), \dots, u(t_N)$$

A reduced basis is extracted via SVD:

$$S = U \Sigma V^T$$

The first n_modes dominant modes are kept, capturing most of the system’s energy.

The solution is approximated as:

$$u(x,t) \approx \Phi a(t) + \bar{u}$$

where:
- $\Phi$ is the reduced basis matrix,
- $a(t)$ are the reduced coordinates,
- $\bar{u}$ is a reference (mean) state.

Hyper-reduction with DEIM

For nonlinear problems, evaluating the full nonlinear flux still scales with the original dimension.
To overcome this, the Discrete Empirical Interpolation Method (DEIM) is used.

DEIM approximates the nonlinear term as:


$$f(u) \approx \Phi_f (P^T \Phi_f)^{-1} P^T f(u)$$

where:

$\Phi_f$ is a reduced basis for the nonlinear term, $P$ is a sampling matrix selecting interpolation indices.

This reduces the computational complexity from:

$$\mathcal{O}(N) \;\longrightarrow\; \mathcal{O}(r)$$

making real-time or many-query simulations feasible.

Reduced dynamical system

Without hyper-reduction:

$$\dot{a} = \Phi^T F(\Phi a + \bar{u})$$


With DEIM, it becomes:

$$\dot{a} = \Phi^T \Pi_{\mathrm{DEIM}} F(\Phi a + \bar{u})$$,

$$\Pi_{\mathrm{DEIM}} = \Phi_f (P^T \Phi_f)^{-1} P^T$$.

The reduced system is then integrated in time using the same explicit Runge–Kutta scheme.

Example result

Below is a typical outcome of the code:

-Red curve: full-order solution (reference)

-Blue curve: reduced-order reconstruction.

Only a small number of modes (e.g. 5–10) is sufficient to accurately recover the solution shape

<img width="1400" height="500" alt="Figure_1" src="https://github.com/user-attachments/assets/bcd051f6-d41e-4dcb-8c36-bf27d2fd8e0e" />


This illustrates that:

✅ the dominant dynamics lie in a low-dimensional subspace

✅ reduced models can reproduce nonlinear PDE behavior accurately

✅ large computational savings are achievable

Why this matters in industry ?

This type of methodology is widely used in:

- digital twins
- real-time simulation and control
- parametric studies and design optimization
- uncertainty quantification
- reduced CFD models
- physics-informed AI pipelines
- surrogate modeling
- hybrid physics–ML systems

ROM techniques allow engineers to replace expensive simulations with fast surrogate models while preserving physical consistency.
