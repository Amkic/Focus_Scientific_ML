POD, Greedy and DEIM applied to the 1D Burgers Equation
📌 Overview

This repository presents a complete implementation of projection-based reduced order modeling (ROM) techniques applied to a nonlinear partial differential equation: the 1D Burgers equation.

The project combines:

a high-fidelity finite volume solver (FOM),

snapshot-based model reduction,

Proper Orthogonal Decomposition (POD),

Greedy reduced basis construction,

hyper-reduction using DEIM,

and a quantitative comparison between full and reduced models.

The goal is to achieve significant computational speed-up while preserving accuracy, a key challenge in scientific computing, digital twins, and physics-informed AI.

🧠 Mathematical model

We consider the nonlinear conservation law:

$$\begin{equation}
\partial_t \rho + \partial_x f(\rho)
= \frac{1}{Re}\,\partial_{xx}\rho,
\end{equation}$$

defined on a one-dimensional spatial domain with suitable boundary conditions.

Flux functions

Two fluxes are supported:

$$\begin{equation}
f(u) = u
\end{equation}$$

$$\begin{equation}
f(u) = \frac{1}{2}u^2
\end{equation}$$

🔢 Full Order Model (FOM)
Spatial discretization

The equation is discretized using a finite volume method on a uniform grid:

- cell-centered unknowns
- ghost cells for boundary conditions
- MUSCL-type reconstruction
- slope limiter for stability
- numerical diffusion for robustness
- The numerical flux has the general form:

$$\begin{equation}
F_{i+\frac12}
=
\frac{f(u_L)+f(u_R)}{2}
- \frac{\lambda}{2}(u_L-u_R)
- \nu \nabla u,
\end{equation}$$

where:
$$\begin{itemize}
    \item $u_L, u_R$ are reconstructed interface states,
    \item $\lambda = \max |f'(u)|$,
    \item $\nu = \frac{1}{Re}$.
\end{itemize}$$

Time integration

A second-order explicit Runge–Kutta scheme (Heun method) is used:

$$\begin{align}
u^{*} &= u^n + \frac{\Delta t}{2} F(u^n), \\
u^{n+1} &= u^n + \Delta t\, F(u^{*}).
\end{align}$$

The timestep satisfies a CFL-like condition:

$$\begin{equation}
\Delta t = 0.4 \min(h, Re\, h^2).
\end{equation}$$

📸 Snapshot generation

During the full-order simulation, solution snapshots are collected:

$$\begin{equation}
S = [u(t_1), u(t_2), \dots, u(t_N)] \in \mathbb{R}^{N_x \times N_t}.
\end{equation}$$


These snapshots form the basis for reduced-order modeling.

📉 Proper Orthogonal Decomposition (POD)

Snapshots are decomposed using Singular Value Decomposition:

$$\begin{equation}
S = U \Sigma V^T.
\end{equation}$$

The reduced basis is defined as:
$$\begin{equation}
\Phi = U_{(:,1:r)},
\end{equation}
where $r \ll N_x$.$$

⚡ Hyper-reduction with DEIM

For nonlinear problems, evaluating the full nonlinear flux remains expensive.

To address this, the Discrete Empirical Interpolation Method (DEIM) is used.

DEIM approximation

For nonlinear problems, evaluating the full nonlinear term is computationally expensive.
The \textbf{Discrete Empirical Interpolation Method (DEIM)} alleviates this cost.

Let $\Phi_f$ be POD modes of the nonlinear flux. The DEIM approximation reads:
$$\begin{equation}
f(u) \approx \Phi_f (P^T \Phi_f)^{-1} P^T f(u),
\end{equation}$$
where $P$ is a sparse selection matrix extracting a few spatial entries.

This reduces the complexity of nonlinear evaluations from $\mathcal{O}(N)$ to $\mathcal{O}(r)$.

🧩 Reduced-order dynamical system.

The reduced solution is written as:
$$\begin{equation}
u(x,t) \approx \Phi a(t) + \bar{u}.
\end{equation}$$

Without hyper-reduction:
$$\begin{equation}
\dot{a} = \Phi^T F(\Phi a + \bar{u}).
\end{equation}$$

With DEIM : 
$$\begin{equation}
\dot{a} = \Phi^T \Pi_{\mathrm{DEIM}} F(\Phi a + \bar{u}),
\end{equation}$$
with
$$\begin{equation}
\Pi_{\mathrm{DEIM}} = \Phi_f (P^T \Phi_f)^{-1} P^T.
\end{equation}$$

Evaluation

The code automatically:

- computes the full-order solution,
- builds reduced bases (POD or Greedy)
- applies DEIM hyper-reduction
- solves the reduced model
- measures execution time
- computes the final-time error

$$|u_{\text{ROM}}(T) - u_{\text{FOM}}(T)\|$$

- plots the reduced and full solutions
- visualizes POD modes


This project demonstrates:

- solid background in numerical PDEs
- projection-based reduced order modeling
- POD–Galerkin and DEIM methods
- scientific Python implementation
- understanding of accuracy–performance trade-offs
- relevance for real-time simulation and digital twins

Technologies :

- Python
- NumPy / SciPy
- Matplotlib
- Numerical linear algebra
- Model Order Reduction
- Possible extensions
- Stabilized POD–Galerkin
- GNAT hyper-reduction
- Autoencoder-based ROM
- Neural operators
- Parametric ROM
- 2D Burgers or Navier–Stokes
- A posteriori error estimation


