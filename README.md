# Focus Scientific ML
An advanced Python implementation of dimensionality reduction and hyper-reduction  techniques applied to non-linear Partial Differential Equations (PDEs).  


🚀 Hyper-Reduced PDE Solver: ROM & DEIM Optimization
Accelerating Fluid Dynamics from Seconds to Milliseconds
This repository implements an advanced Model Order Reduction (ROM) framework to solve non-linear Partial Differential Equations (PDEs). By leveraging Proper Orthogonal Decomposition (POD) and Discrete Empirical Interpolation (DEIM), this solver transforms high-dimensional physical problems into low-rank systems capable of real-time execution.

💡 Why This Project Matters
In industries such as Aerospace, Automotive, and Energy, high-fidelity simulations (CFD) are often too slow for real-time control, digital twins, or large-scale optimization. The Solution: This project bridges the gap by projecting the physics onto a reduced subspace, maintaining 99%+ accuracy while reducing computational costs by several orders of magnitude.

🛠️ Technical Stack & Algorithms

1. Full Order Model (FOM):
- Finite Volume Method (FVM) discretization.
- Slope Limiters implementation for shock-capturing without numerical oscillations.
- Time integration via 2nd-order Runge-Kutta (Heun’s Method).

2. Dimensionality Reduction:
- POD (Proper Orthogonal Decomposition) via SVD to extract dominant physical modes.

Greedy Basis Algorithm for adaptive snapshot selection.

3. Hyper-Reduction (The Core):

- DEIM (Discrete Empirical Interpolation Method): An advanced technique to evaluate non-linear fluxes at a few selected interpolation points. This ensures the computational cost scales with the reduced dimension rather than the full mesh size.

📊 Performance Benchmark
The reduced model successfully captures complex physics (like shock wave propagation) while significantly decreasing CPU time.

**TODO

Quick Start
Python

from main import plot_reduced

# Simulate Burgers' equation with 10 POD modes and DEIM acceleration
plot_reduced(
    model="Burgers", 
    Re=100, 
    n_modes=10, 
    red="POD", 
    hyp="DEIM", 
    affine=True
)

🧠 Demonstrated ExpertiseNumerical Analysis: Spatio-temporal discretization, scheme stability, and flux limiters.Advanced Linear Algebra: Low-rank approximations, SVD projection, and empirical interpolation.Optimization: Reducing algorithmic complexity from $O(N)$ to $O(K)$, where $K \ll N$.Scientific Programming: Expert use of NumPy, SciPy, and Matplotlib for high-performance computing.

Dimensionality Reduction (POD)The goal is to approximate the high-dimensional state $u(t) \in \mathbb{R}^N$ by a linear combination of $K$ basis vectors (modes) $\Phi \in \mathbb{R}^{N \times K}$:$$u(t) \approx \bar{u} + \Phi \hat{u}(t)$$where $\Phi$ is obtained via Singular Value Decomposition (SVD) of the snapshot matrix $S$:$$S = U \Sigma V^T \implies \Phi = U_{1:K}$$2. Discrete Empirical Interpolation Method (DEIM)For non-linear systems $\dot{u} = F(u)$, projecting the equations $ \dot{\hat{u}} = \Phi^T F(\bar{u} + \Phi \hat{u})$ still requires evaluating the full vector $F \in \mathbb{R}^N$. DEIM approximates the non-linear term using a secondary basis $\Psi$ and a selection matrix $P$:$$\hat{F}(u) \approx \Psi (P^T \Psi)^{-1} P^T F(u)$$This allows evaluating $F(u)$ only at $M$ selected points ($M \ll N$), providing the massive speed-up.3. Governing Equation (Burgers)The solver focuses on the viscous Burgers' equation:$$\frac{\partial \rho}{\partial t} + \rho \frac{\partial \rho}{\partial x} = \frac{1}{Re} \frac{\partial^2 \rho}{\partial x^2}$$discretized using a Finite Volume scheme:$$\rho_i^{n+1} = \rho_i^n - \frac{\Delta t}{h} (G_{i+1/2} - G_{i-1/2})$$
