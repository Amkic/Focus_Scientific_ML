This project presents a complete implementation of \textbf{projection-based reduced order modeling (ROM)} techniques applied to a nonlinear partial differential equation, namely the \textbf{1D Burgers equation}.

The implementation combines:
\begin{itemize}
    \item a high-fidelity finite volume solver (Full Order Model -- FOM),
    \item snapshot-based dimensionality reduction,
    \item Proper Orthogonal Decomposition (POD),
    \item Greedy reduced basis construction,
    \item hyper-reduction using the Discrete Empirical Interpolation Method (DEIM),
    \item quantitative comparison between full and reduced models.
\end{itemize}

Such techniques are widely used in scientific computing, reduced-order modeling, digital twins, and physics-informed machine learning.

\section{Mathematical model}

We consider the nonlinear conservation law
\begin{equation}
\partial_t \rho + \partial_x f(\rho)
= \frac{1}{Re}\,\partial_{xx}\rho,
\end{equation}
defined on a one-dimensional spatial domain with suitable boundary conditions.

\subsection*{Flux functions}

Two flux functions are implemented:

\paragraph{Linear advection}
\begin{equation}
f(u) = u
\end{equation}

\paragraph{Burgers equation}
\begin{equation}
f(u) = \frac{1}{2}u^2
\end{equation}

\section{Full Order Model (FOM)}

\subsection{Spatial discretization}

The equation is discretized using a \textbf{finite volume method} on a uniform mesh:

\begin{itemize}
    \item cell-centered unknowns,
    \item ghost cells for boundary conditions,
    \item MUSCL-type reconstruction,
    \item slope limiter for stability,
    \item numerical viscosity for robustness.
\end{itemize}

The numerical flux at an interface reads
\begin{equation}
F_{i+\frac12}
=
\frac{f(u_L)+f(u_R)}{2}
- \frac{\lambda}{2}(u_L-u_R)
- \nu \nabla u,
\end{equation}
where:
\begin{itemize}
    \item $u_L, u_R$ are reconstructed interface states,
    \item $\lambda = \max |f'(u)|$,
    \item $\nu = \frac{1}{Re}$.
\end{itemize}

\subsection{Time integration}

A second-order explicit Runge--Kutta (Heun) scheme is used:
\begin{align}
u^{*} &= u^n + \frac{\Delta t}{2} F(u^n), \\
u^{n+1} &= u^n + \Delta t\, F(u^{*}).
\end{align}

The time step satisfies the stability condition:
\begin{equation}
\Delta t = 0.4 \min(h, Re\, h^2).
\end{equation}

\section{Snapshot generation}

During the time integration, solution snapshots are collected:
\begin{equation}
S = [u(t_1), u(t_2), \dots, u(t_N)] \in \mathbb{R}^{N_x \times N_t}.
\end{equation}

These snapshots form the database for reduced-order modeling.

\section{Proper Orthogonal Decomposition (POD)}

Snapshots are decomposed using Singular Value Decomposition:
\begin{equation}
S = U \Sigma V^T.
\end{equation}

The reduced basis is defined as:
\begin{equation}
\Phi = U_{(:,1:r)},
\end{equation}
where $r \ll N_x$.

\subsection{Affine decomposition}

To improve approximation quality, an affine offset is introduced:
\begin{equation}
\bar{u} = \frac{1}{N_t} \sum_{k=1}^{N_t} u(t_k).
\end{equation}

The reduced representation becomes:
\begin{equation}
u(x,t) \approx \Phi a(t) + \bar{u}.
\end{equation}

\section{Greedy reduced basis construction}

A greedy algorithm is also implemented:

\begin{enumerate}
    \item Initialize the basis with one snapshot.
    \item Project all snapshots onto the current basis.
    \item Compute the projection error.
    \item Select the snapshot with maximal error.
    \item Orthonormalize and enrich the basis.
    \item Repeat until the desired dimension is reached.
\end{enumerate}

This approach is commonly used in classical reduced basis methods.

\section{Hyper-reduction with DEIM}

For nonlinear problems, evaluating the full nonlinear term is computationally expensive.
The \textbf{Discrete Empirical Interpolation Method (DEIM)} alleviates this cost.

Let $\Phi_f$ be POD modes of the nonlinear flux. The DEIM approximation reads:
\begin{equation}
f(u) \approx \Phi_f (P^T \Phi_f)^{-1} P^T f(u),
\end{equation}
where $P$ is a sparse selection matrix extracting a few spatial entries.

This reduces the complexity of nonlinear evaluations from $\mathcal{O}(N)$ to $\mathcal{O}(r)$.

\section{Reduced-order dynamical system}

The reduced solution is written as:
\begin{equation}
u(x,t) \approx \Phi a(t) + \bar{u}.
\end{equation}

\subsection*{Without hyper-reduction}
\begin{equation}
\dot{a} = \Phi^T F(\Phi a + \bar{u}).
\end{equation}

\subsection*{With DEIM}
\begin{equation}
\dot{a} = \Phi^T \Pi_{\mathrm{DEIM}} F(\Phi a + \bar{u}),
\end{equation}
with
\begin{equation}
\Pi_{\mathrm{DEIM}} = \Phi_f (P^T \Phi_f)^{-1} P^T.
\end{equation}

\section{Time integration of the reduced system}

The reduced system is integrated using the same second-order Runge--Kutta scheme:
\begin{equation}
a^{n+1} = a^n + \Delta t\, \Phi^T F(\cdot).
\end{equation}

This ensures a fair comparison between the full and reduced models.

\section{Evaluation and outputs}

The code automatically:
\begin{itemize}
    \item computes the full-order solution,
    \item builds reduced bases (POD or Greedy),
    \item applies DEIM hyper-reduction,
    \item solves the reduced system,
    \item measures execution time,
    \item computes the final-time error
    \[
    \|u_{\text{ROM}}(T) - u_{\text{FOM}}(T)\|,
    \]
    \item plots solution comparisons and POD modes.
\end{itemize}

\section{Typical performance}

Example output:
\begin{verbatim}
Time FOM: 2.31 s
Time ROM: 0.08 s
SPEED-UP: 28.9x
\end{verbatim}

\section{Why this project matters}

This project demonstrates:
\begin{itemize}
    \item strong foundations in numerical analysis and PDEs,
    \item mastery of projection-based reduced order models,
    \item POD--Galerkin and DEIM techniques,
    \item efficient scientific Python implementation,
    \item understanding of stability and accuracy trade-offs,
    \item relevance to real-time simulation and digital twins.
\end{itemize}

\section{Technologies}

\begin{itemize}
    \item Python
    \item NumPy / SciPy
    \item Matplotlib
    \item Numerical linear algebra
    \item Model Order Reduction
\end{itemize}

\section{Possible extensions}

\begin{itemize}
    \item Stabilized POD--Galerkin methods
    \item GNAT hyper-reduction
    \item Autoencoder-based ROMs
    \item Neural operators
    \item Parametric ROMs
    \item 2D Burgers or Navier--Stokes equations
    \item A posteriori error estimation
    \item Comparison with Deep Galerkin Method (DGM)
\end{itemize}

\section*{Author}

\textbf{Amer Jukic}\\
MSc in Mathematics -- PDEs \& Deep Learning\\
Software Engineer (Python, C\#, scientific computing)
