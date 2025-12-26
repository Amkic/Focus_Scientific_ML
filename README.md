🚀 Reduced Order Modeling for Nonlinear PDEs
POD, Greedy and DEIM applied to the 1D Burgers equation
Overview

This repository presents a complete implementation of projection-based reduced order modeling (ROM) techniques applied to a nonlinear partial differential equation: the 1D Burgers equation.

The project combines:

a full-order finite volume solver (FOM),

snapshot-based dimensionality reduction,

Proper Orthogonal Decomposition (POD),

Greedy reduced basis construction,

hyper-reduction using the Discrete Empirical Interpolation Method (DEIM),

comparison between accuracy and computational cost.

These techniques are widely used in scientific computing, reduced-order modeling, digital twins and physics-informed AI.

Mathematical model

We consider the nonlinear conservation law

∂
𝑡
𝜌
+
∂
𝑥
𝑓
(
𝜌
)
=
1
𝑅
𝑒
 
∂
𝑥
𝑥
𝜌
∂
t
	​

ρ+∂
x
	​

f(ρ)=
Re
1
	​

∂
xx
	​

ρ

defined on a one-dimensional spatial domain.

Flux functions

Linear advection

𝑓
(
𝑢
)
=
𝑢
f(u)=u

Burgers equation

𝑓
(
𝑢
)
=
1
2
𝑢
2
f(u)=
2
1
	​

u
2
Full Order Model (FOM)
Spatial discretization

The equation is discretized using a finite volume method with:

cell-centered unknowns

ghost cells for boundary conditions

MUSCL-type reconstruction

slope limiter for stability

numerical viscosity

The numerical flux reads

𝐹
𝑖
+
1
2
=
𝑓
(
𝑢
𝐿
)
+
𝑓
(
𝑢
𝑅
)
2
−
𝜆
2
(
𝑢
𝐿
−
𝑢
𝑅
)
−
𝜈
∇
𝑢
F
i+
2
1
	​

	​

=
2
f(u
L
	​

)+f(u
R
	​

)
	​

−
2
λ
	​

(u
L
	​

−u
R
	​

)−ν∇u

where:

$u_L, u_R$ are reconstructed interface values

$\lambda = \max |f'(u)|$

$\nu = \frac{1}{Re}$

Time integration

A second-order explicit Runge–Kutta (Heun) scheme is used:

𝑢
∗
=
𝑢
𝑛
+
Δ
𝑡
2
𝐹
(
𝑢
𝑛
)
u
∗
=u
n
+
2
Δt
	​

F(u
n
)
𝑢
𝑛
+
1
=
𝑢
𝑛
+
Δ
𝑡
 
𝐹
(
𝑢
∗
)
u
n+1
=u
n
+ΔtF(u
∗
)

The time step is chosen as

Δ
𝑡
=
0.4
min
⁡
(
ℎ
,
𝑅
𝑒
 
ℎ
2
)
Δt=0.4min(h,Reh
2
)
Snapshot generation

During the full-order simulation, solution snapshots are collected:

𝑆
=
[
𝑢
(
𝑡
1
)
,
𝑢
(
𝑡
2
)
,
…
,
𝑢
(
𝑡
𝑁
)
]
∈
𝑅
𝑁
𝑥
×
𝑁
𝑡
S=[u(t
1
	​

),u(t
2
	​

),…,u(t
N
	​

)]∈R
N
x
	​

×N
t
	​


These snapshots are used to construct reduced bases.

Proper Orthogonal Decomposition (POD)

Snapshots are decomposed using Singular Value Decomposition:

𝑆
=
𝑈
Σ
𝑉
𝑇
S=UΣV
T

The reduced basis is defined as:

Φ
=
𝑈
(
:
,
1
:
𝑟
)
Φ=U
(:,1:r)
	​


where $r \ll N_x$.

Affine decomposition (optional)

To improve accuracy, an affine offset is introduced:

𝑢
ˉ
=
1
𝑁
𝑡
∑
𝑘
=
1
𝑁
𝑡
𝑢
(
𝑡
𝑘
)
u
ˉ
=
N
t
	​

1
	​

k=1
∑
N
t
	​

	​

u(t
k
	​

)

The reduced approximation becomes:

𝑢
(
𝑥
,
𝑡
)
≈
Φ
𝑎
(
𝑡
)
+
𝑢
ˉ
u(x,t)≈Φa(t)+
u
ˉ
Greedy reduced basis

A greedy algorithm is also implemented:

Initialize the basis with one snapshot

Project all snapshots onto the current basis

Compute the projection error

Select the snapshot with maximal error

Orthonormalize and enrich the basis

Repeat until the desired dimension is reached

Hyper-reduction with DEIM

For nonlinear problems, evaluating the full nonlinear term is expensive.
The Discrete Empirical Interpolation Method (DEIM) is therefore used.

Let $\Phi_f$ be POD modes of the nonlinear flux. The DEIM approximation is

𝑓
(
𝑢
)
≈
Φ
𝑓
(
𝑃
𝑇
Φ
𝑓
)
−
1
𝑃
𝑇
𝑓
(
𝑢
)
f(u)≈Φ
f
	​

(P
T
Φ
f
	​

)
−1
P
T
f(u)

where $P$ is a sparse selection matrix.

This reduces the computational cost of nonlinear evaluations from
$\mathcal{O}(N)$ to $\mathcal{O}(r)$.

Reduced-order dynamical system

The reduced solution is written as:

𝑢
(
𝑥
,
𝑡
)
≈
Φ
𝑎
(
𝑡
)
+
𝑢
ˉ
u(x,t)≈Φa(t)+
u
ˉ
Without hyper-reduction
𝑎
˙
=
Φ
𝑇
𝐹
(
Φ
𝑎
+
𝑢
ˉ
)
a
˙
=Φ
T
F(Φa+
u
ˉ
)
With DEIM
𝑎
˙
=
Φ
𝑇
Π
D
E
I
M
𝐹
(
Φ
𝑎
+
𝑢
ˉ
)
a
˙
=Φ
T
Π
DEIM
	​

F(Φa+
u
ˉ
)

with

Π
D
E
I
M
=
Φ
𝑓
(
𝑃
𝑇
Φ
𝑓
)
−
1
𝑃
𝑇
Π
DEIM
	​

=Φ
f
	​

(P
T
Φ
f
	​

)
−1
P
T
Time integration of the reduced system

The reduced system is integrated using the same RK2 scheme:

𝑎
𝑛
+
1
=
𝑎
𝑛
+
Δ
𝑡
 
Φ
𝑇
𝐹
(
⋅
)
a
n+1
=a
n
+ΔtΦ
T
F(⋅)
Evaluation

The code automatically:

computes the full-order solution,

builds reduced bases (POD or Greedy),

applies DEIM hyper-reduction,

solves the reduced model,

measures execution time,

computes the final-time error

∥
𝑢
ROM
(
𝑇
)
−
𝑢
FOM
(
𝑇
)
∥
∥u
ROM
	​

(T)−u
FOM
	​

(T)∥

plots the reduced and full solutions

visualizes POD modes

Why this project matters

This project demonstrates:

solid background in numerical PDEs

projection-based reduced order modeling

POD–Galerkin and DEIM methods

scientific Python implementation

understanding of accuracy–performance trade-offs

relevance for real-time simulation and digital twins

Technologies

Python

NumPy / SciPy

Matplotlib

Numerical linear algebra

Model Order Reduction

Possible extensions

Stabilized POD–Galerkin

GNAT hyper-reduction

Autoencoder-based ROM

Neural operators

Parametric ROM

2D Burgers or Navier–Stokes

A posteriori error estimation

Comparison with Deep Galerkin Method

Author

Amer Jukic
MSc in Mathematics — PDEs & Deep Learning
Software Engineer (Python, C#, scientific computing)
