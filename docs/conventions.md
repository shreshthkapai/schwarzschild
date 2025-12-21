# Physical Conventions & Units

## Geometric Units
We use **geometric units** where:
- G = 1 (gravitational constant)
- c = 1 (speed of light)
- M = 1 (black hole mass)

All quantities are dimensionless in these units.

## Schwarzschild Metric
Line element in Schwarzschild coordinates (t, r, θ, φ):

ds² = -(1 - 2M/r) dt² + (1 - 2M/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²

In our units (M=1):

ds² = -(1 - 2/r) dt² + (1 - 2/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²

## Critical Radii
- **Event Horizon**: r_s = 2M = 2
- **Photon Sphere**: r_ph = 3M = 3
- **ISCO** (timelike): r_isco = 6M = 6

## Coordinate System
We evolve geodesics in **Schwarzschild coordinates** (t, r, θ, φ):
- t: Schwarzschild time
- r: radial coordinate (areal radius)
- θ: polar angle [0, π]
- φ: azimuthal angle [0, 2π]

For **visualization**, we convert to Cartesian:
- x = r sin(θ) cos(φ)
- y = r sin(θ) sin(φ)
- z = r cos(θ)

## Affine Parameter
Null geodesics are parameterized by **affine parameter λ**:
- Not physical time
- Monotonically increases along ray
- Natural parameter for geodesic equation

## Hamiltonian Formulation
Phase space: (x^μ, p_μ) where μ ∈ {t, r, θ, φ}

Hamiltonian constraint (null condition):
H = (1/2) g^μν p_μ p_ν = 0

For null geodesics in Schwarzschild:
H = -(1 - 2/r)⁻¹ p_t² + (1 - 2/r) p_r² + (1/r²) p_θ² + (1/(r² sin²θ)) p_φ²

## Constants of Motion
Due to symmetries:
- **Energy**: E = -p_t (time translation symmetry)
- **Angular momentum**: L = p_φ (axial symmetry)

These are conserved along geodesics.

## Numerical Conventions
- Default integration: RK4 with adaptive step size consideration
- Step size: Δλ ~ 0.01 (tunable)
- Max affine parameter: λ_max ~ 100 (tunable)
- Constraint violation tolerance: |H| < 10⁻⁶

## Sign Conventions
- Metric signature: (-,+,+,+)
- Indices: Greek letters μ,ν run over {0,1,2,3} = {t,r,θ,φ}
- Lowering/raising: p^μ = g^μν p_ν