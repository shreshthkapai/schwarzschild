# Hamiltonian Formulation of Null Geodesics

## Phase Space
We work in 8-dimensional phase space: (x^μ, p_μ) where μ ∈ {t, r, θ, φ}

## Hamiltonian for Geodesics
The Hamiltonian for geodesic motion is:

H = (1/2) g^μν p_μ p_ν

For **null geodesics** (photons), the constraint is:
H = 0

This is equivalent to the null condition: g_μν dx^μ/dλ dx^ν/dλ = 0

## Hamilton's Equations
The geodesic equations are derived from:

**Position evolution:**
dx^μ/dλ = ∂H/∂p_μ = g^μν p_ν

**Momentum evolution:**
dp_μ/dλ = -∂H/∂x^μ = -(1/2) (∂g^αβ/∂x^μ) p_α p_β

## Schwarzschild Specifics
For Schwarzschild metric (diagonal):
- g^μν only has diagonal components
- Position derivatives simplify: dx^μ/dλ = g^μμ p_μ (no sum!)

## Explicit Equations (Schwarzschild)

Let f(r) = 1 - 2M/r. Then:

**Position derivatives:**
```
dt/dλ = g^tt p_t = -(1/f) p_t
dr/dλ = g^rr p_r = f p_r
dθ/dλ = g^θθ p_θ = (1/r²) p_θ
dφ/dλ = g^φφ p_φ = (1/(r² sin²θ)) p_φ
```

**Momentum derivatives:**
```
dp_t/dλ = 0  (time translation symmetry → E conserved)

dp_r/dλ = -(1/2) [∂g^tt/∂r p_t² + ∂g^rr/∂r p_r² + ∂g^θθ/∂r p_θ² + ∂g^φφ/∂r p_φ²]

dp_θ/dλ = -(1/2) ∂g^φφ/∂θ p_φ²

dp_φ/dλ = 0  (axial symmetry → L conserved)
```

## Constants of Motion
Due to Killing vectors:
- **Energy**: E = -p_t (from ∂_t symmetry)
- **Angular momentum**: L = p_φ (from ∂_φ symmetry)

These are **exactly conserved** along geodesics.

## Equatorial Plane Reduction
For motion starting in equatorial plane (θ = π/2, p_θ = 0):
- Motion remains in equatorial plane (p_θ = 0 preserved)
- θ = π/2 for entire trajectory
- Reduces 8D system → 6D effective system

This is what we'll implement for the initial version.

## Numerical Integration
We integrate the coupled system:
- State vector: [t, r, θ, φ, p_t, p_r, p_θ, p_φ]
- RHS given by Hamilton's equations above
- Use RK4 integrator with affine parameter λ

## Constraint Monitoring
At each step, verify:
- |H| < tolerance (null condition)
- |E - E_0| < tolerance (energy conservation)
- |L - L_0| < tolerance (angular momentum conservation)

If constraints violated → numerical error accumulating