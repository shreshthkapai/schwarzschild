# Physical Conventions & Units

## Geometric Units
We use **geometric units** where:
- $G = 1$ (gravitational constant)
- $c = 1$ (speed of light)
- $M = 1$ (black hole mass)

Consequently, all physical quantities (length, time, mass) are dimensionless. For example, a radius of $r=2$ corresponds to $2GM/c^2$ in SI units.

## Schwarzschild Metric
The line element in Schwarzschild coordinates $(t, r, \theta, \phi)$ is:

$$ ds^2 = -\left(1 - \frac{2}{r}\right)dt^2 + \left(1 - \frac{2}{r}\right)^{-1}dr^2 + r^2d\theta^2 + r^2\sin^2\theta d\phi^2 $$

## Critical Radii
- **Event Horizon**: $r_s = 2M = 2$
- **Photon Sphere**: $r_{ph} = 3M = 3$ (unstable circular photon orbit)
- **ISCO** (Massive particle): $r_{isco} = 6M = 6$
- **Accretion Disk Inner Edge**: Defaults to ISCO ($r=6$).

## Coordinate System
Geodesics are evolved in **Schwarzschild coordinates**:
- $t$: Schwarzschild time (non-affine coordinate)
- $r$: Radial coordinate (areal radius)
- $\theta$: Polar angle $[0, \pi]$ (Equatorial plane is $\theta = \pi/2$)
- $\phi$: Azimuthal angle $[0, 2\pi]$

### Cartesian Conversion for Visualization:
- $x = r \sin\theta \cos\phi$
- $y = r \cos\theta$
- $z = r \sin\theta \sin\phi$

## Affine Parameter
Null geodesics are parameterized by the **affine parameter $\lambda$**. Unlike timelike geodesics where $\lambda$ usually corresponds to proper time $\tau$, for photons $\lambda$ is an arbitrary monotonically increasing parameter that allows the geodesic equation to be written as $d^2x^\mu/d\lambda^2 + \Gamma^\mu_{\alpha\beta} dx^\alpha/d\lambda dx^\beta/d\lambda = 0$.

## Hamiltonian Formulation
Phase space consists of coordinates $x^\mu$ and conjugate momenta $p_\mu$:
- $(x^\mu, p_\mu) = (t, r, \theta, \phi, p_t, p_r, p_\theta, p_\phi)$

The null constraint is:
$$ H = \frac{1}{2}g^{\mu\nu} p_\mu p_\nu = 0 $$

## Constants of Motion
Due to symmetries in the Schwarzschild metric:
- **Energy**: $E = -p_t$ (strictly conserved)
- **Angular Momentum**: $L = p_\phi$ (strictly conserved)

## Numerical Parameters
- **Integrator**: 4th-order Runge Kutta (RK4)
- **Step Size**: $\Delta\lambda \sim 0.05$ (Scaled by proximity to horizon)
- **Max Path**: $\lambda_{max} = 100$
- **Constraint Tolerance**: $|H| < 10^{-6}$
- **Metric Signature**: $(-,+,+,+)$