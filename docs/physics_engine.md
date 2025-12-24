# Physics Engine Technical Reference

The core of the simulation is the relativistic integration of null geodesics (light paths) around a Schwarzschild black hole.

## 1. Schwarzschild Spacetime

The metric $g_{\mu\nu}$ for a non-rotating, uncharged mass $M$ in Schwarzschild coordinates $(t, r, \theta, \phi)$ is:

$$ ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2d\theta^2 + r^2\sin^2\theta d\phi^2 $$

In this simulation, we use geometric units where $G = c = M = 1$.

## 2. Hamiltonian Formulation

To evolve rays, we use the Hamiltonian $H$:

$$ H = \frac{1}{2}g^{\mu\nu} p_\mu p_\nu $$

For photons (null geodesics), the physical constraint is $H = 0$.

### Equations of Motion

Hamilton's equations provide the derivatives with respect to the affine parameter $\lambda$:

1.  **Coordinate evolution**:
    $$\frac{dx^\mu}{d\lambda} = \frac{\partial H}{\partial p_\mu} = g^{\mu\nu} p_\nu$$

2.  **Momentum evolution**:
    $$\frac{dp_\mu}{d\lambda} = -\frac{\partial H}{\partial x^\mu} = -\frac{1}{2}\left(\frac{\partial g^{\alpha\beta}}{\partial x^\mu}\right)p_\alpha p_\beta$$

## 3. Implementation Details

### Analytic Derivatives

To avoid the overhead and instability of finite differences, we implement explicit analytic derivatives for the contravariant metric $\partial_\mu g^{\alpha\beta}$.

For the Schwarzschild metric, the non-zero derivatives are primarily with respect to $r$ and $\theta$:

*   $\partial_r g^{tt} = \frac{2M}{r^2 (1-2M/r)^2}$
*   $\partial_r g^{rr} = \frac{2M}{r^2}$
*   $\partial_r g^{\theta\theta} = -\frac{2}{r^3}$
*   $\partial_r g^{\phi\phi} = -\frac{2}{r^3 \sin^2\theta}$
*   $\partial_\theta g^{\phi\phi} = -\frac{2 \cot\theta}{r^2 \sin^2\theta}$

### Constraint Stabilization

Numerical integration (RK4) accumulates error over time, causing $|H|$ to drift from zero. We apply a projection method at each step:

1.  Compute the current kinetic part of the Hamiltonian.
2.  Rescale the momentum components $p_r$ and $p_\theta$ to force $H=0$ while strictly preserving the conserved quantities $E$ and $L$.

### Singularity Handling

The Schwarzschild coordinate system is singular at the event horizon ($r=2M$). To avoid numerical blow-up, the integrator terminates if $r < 2.1M$.

## 4. Conserved Quantities

*   **Energy ($E$)**: Since the metric is static ($\partial_t g_{\mu\nu} = 0$), $E = -p_t$ is strictly conserved.
*   **Angular Momentum ($L$)**: Since the metric is axially symmetric ($\partial_\phi g_{\mu\nu} = 0$), $L = p_\phi$ is strictly conserved.
