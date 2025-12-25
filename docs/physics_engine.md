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

Explicit analytic derivatives for the contravariant metric $\partial_\mu g^{\alpha\beta}$ are used to ensure stability and performance. For the Schwarzschild metric, the non-zero derivatives computed are:

*   $\partial_r g^{tt} = \frac{2M}{r^2 (1-2M/r)^2}$
*   $\partial_r g^{rr} = \frac{2M}{r^2}$
*   $\partial_r g^{\theta\theta} = -\frac{2}{r^3}$
*   $\partial_r g^{\phi\phi} = -\frac{2}{r^3 \sin^2\theta}$
*   $\partial_\theta g^{\phi\phi} = -\frac{2 \cot\theta}{r^2 \sin^2\theta}$

### Step-Size Scaling

The integrator uses a proximity-based scaling for the step size $\Delta\lambda$. As a ray approaches a critical radius (the horizon at $r=2M$ or the photon sphere at $r=3M$), the step size is reduced to maintain numerical precision. Note that this is a heuristic scaling and not a formal error-estimated adaptive step size.

### Constraint Stabilization

To counteract numerical drift in $|H|$, a stabilization pass is applied after each integration step. The momentum components $p_r$ and $p_\theta$ are rescaled to satisfy the constraint $H=0$ while strictly maintaining the constants of motion $E$ and $L$.

## 4. Conserved Quantities

*   **Energy ($E$)**: $E = -p_t$ is conserved due to time translation symmetry.
*   **Angular Momentum ($L$)**: $L = p_\phi$ is conserved due to axial symmetry.
