# Visualization and Rendering

The visualization engine uses WebGL 2.0 to render photon trajectories and scene geometry in real-time.

## 1. WebGL Pipeline

The renderer operates on a scene composed of static meshes (Spheres, Starfield) and dynamic line primitives (Geodesics).

### Coordinate Transformation
Physics calculations are performed in Schwarzschild spherical coordinates $(r, \theta, \phi)$. These are converted to Cartesian $(x, y, z)$ for rendering:

*   $x = r \sin\theta \cos\phi$
*   $y = r \cos\theta$
*   $z = r \sin\theta \sin\phi$

> [!NOTE]
> The vertical axis is mapped to $y$ (Cartesian) corresponding to the polar angle $\theta$.

## 2. Bloom Post-processing

To simulate the high intensity of light near the black hole, a bloom pass is implemented:

1.  **Brightness Extraction**: A shader isolates fragments above a certain luminance threshold.
2.  **Gaussian Blur**: The isolated "hot spots" are blurred horizontally and then vertically in a ping-pong buffer.
3.  **Composite**: The original scene is additive-blended with the blurred bright texture to create a glow effect.

## 3. Color Modes

The renderer supports multiple modes to visualize different aspects of the physics:

| Mode | Description |
| :--- | :--- |
| **Termination** | Colors rays based on their final state: Red (captured), Blue/Cyan (escaped), Yellow (unfinished). |
| **Doppler** | Visualizes relativistic redshift and kinetic Doppler. Redshifted photons appear redder; blueshifted appear violet/blue. |
| **Error** | Maps the numerical Hamiltonian error ($|H|$) to a color gradient, aiding in debugging integrator instability. |
| **Lensing Grid** | Projects a checkerboard pattern onto the celestial sphere to visualize gravitational lensing distortions. |

## 4. Scene Elements

*   **Event Horizon**: Rendered as a wireframe sphere at $r = 2M$ (black tint).
*   **Photon Sphere**: Rendered at $r = 3M$ (gold/yellow transparency).
*   **Einstein Ring**: A theoretical overlay at $b = \sqrt{27}M \approx 5.2M$.
*   **Accretion Disk**: A procedural disk starting at $r = 6M$ (ISCO) with a heat gradient (White to Red).
*   **Starfield**: A point cloud at large radius used as a reference background for lensing effects.
