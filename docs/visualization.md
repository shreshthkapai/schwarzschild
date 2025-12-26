# Visualization and Rendering

The visualization engine uses WebGL 2.0 to render photon trajectories and scene geometry in real-time.

## 1. WebGL Pipeline

The renderer translates physics results into visual primitives using a standard GLES3/WebAssembly bridge.

### Coordinate Transformation
Physics calculations are performed in Schwarzschild spherical coordinates $(r, \theta, \phi)$. These are converted to Cartesian $(x, y, z)$ for rendering:

*   $x = r \sin\theta \cos\phi$
*   $y = r \cos\theta$
*   $z = r \sin\theta \sin\phi$

## 2. Rendering Layers

### Scene Geometry
The renderer manages several geometric layers to provide context for the geodesics:
- **Event Horizon**: A wireframe sphere at $r = 2M$.
- **Photon Sphere**: A semi-transparent sphere at $r = 3M$.
- **Accretion Disk**: A procedural mesh representing the equatorial plane at $r \ge 6M$.
- **Starfield**: A background point cloud for lensing reference.

### Color Mapping
Trajectory colors can be cycled via the `C` key:
- **Termination Status**: Visualizes capture (red) vs escape (cyan/blue).
- **Redshift/Doppler**: Maps frequency shifts to the visible spectrum.
- **Hamiltonian Error**: Debug mode for numerical stability monitoring.
- **Lensing Grid**: Checkerboard projection on the sky sphere.

## 3. Implementation Details

### Geometry Generation
All scene geometry is generated procedurally:
- Spheres use latitude/longitude tessellation
- The accretion disk uses radial segments with gradient coloring
- Starfield uses random spherical distribution

### Adaptive Rendering
The renderer optimizes performance through:
- **Vertex stride adaptation**: Reduces vertex count for distant geodesics
- **Batched rendering**: Separates captured, escaped, and other trajectories
- **Alpha blending**: Proper transparency for overlapping paths

### Einstein Ring Overlay
The critical impact parameter ($b = \sqrt{27}M \approx 5.196M$) is visualized as a yellow wireframe sphere, showing where photons orbit the black hole before escaping.