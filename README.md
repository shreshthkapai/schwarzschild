# Schwarzschild Geodesics Visualization

Technical implementation of null geodesic integration in Schwarzschild spacetime using a Hamiltonian formulation. This project provides a real-time WebGL visualization of photon trajectories around a non-rotating black hole.

## Features

- **Geodesic Integration**: Performance-oriented C++ implementation compiled to WebAssembly.
- **Physics Implementation**: Explicit Schwarzschild metric tensors with analytical derivatives and Hamiltonian constraint stabilization.
- **Adaptive Step-Size**: Numerical integration uses 4th-order Runge-Kutta with heuristic step-size scaling.
- **Visualization Suite**: WebGL 2.0 renderer with bloom post-processing and coordinate mapping for termination, redshift, and lensing effects.
- **Interactive Observer**: Spherical orbital camera with runtime parameter adjustment and ray-firing feedback.

## Physics and Mathematics

The system solves the equations of motion for mass-less particles (photons) governed by the Hamiltonian:

$$ H = \frac{1}{2}g^{\mu\nu} p_\mu p_\nu = 0 $$

Numerical evolution is performed in Schwarzschild coordinates, with coordinate values transformed to Cartesian space for rendering. For comprehensive technical details, refer to the documentation:

- [System Architecture](docs/architecture.md)
- [Physics Engine and Metrics](docs/physics_engine.md)
- [Visualization Pipeline](docs/visualization.md)
- [Physical Conventions](docs/conventions.md)

## Controls

The simulation is controlled via keyboard and mouse input.

| Key | Function |
| :--- | :--- |
| `Up` / `Down` | Increment/Decrement observer radius |
| `Left` / `Right` | Adjust ray count in bundle |
| `[` / `]` | Adjust impact parameter range |
| `R` | Re-initialize and fire rays (Cache flush) |
| `H` | Toggle event horizon visibility |
| `P` | Toggle photon sphere visibility |
| `D` | Toggle accretion disk visibility |
| `S` | Toggle starfield visibility |
| `C` | Cycle color mapping mode |
| `Mouse Drag` | Rotate camera |
| `Scroll` | Adjust zoom level |

## Building and Deployment

### Prerequisites
- Emscripten SDK (emsdk)

### Compilation
The project can be built using the provided shell script or CMake:

```bash
# Via build script
chmod +x build.sh
./build.sh

# Via CMake
mkdir build && cd build
emcmake cmake ..
emmake make
```

### Execution
Serve the project root using any static file server:

```bash
python3 -m http.server 8000
```
Navigate to `http://localhost:8000/schwarzschild.html` in a WebGL2 compliant browser.

## References

1. Misner, C. W., Thorne, K. S., & Wheeler, J. A. (1973). *Gravitation*.
2. Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes*.
3. Luminet, J. P. (1979). *Image of a spherical black hole with any thin accretion disk*.

## License

MIT License. See LICENSE for details.
