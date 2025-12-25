# Schwarzschild Geodesics Visualization

Technical implementation of null geodesic integration in Schwarzschild spacetime. This project provides a WebGL-based visualization of photon trajectories around a static black hole.

**Live Demo**: [schwarzschild-vercel.vercel.app](https://schwarzschild-vercel.vercel.app/)

## Features

- **Geodesic Solver**: C++/WebAssembly implementation of the geodesic equation using a Hamiltonian formulation.
- **Physical Model**: Schwarzschild metric with analytical derivatives and Hamiltonian constraint stabilization.
- **Numerical Integration**: 4th-order Runge-Kutta (RK4) with proximity-based step-size scaling.
- **Visualization Suite**: WebGL 2.0 renderer with multiple analysis modes (termination, error gradient, Doppler shift).
- **Interactive Observation**: Orbital camera for scene navigation and a Web-based control panel for parameter adjustment.

## Technical Documentation

Detailed information on the project's components and physics:

- [System Architecture](docs/architecture.md)
- [Physics Engine](docs/physics_engine.md)
- [Visualization Pipeline](docs/visualization.md)
- [Physical Conventions](docs/conventions.md)

## Controls

The simulation is controlled via a combination of the web interface, keyboard, and mouse.

### UI Panel
The sidebar provides real-time control over:
- **Observer Radius**: Distance from the center of the black hole.
- **Ray Bundle**: Impact parameter range and angular resolution.
- **Geometry Toggles**: Event Horizon, Photon Sphere, Accretion Disk, and Starfield.

### Keyboard Shortcuts
| Key | Function |
| :--- | :--- |
| `H` | Toggle event horizon visibility |
| `P` | Toggle photon sphere visibility |
| `C` | Cycle color mapping mode |

### Mouse Interactions
- **Drag**: Rotate camera orbit |
- **Scroll**: Adjust zoom distance |

## Building and Running

### Prerequisites
- Emscripten SDK (emsdk)

### Compilation
Build using the provided shell script:

```bash
chmod +x build.sh
./build.sh
```

### Execution
Serve the project root locally:

```bash
python3 -m http.server 8000
```
Open `schwarzschild.html` in a WebGL2-compatible browser.

## References

1. Misner, C. W., Thorne, K. S., & Wheeler, J. A. (1973). *Gravitation*.
2. Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes*.

## License

MIT License.
