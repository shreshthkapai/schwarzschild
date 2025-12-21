# Schwarzschild Geodesics Visualization

Interactive WebGL visualization of photon trajectories around a Schwarzschild black hole, computed using the Hamiltonian formulation of general relativity.

![Demo](docs/demo.png)

## 🎮 Live Demo

[**Launch Demo →**](https://yourusername.github.io/schwarzschild-geodesics/)

## Physics

### Geometric Units
We use **geometric units** where:
- G = c = M = 1

All distances are measured in units of black hole mass M.

### Schwarzschild Metric
The spacetime around a non-rotating black hole:

```
ds² = -(1 - 2M/r) dt² + (1 - 2M/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²
```

### Hamiltonian Formulation
We evolve photon trajectories using the Hamiltonian:

```
H = ½ gᵘᵛ pᵤ pᵥ = 0    (null geodesic constraint)
```

Hamilton's equations give the geodesic evolution:

```
dx^μ/dλ =  ∂H/∂pᵤ = gᵘᵛ pᵥ
dpᵤ/dλ = -∂H/∂x^μ = -½ (∂gᵅᵝ/∂x^μ) pᵅ pᵝ
```

### Constants of Motion
Conserved due to spacetime symmetries:
- **Energy**: E = -pₜ (time translation symmetry)
- **Angular momentum**: L = pᵩ (axial symmetry)

### Critical Radii
| Radius | Value | Physical Meaning |
|--------|-------|------------------|
| Event Horizon | r = 2M | Point of no return |
| Photon Sphere | r = 3M | Unstable circular photon orbit |
| ISCO | r = 6M | Innermost stable circular orbit (massive particles) |

## Controls

| Key | Action |
|-----|--------|
| `↑` / `↓` | Adjust observer distance (±2M) |
| `←` / `→` | Adjust number of rays (±5) |
| `[` / `]` | Adjust impact parameter range (±0.5) |
| `R` | Refire rays with current parameters |
| `H` | Toggle event horizon sphere |
| `P` | Toggle photon sphere |
| `C` | Cycle color mode (termination → error → solid) |
| `I` | Print current parameters to console |
| Mouse drag | Rotate camera |
| Scroll | Zoom in/out |

## Color Coding

- **Red rays**: Captured by black hole (crossed horizon)
- **Green rays**: Escaped to infinity
- **Blue rays**: Still integrating (hit λ_max)

## Building

### Prerequisites
- [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)

### Compile
```bash
source /path/to/emsdk/emsdk_env.sh
./build.sh
```

### Run
```bash
emrun schwarzschild.html
# or
python3 -m http.server 8000
# then open http://localhost:8000/schwarzschild.html
```

## Project Structure

```
src/
├── app/
│   ├── main.cpp              # Entry point, render loop
│   ├── controls.cpp/.h       # Keyboard/mouse input
│   └── simulation_params.h   # Runtime parameters
├── physics/
│   ├── constants.h           # G=c=M=1, coordinate indices
│   ├── schwarzschild_metric.cpp/.h  # Metric tensors
│   └── hamiltonian.cpp/.h    # Hamilton's equations
├── numerics/
│   ├── rk4.h                 # 4th-order Runge-Kutta
│   └── integrator.cpp/.h     # Geodesic integration + diagnostics
├── rays/
│   ├── ray_state.h           # Phase space state
│   ├── ray_initializer.cpp/.h # Null geodesic initial conditions
│   └── ray_bundle.cpp/.h     # Bundle of rays
└── render/
    ├── camera.cpp/.h         # View/projection matrices
    ├── renderer.cpp/.h       # WebGL rendering
    └── geometry.cpp/.h       # Sphere/line generation
```

## Assumptions & Limitations

### Physical Assumptions
- **Schwarzschild spacetime**: Non-rotating, uncharged black hole
- **Null geodesics only**: Massless particles (photons)
- **Equatorial plane**: Rays initialized and confined to θ = π/2
- **Geometric optics**: No wave effects, diffraction, or gravitational lensing of extended sources

### Numerical Limitations
- **RK4 integrator**: Fixed-step, 4th-order accurate. Not symplectic.
- **Constraint drift**: H ≈ 0 monitored but not enforced
- **No adaptive stepping**: May miss fine details near critical radius
- **Horizon cutoff**: Terminates at r = 2.1M (slightly outside horizon)

### Visualization Limitations
- **Schwarzschild coordinates**: Not horizon-penetrating
- **No accretion disk**: Just geodesics + spheres
- **No shadows**: Not a raytraced image of a black hole

## Future Improvements
- [ ] Kerr metric (rotating black holes)
- [ ] Adaptive step size
- [ ] Accretion disk visualization
- [ ] Shadow computation
- [ ] VR support

## References

1. Misner, Thorne, Wheeler — *Gravitation* (1973)
2. Chandrasekhar — *The Mathematical Theory of Black Holes* (1983)
3. Luminet — "Image of a Spherical Black Hole..." (1979)

## License

MIT License — See [LICENSE](LICENSE) for details.
