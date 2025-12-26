# System Architecture

The Schwarzschild Geodesic Visualization is a C++ application targeting WebAssembly/WebGL2 for real-time simulation of light trajectories in curved spacetime.

## Component Overview

### 1. Application Layer (`src/app`)
*   **Main Loop**: `main.cpp` manages the Emscripten main loop and coordinates between physics and rendering.
*   **Controls**: `controls.cpp` implements mouse interaction and keyboard shortcuts.
*   **Web Interface**: `index.html` provides the primary control panel for simulation parameters.

### 2. Physics Engine (`src/physics`)
*   **Metric**: `schwarzschild_metric.cpp` provides the metric tensor components and analytic derivatives.
*   **Hamiltonian**: `hamiltonian.cpp` defines the relativistic equations of motion.

### 3. Numerical Integration (`src/numerics`)
*   **RK4 Solver**: Standard 4th-order Runge-Kutta implementation.
*   **Geodesic Integrator**: Manages state evolution, constraint monitoring, and termination heuristics.

### 4. Ray Management (`src/rays`)
*   **Ray Initializer**: Calculates initial phase-space $(x, p)$ vectors satisfying the null condition $H=0$.
*   **Ray Bundle**: Orchestrates batch initialization for equatorial grids or spherical shells.

### 5. Rendering Engine (`src/render`)
*   **WebGL Renderer**: Handles GLES3 primitive dispatch and shader management.
*   **Geometry**: Generates scene meshes including the horizon, photon sphere, and accretion disk.

---

## Data Flow

The application follows a straightforward pipeline:

1. **User Input** → UI controls and keyboard/mouse events
2. **Parameter Update** → Simulation parameters adjusted
3. **Ray Generation** → Initial conditions computed from impact parameters
4. **Geodesic Integration** → RK4 solver computes trajectories
5. **Geometry Update** → Results converted to renderable vertices
6. **WebGL Rendering** → Scene drawn to canvas

All physics calculations occur in WebAssembly on the main thread, with results passed directly to the WebGL renderer for display.

## Implementation Notes

- **Integration**: The system uses a fixed-step solver with proximity-aware scaling rather than a formal adaptive error-estimated integrator.
- **Single-threaded**: All computations run on the main thread in WebAssembly for simplicity and reliability.
- **Real-time Updates**: Parameter changes trigger immediate recomputation and visualization refresh.

## Performance Considerations

The application maintains 60 FPS for typical ray counts (< 500 rays). For larger simulations:
- Adaptive stride reduces vertex count for distant geodesics
- LRU caching minimizes redundant computations
- Proximity-based step scaling optimizes integration near critical radii
