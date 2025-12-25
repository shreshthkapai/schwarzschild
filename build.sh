#!/bin/bash

echo "Building Schwarzschild Geodesics Visualization..."
source /home/shreshth/emsdk/emsdk_env.sh

# Build main app
echo "Building main application..."
em++ -std=c++17 -Isrc \
  src/app/main.cpp \
  src/app/controls.cpp \
  src/physics/schwarzschild_metric.cpp \
  src/physics/hamiltonian.cpp \
  src/numerics/integrator.cpp \
  src/rays/ray_initializer.cpp \
  src/rays/ray_bundle.cpp \
  src/render/geometry.cpp \
  src/render/renderer.cpp \
  src/render/camera.cpp \
  -s WASM=1 \
  -s USE_WEBGL2=1 \
  -s FULL_ES3=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s NO_EXIT_RUNTIME=1 \
  -s INVOKE_RUN=1 \
  -s 'EXPORTED_RUNTIME_METHODS=["callMain","ccall","cwrap"]' \
  --bind \
  -o schwarzschild.js

echo ""
echo "Build complete!"
echo "Output files:"
echo "  - schwarzschild.js"
echo "  - schwarzschild.wasm"
echo ""
echo "To run:"
echo "  python3 -m http.server 8000"
echo "  Open http://localhost:8000/index.html"
echo ""
echo "Controls:"
echo "  Mouse drag: Rotate camera"
echo "  Scroll: Zoom"
echo "  H: Toggle horizon"
echo "  P: Toggle photon sphere"
echo "  C: Cycle color mode"