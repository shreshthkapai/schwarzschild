#!/bin/bash

echo "Building with direct em++ compilation..."
source /home/shreshth/emsdk/emsdk_env.sh

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
  --bind \
  --shell-file shell.html \
  -o schwarzschild.html

echo "Built: schwarzschild.html"
echo "Run: emrun schwarzschild.html"
echo ""
echo "Controls:"
echo "  Mouse drag: Rotate camera"
echo "  Scroll: Zoom"
echo "  Up/Down: Observer distance"
echo "  Left/Right: Number of rays"
echo "  [ / ]: Impact parameter range"
echo "  R: Refire rays"
echo "  H: Toggle horizon"
echo "  P: Toggle photon sphere"
echo "  C: Cycle color mode"