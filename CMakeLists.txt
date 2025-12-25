cmake_minimum_required(VERSION 3.15)
project(schwarzschild_geodesics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/src)

# Source files
set(SOURCES
    src/app/main.cpp
    src/app/controls.cpp
    src/physics/schwarzschild_metric.cpp
    src/physics/hamiltonian.cpp
    src/numerics/integrator.cpp
    src/rays/ray_initializer.cpp
    src/rays/ray_bundle.cpp
    src/render/geometry.cpp
    src/render/renderer.cpp
    src/render/camera.cpp
)

# Executable
add_executable(schwarzschild ${SOURCES})

# Emscripten-specific settings
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".js")
    
    set_target_properties(schwarzschild PROPERTIES
        LINK_FLAGS "\
            -s WASM=1 \
            -s USE_WEBGL2=1 \
            -s FULL_ES3=1 \
            -s ALLOW_MEMORY_GROWTH=1 \
            -s NO_EXIT_RUNTIME=1 \
            -s INVOKE_RUN=1 \
            -s 'EXPORTED_RUNTIME_METHODS=[\"callMain\",\"ccall\",\"cwrap\"]' \
            -s ASSERTIONS=1 \
            --bind"
    )
endif()