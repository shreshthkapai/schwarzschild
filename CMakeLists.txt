cmake_minimum_required(VERSION 3.15)
project(schwarzschild_geodesics)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR})

# Source files
set(SOURCES
    app/main.cpp
    physics/schwarzschild_metric.cpp
)

# Executable
add_executable(schwarzschild ${SOURCES})

# Emscripten-specific settings
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    
    set_target_properties(schwarzschild PROPERTIES
        LINK_FLAGS "\
            -s WASM=1 \
            -s USE_WEBGL2=1 \
            -s FULL_ES3=1 \
            -s ALLOW_MEMORY_GROWTH=1 \
            -s NO_EXIT_RUNTIME=1 \
            -s ASSERTIONS=1 \
            --shell-file ${CMAKE_SOURCE_DIR}/shell.html"
    )
endif()