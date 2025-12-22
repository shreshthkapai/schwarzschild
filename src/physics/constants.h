#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace Physics {

// ========================================
// GEOMETRIC UNITS: G = c = 1 (M is scalable, currently 1.0)
// ========================================

constexpr double G = 1.0;  // Gravitational constant
constexpr double c = 1.0;  // Speed of light
constexpr double M = 1.0;  // Black hole mass

// ========================================
// CRITICAL RADII (in units of M)
// ========================================

constexpr double R_SCHWARZSCHILD = 2.0 * M;  // Event horizon
constexpr double R_PHOTON_SPHERE = 3.0 * M;  // Photon sphere (unstable circular orbit)
constexpr double R_ISCO = 6.0 * M;           // Innermost stable circular orbit (timelike)

// ========================================
// NUMERICAL PARAMETERS
// ========================================

constexpr double LAMBDA_STEP_DEFAULT = 0.05;    // Default affine parameter step (optimized)
constexpr double LAMBDA_MAX_DEFAULT = 100.0;    // Default max affine parameter
constexpr double CONSTRAINT_TOLERANCE = 1e-6;   // Hamiltonian constraint |H| < tol
constexpr double SINGULARITY_THRESHOLD = 2.1;   // Stop integration if r < this
constexpr int STORE_INTERVAL = 15;              // Store every Nth point for rendering

// ========================================
// INITIAL CONDITIONS
// ========================================

constexpr double OBSERVER_RADIUS_DEFAULT = 20.0;  // Default observer position
constexpr double OBSERVER_THETA = M_PI / 2.0;     // Equatorial plane
constexpr double OBSERVER_PHI = 0.0;

// Impact parameter range for ray bundle
constexpr double IMPACT_PARAM_MIN = 0.0;
constexpr double IMPACT_PARAM_MAX = 10.0;
constexpr int NUM_RAYS_DEFAULT = 50;

// ========================================
// VISUAL PARAMETERS
// ========================================
constexpr float DISK_INNER_R = 4.0f;
constexpr float DISK_OUTER_R = 12.0f;
constexpr int DISK_SEGMENTS = 60;
constexpr int SPHERE_SEGMENTS = 40;
constexpr int PHOTON_SPHERE_SEGMENTS = 30;
constexpr int STAR_COUNT = 500;
constexpr float STAR_DIST = 800.0f;

// ========================================
// INTEGRATION PARAMETERS
// ========================================
constexpr double ADAPTIVE_STEP_RADIUS = 12.0;  // Start adaptive stepping below this radius
constexpr double ADAPTIVE_STEP_DISTANCE = 10.0;  // Distance scaling factor
constexpr double ADAPTIVE_STEP_MIN = 0.05;      // Minimum step scale
constexpr double ADAPTIVE_STEP_MAX = 1.0;       // Maximum step scale
constexpr int CONSTRAINT_STABILIZE_INTERVAL = 10; // Stabilize constraints every N steps
constexpr double CONSTRAINT_TOL = 1e-10;        // Constraint tolerance for stabilization

// ========================================
// RENDERING PARAMETERS
// ========================================
constexpr float COLOR_RED_R = 1.0f;              // Red color component for captured rays
constexpr float COLOR_RED_G = 0.266f;            // Green component (#ff4444)
constexpr float COLOR_RED_B = 0.266f;             // Blue component
constexpr float BLOOM_THRESHOLD = 0.8f;          // Brightness threshold for bloom
constexpr int BLOOM_BLUR_PASSES = 10;            // Number of blur passes

// ========================================
// COORDINATE CONVENTIONS
// ========================================

// Indices for phase space arrays
enum CoordIndex {
    T = 0,    // Schwarzschild time
    R = 1,    // Radial coordinate
    THETA = 2,// Polar angle
    PHI = 3   // Azimuthal angle
};

// Metric signature: (-,+,+,+)
constexpr int METRIC_SIGNATURE[4] = {-1, 1, 1, 1};

} // namespace Physics

#endif // CONSTANTS_H