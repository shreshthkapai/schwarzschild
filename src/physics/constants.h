#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace Physics {

// ========================================
// GEOMETRIC UNITS: G = c = M = 1
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

constexpr double LAMBDA_STEP_DEFAULT = 0.01;    // Default affine parameter step
constexpr double LAMBDA_MAX_DEFAULT = 100.0;    // Default max affine parameter
constexpr double CONSTRAINT_TOLERANCE = 1e-6;   // Hamiltonian constraint |H| < tol
constexpr double SINGULARITY_THRESHOLD = 2.1;   // Stop integration if r < this

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