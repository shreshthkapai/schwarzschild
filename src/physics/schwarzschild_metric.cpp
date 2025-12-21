#include "physics/schwarzschild_metric.h"
#include <cmath>
#include <cstring>
#include <iostream>

namespace Physics {

void SchwarzschildMetric::compute_metric_covariant(const double x[4], double g[4][4]) const {
    // Zero out all components
    std::memset(g, 0, 16 * sizeof(double));
    
    // Extract coordinates
    const double r = x[R];
    const double theta = x[THETA];
    
    // Compute metric function f(r) = 1 - 2M/r
    const double f_r = f(r);
    const double sin_theta = std::sin(theta);
    
    // Schwarzschild metric in coordinates (t, r, θ, φ)
    // ds² = -f(r)dt² + f(r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
    // Signature: (-,+,+,+)
    
    g[T][T] = -f_r;                                    // g_tt
    g[R][R] = 1.0 / f_r;                               // g_rr
    g[THETA][THETA] = r * r;                           // g_θθ
    g[PHI][PHI] = r * r * sin_theta * sin_theta;       // g_φφ
}

void SchwarzschildMetric::compute_metric_contravariant(const double x[4], double g[4][4]) const {
    // Zero out all components
    std::memset(g, 0, 16 * sizeof(double));
    
    // Extract coordinates
    const double r = x[R];
    const double theta = x[THETA];
    
    // Compute metric function
    const double f_r = f(r);
    const double sin_theta = std::sin(theta);
    
    // Contravariant metric (inverse of covariant)
    // For diagonal metric: g^μν = 1/g_μν
    
    g[T][T] = -1.0 / f_r;                              // g^tt
    g[R][R] = f_r;                                     // g^rr
    g[THETA][THETA] = 1.0 / (r * r);                   // g^θθ
    g[PHI][PHI] = 1.0 / (r * r * sin_theta * sin_theta);  // g^φφ
}

bool SchwarzschildMetric::is_valid_position(const double x[4]) const {
    const double r = x[R];
    
    // Position is valid if outside the singularity threshold
    // We use 2.1M as threshold (slightly outside event horizon at 2M)
    return r > SINGULARITY_THRESHOLD;
}

bool SchwarzschildMetric::validate_symmetry(const double g[4][4]) const {
    // Check that metric is symmetric: g_μν = g_νμ
    const double tol = 1e-10;
    
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            if (std::abs(g[mu][nu] - g[nu][mu]) > tol) {
                std::cerr << "Symmetry violation at (" << mu << "," << nu << "): "
                          << g[mu][nu] << " != " << g[nu][mu] << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool SchwarzschildMetric::validate_signature(const double g[4][4]) const {
    // Check signature (-,+,+,+)
    // For diagonal metric, just check signs of diagonal elements
    
    if (g[T][T] >= 0) {
        std::cerr << "Signature error: g_tt should be negative, got " << g[T][T] << std::endl;
        return false;
    }
    
    if (g[R][R] <= 0) {
        std::cerr << "Signature error: g_rr should be positive, got " << g[R][R] << std::endl;
        return false;
    }
    
    if (g[THETA][THETA] <= 0) {
        std::cerr << "Signature error: g_θθ should be positive, got " << g[THETA][THETA] << std::endl;
        return false;
    }
    
    if (g[PHI][PHI] <= 0) {
        std::cerr << "Signature error: g_φφ should be positive, got " << g[PHI][PHI] << std::endl;
        return false;
    }
    
    return true;
}

} // namespace Physics