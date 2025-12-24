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
    
    const double f_r = f(r);
    const double sin_theta = std::sin(theta);
    
    // Schwarzschild metric in coordinates (t, r, θ, φ)
    // ds² = -f(r)dt² + f(r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
    // Signature: (-,+,+,+)
    
    g[T][T] = -f_r;
    g[R][R] = 1.0 / f_r;
    g[THETA][THETA] = r * r;
    g[PHI][PHI] = r * r * sin_theta * sin_theta;
}

void SchwarzschildMetric::compute_metric_contravariant(const double x[4], double g[4][4]) const {
    // Zero out all components
    std::memset(g, 0, 16 * sizeof(double));
    
    // Extract coordinates
    const double r = x[R];
    const double theta = x[THETA];
    
    const double f_r = f(r);
    
    // Pole regularization
    const double min_sin = 1e-6;
    double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < min_sin) {
        sin_theta = (sin_theta >= 0) ? min_sin : -min_sin;
    }
    
    
    g[T][T] = -1.0 / f_r;
    g[R][R] = f_r;
    g[THETA][THETA] = 1.0 / (r * r);
    g[PHI][PHI] = 1.0 / (r * r * sin_theta * sin_theta);
}

bool SchwarzschildMetric::is_valid_position(const double x[4]) const {
    const double r = x[R];
    
    return r > SINGULARITY_THRESHOLD;
}

bool SchwarzschildMetric::validate_symmetry(const double g[4][4]) const {
    // Symmetry check: g_μν = g_νμ
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
    // Signature check (-,+,+,+)
    
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