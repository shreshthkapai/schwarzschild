#include "rays/ray_initializer.h"
#include "physics/constants.h"
#include <cmath>
#include <iostream>
#include <limits>

namespace Rays {

RayInitializer::RayInitializer(const Physics::SchwarzschildMetric* metric)
    : metric_(metric) {}

RayState RayInitializer::initialize_ray(double observer_r, 
                                       double observer_theta,
                                       double observer_phi,
                                       double impact_param,
                                       double azimuthal_angle) const {
    using namespace Physics;
    
    RayState ray;
    
    // Set initial position (observer location)
    ray.x[T] = 0.0;
    ray.x[R] = observer_r;
    ray.x[THETA] = observer_theta;
    ray.x[PHI] = observer_phi;
    
    // For distant observer with photon aimed at black hole:
    // The impact parameter b relates to angular momentum: L = b * E
    // We set E = 1 (energy normalization), so L = b
    
    // E = 1.0 / sqrt(-g_tt) for a photon with E=1 at infinity.
    // g_tt = -(1 - 2M/r). M=1.
    double E = 1.0 / std::sqrt(1.0 - 2.0 * M / observer_r);
    double L = impact_param * E;  // Angular momentum L = b * E
    
    // Set conserved momenta
    ray.p[T] = -E;      // p_t = -E (energy)
    ray.p[PHI] = L;     // p_φ = L (angular momentum)
    
    // Equatorial motion (for now - can generalize later)
    ray.p[THETA] = 0.0;
    
    // Compute p_r from null condition: H = (1/2) g^μν p_μ p_ν = 0
    ray.p[R] = compute_radial_momentum(ray.x, ray.p);
    
    // For ingoing rays from distant observer, p_r should be negative
    if (ray.p[R] > 0) {
        ray.p[R] = -ray.p[R];
    }
    
    return ray;
}

double RayInitializer::compute_radial_momentum(const double x[4], const double p[4]) const {
    using namespace Physics;
    
    // Get contravariant metric
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    // Null condition: (1/2) g^μν p_μ p_ν = 0
    // Fix for Caveat 1: Solve general quadratic equation for p_r
    // A (p_r)^2 + B (p_r) + C = 0
    
    double A = g_up[R][R];
    double B = 0.0;
    double C = 0.0;
    
    // Calculate B and C
    // B = 2 * sum_{ν != r} g^{rν} p_ν
    // C = sum_{μ != r} sum_{ν != r} g^{μν} p_μ p_ν
    
    for (int mu = 0; mu < 4; ++mu) {
        if (mu == R) continue;
        B += 2.0 * g_up[R][mu] * p[mu];
        
        for (int nu = 0; nu < 4; ++nu) {
            if (nu == R) continue;
            C += g_up[mu][nu] * p[mu] * p[nu];
        }
    }
    
    // Solve Quadratic: A*x^2 + B*x + C = 0
    // Discriminant
    double disc = B*B - 4.0*A*C;
    
    if (disc < 0) {
        // Return NaN to signal invalid initial conditions
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    // We want the magnitude, direction is handled by caller (who flips it to negative)
    // We compute one root, take absolute value, and let caller handle sign.
    
    double root = (-B + std::sqrt(disc)) / (2.0 * A);
    return std::abs(root);
}

bool RayInitializer::verify_null_condition(const RayState& ray, double tolerance) const {
    using namespace Physics;
    
    // Compute H = (1/2) g^μν p_μ p_ν
    double g_up[4][4];
    metric_->compute_metric_contravariant(ray.x, g_up);
    
    double H = 0.0;
    // Fix for Caveat 1: Full sum
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            H += 0.5 * g_up[mu][nu] * ray.p[mu] * ray.p[nu];
        }
    }
    
    return std::abs(H) < tolerance;
}

double RayInitializer::get_energy(const RayState& ray) const {
    // E = -p_t
    return -ray.p[Physics::T];
}

double RayInitializer::get_angular_momentum(const RayState& ray) const {
    // L = p_φ
    return ray.p[Physics::PHI];
}

} // namespace Rays
