#include "rays/ray_initializer.h"
#include "physics/constants.h"
#include <cmath>
#include <iostream>

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
    
    double E = 1.0;  // Energy (normalized)
    double L = impact_param;  // Angular momentum = impact parameter
    
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
    // For Schwarzschild (diagonal metric):
    // (1/2)[g^tt p_t² + g^rr p_r² + g^θθ p_θ² + g^φφ p_φ²] = 0
    
    // Solve for p_r²:
    // g^rr p_r² = -(g^tt p_t² + g^θθ p_θ² + g^φφ p_φ²)
    
    double sum_other = g_up[T][T] * p[T] * p[T] +
                       g_up[THETA][THETA] * p[THETA] * p[THETA] +
                       g_up[PHI][PHI] * p[PHI] * p[PHI];
    
    double pr_squared = -sum_other / g_up[R][R];
    
    if (pr_squared < 0) {
        std::cerr << "Warning: pr_squared < 0 in compute_radial_momentum" << std::endl;
        std::cerr << "  This should not happen for valid initial conditions" << std::endl;
        return 0.0;
    }
    
    return std::sqrt(pr_squared);
}

bool RayInitializer::verify_null_condition(const RayState& ray, double tolerance) const {
    using namespace Physics;
    
    // Compute H = (1/2) g^μν p_μ p_ν
    double g_up[4][4];
    metric_->compute_metric_contravariant(ray.x, g_up);
    
    double H = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        H += 0.5 * g_up[mu][mu] * ray.p[mu] * ray.p[mu];
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