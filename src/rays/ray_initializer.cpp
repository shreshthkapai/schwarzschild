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
    
    // Conserved quantities: E and L
    // L = impact_parameter * E
    
    // Local energy at observer (E_at_inf = 1)
    double E = 1.0 / std::sqrt(1.0 - 2.0 * M / observer_r);
    double L = impact_param * E;  // Angular momentum L = b * E
    
    // Set conserved momenta
    ray.p[T] = -E;      // p_t = -E (energy)
    ray.p[PHI] = L;     // p_φ = L (angular momentum)
    // Solve for p_r from H = (1/2) g^μν p_μ p_ν = 0
    ray.p[THETA] = 0.0;
    ray.p[R] = compute_radial_momentum(ray.x, ray.p);
    
    // Ingoing direction
    if (ray.p[R] > 0) ray.p[R] = -ray.p[R];
    
    return ray;
}

double RayInitializer::compute_radial_momentum(const double x[4], const double p[4]) const {
    using namespace Physics;
    
    // Get contravariant metric
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    // Solve quadratic A(p_r)^2 + B(p_r) + C = 0 for p_r
    
    double A = g_up[R][R];
    double B = 0.0;
    double C = 0.0;
    
    // B = 2 * g^{rν} p_{ν!=r}, C = g^{μν} p_{μ!=r} p_{ν!=r}
    
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
    
    if (disc < 0) return std::numeric_limits<double>::quiet_NaN();
    
    // Radial root magnitude
    
    double root = (-B + std::sqrt(disc)) / (2.0 * A);
    return std::abs(root);
}

bool RayInitializer::verify_null_condition(const RayState& ray, double tolerance) const {
    using namespace Physics;
    
    double g_up[4][4];
    metric_->compute_metric_contravariant(ray.x, g_up);
    double H = 0.0;
    // H = (1/2) g^μν p_μ p_ν
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
