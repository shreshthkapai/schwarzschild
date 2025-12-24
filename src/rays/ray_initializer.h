#ifndef RAY_INITIALIZER_H
#define RAY_INITIALIZER_H

#include "rays/ray_state.h"
#include "physics/schwarzschild_metric.h"

namespace Rays {

// Ray initialization from observer coordinates and impact parameters
class RayInitializer {
public:
    RayInitializer(const Physics::SchwarzschildMetric* metric);
    
    // Initialize ray with energy normalization E=1 at infinity
    RayState initialize_ray(double observer_r, 
                           double observer_theta,
                           double observer_phi,
                           double impact_param,
                           double azimuthal_angle) const;
    
    // Verify null condition (H=0 constraint)
    bool verify_null_condition(const RayState& ray, double tolerance = 1e-10) const;
    
    // Get conserved quantities for this ray
    double get_energy(const RayState& ray) const;
    double get_angular_momentum(const RayState& ray) const;
    
private:
    const Physics::SchwarzschildMetric* metric_;
    
    // Compute radial momentum from null condition
    double compute_radial_momentum(const double x[4], const double p[4]) const;
};

} // namespace Rays

#endif // RAY_INITIALIZER_H