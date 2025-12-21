#ifndef RAY_BUNDLE_H
#define RAY_BUNDLE_H

#include "rays/ray_state.h"
#include "rays/ray_initializer.h"
#include <vector>

namespace Rays {

// Collection of rays with consistent initialization
class RayBundle {
public:
    RayBundle(const RayInitializer* initializer);
    
    // Generate bundle of rays with uniform impact parameter sampling (equatorial only)
    void generate_uniform_bundle(double observer_r,
                                double impact_min,
                                double impact_max,
                                int num_rays);
    
    // NEW: Generate spherical ray bundle - rays from many angles around the black hole
    // Creates a grid in (θ, φ) space with varying impact parameters
    // Much more realistic visualization of light bending from all directions
    void generate_spherical_bundle(double observer_r,
                                   double impact_min,
                                   double impact_max,
                                   int num_theta,      // rings of latitude
                                   int num_phi,        // rays per ring
                                   int num_impact);    // impact samples per ray angle
    
    // Generate bundle with specific impact parameters (for testing)
    void generate_custom_bundle(double observer_r,
                               const std::vector<double>& impact_params);
    
    // Access rays
    const std::vector<RayState>& get_rays() const { return rays_; }
    int size() const { return rays_.size(); }
    
    // Diagnostics
    void print_summary() const;
    bool verify_all_null(double tolerance = 1e-10) const;
    
private:
    const RayInitializer* initializer_;
    std::vector<RayState> rays_;
};

} // namespace Rays

#endif // RAY_BUNDLE_H