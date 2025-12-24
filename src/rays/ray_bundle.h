#ifndef RAY_BUNDLE_H
#define RAY_BUNDLE_H

#include "rays/ray_state.h"
#include "rays/ray_initializer.h"
#include <vector>

namespace Rays {

// Set of rays with common initialization
class RayBundle {
public:
    RayBundle(const RayInitializer* initializer);
    
    // Generate 2D uniform bundle (equatorial)
    void generate_uniform_bundle(double observer_r,
                                double impact_min,
                                double impact_max,
                                int num_rays);
    
    // Generate 3D spherical bundle (θ, φ, impact grid)
    void generate_spherical_bundle(double observer_r,
                                   double impact_min,
                                   double impact_max,
                                   int num_theta,
                                   int num_phi,
                                   int num_impact);
    
    // Stateless ray generation at index
    RayState generate_ray_at_index(int index,
                                  double observer_r,
                                  double impact_min,
                                  double impact_max,
                                  int num_theta,
                                  int num_phi,
                                  int num_impact,
                                  bool use_spherical,
                                  int num_rays_2d) const;
                                  
    // Generate bundle from explicit impact parameters
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