#include "rays/ray_bundle.h"
#include "physics/constants.h"
#include <iostream>
#include <cmath>

namespace Rays {

RayBundle::RayBundle(const RayInitializer* initializer)
    : initializer_(initializer) {}

void RayBundle::generate_uniform_bundle(double observer_r,
                                       double impact_min,
                                       double impact_max,
                                       int num_rays) {
    rays_.clear();
    rays_.reserve(num_rays);
    
    // Uniform sampling in impact parameter (equatorial plane)
    for (int i = 0; i < num_rays; ++i) {
        double t = (num_rays > 1) ? double(i) / (num_rays - 1) : 0.5;
        double impact = impact_min + t * (impact_max - impact_min);
        
        RayState ray = initializer_->initialize_ray(
            observer_r,
            M_PI / 2.0,  // Equatorial
            0.0,         // Azimuthal start
            impact,
            0.0          // Ray direction
        );
        

        
        if (!std::isnan(ray.p[Physics::R])) {
            rays_.push_back(ray);
        }
    }
}

void RayBundle::generate_spherical_bundle(double observer_r,
                                          double impact_min,
                                          double impact_max,
                                          int num_theta,
                                          int num_phi,
                                          int num_impact) {
    rays_.clear();
    
    // Total rays = num_theta * num_phi * num_impact
    rays_.reserve(num_theta * num_phi * num_impact);
    
    std::cout << "Generating spherical bundle: " 
              << num_theta << " theta × " << num_phi << " phi × " 
              << num_impact << " impacts = " 
              << (num_theta * num_phi * num_impact) << " rays" << std::endl;
    
    // Sample theta from slightly above 0 to slightly below pi
    // (avoid exact poles where sin(θ) = 0)
    for (int i_theta = 0; i_theta < num_theta; ++i_theta) {
        double t_theta = (num_theta > 1) ? double(i_theta) / (num_theta - 1) : 0.5;
        double theta = 0.1 + t_theta * (M_PI - 0.2);  // Range [0.1, π-0.1]
        
        // At each latitude, sample different azimuthal angles
        for (int i_phi = 0; i_phi < num_phi; ++i_phi) {
            double t_phi = (num_phi > 1) ? double(i_phi) / num_phi : 0.0;
            double phi = t_phi * 2.0 * M_PI;  // Range [0, 2π)
            
            // For each direction, sample different impact parameters
            for (int i_impact = 0; i_impact < num_impact; ++i_impact) {
                double t_impact = (num_impact > 1) ? double(i_impact) / (num_impact - 1) : 0.5;
                double impact = impact_min + t_impact * (impact_max - impact_min);
                
                RayState ray = initializer_->initialize_ray(
                    observer_r,
                    theta,
                    phi,
                    impact,
                    0.0  // azimuthal angle for ray direction
                );
                
                if (!std::isnan(ray.p[Physics::R])) {
                    rays_.push_back(ray);
                } else {
                    // std::cerr << "Skipping invalid ray..." << std::endl;
                }
            }
        }
    }
    
    std::cout << "Generated " << rays_.size() << " rays in spherical pattern (skipped " 
              << (num_theta * num_phi * num_impact - rays_.size()) << " invalid)" << std::endl;
}

void RayBundle::generate_custom_bundle(double observer_r,
                                      const std::vector<double>& impact_params) {
    rays_.clear();
    rays_.reserve(impact_params.size());
    
    for (double impact : impact_params) {
        RayState ray = initializer_->initialize_ray(
            observer_r,
            M_PI / 2.0,
            0.0,
            impact,
            0.0
        );
        rays_.push_back(ray);
    }
}

void RayBundle::print_summary() const {
    std::cout << "Ray Bundle Summary:" << std::endl;
    std::cout << "  Total rays: " << rays_.size() << std::endl;
    
    if (rays_.empty()) return;
    
    double min_L = initializer_->get_angular_momentum(rays_[0]);
    double max_L = min_L;
    
    for (const auto& ray : rays_) {
        double L = initializer_->get_angular_momentum(ray);
        min_L = std::min(min_L, L);
        max_L = std::max(max_L, L);
    }
    
    std::cout << "  Impact parameter range: [" << min_L << ", " << max_L << "]" << std::endl;
    std::cout << "  Observer r: " << rays_[0].x[Physics::R] << std::endl;
}

bool RayBundle::verify_all_null(double tolerance) const {
    bool all_valid = true;
    
    for (size_t i = 0; i < rays_.size(); ++i) {
        if (!initializer_->verify_null_condition(rays_[i], tolerance)) {
            std::cerr << "Ray " << i << " fails null condition!" << std::endl;
            all_valid = false;
        }
    }
    
    return all_valid;
}

} // namespace Rays