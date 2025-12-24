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
    
    // Equatorial uniform sampling
    for (int i = 0; i < num_rays; ++i) {
        double t = (num_rays > 1) ? double(i) / (num_rays - 1) : 0.5;
        double impact = impact_min + t * (impact_max - impact_min);
        
        RayState ray = initializer_->initialize_ray(
            observer_r,
            M_PI / 2.0,
            0.0,
            impact,
            0.0
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
    int total_rays = num_theta * num_phi * num_impact;
    rays_.reserve(total_rays);
    
    for (int i = 0; i < total_rays; ++i) {
        // Index decomposition
        int idx = i;
        int i_impact = idx % num_impact;
        idx /= num_impact;
        int i_phi = idx % num_phi;
        int i_theta = idx / num_phi;
        
        double t_theta = (num_theta > 1) ? double(i_theta) / (num_theta - 1) : 0.5;
        double theta = 0.1 + t_theta * (M_PI - 0.2);
        
        double t_phi = (num_phi > 1) ? double(i_phi) / num_phi : 0.0;
        double phi = t_phi * 2.0 * M_PI;
        
        double t_impact = (num_impact > 1) ? double(i_impact) / (num_impact - 1) : 0.5;
        double impact = impact_min + t_impact * (impact_max - impact_min);
        
        RayState ray = initializer_->initialize_ray(
            observer_r, theta, phi, impact, 0.0
        );
        
        if (!std::isnan(ray.p[Physics::R])) {
            rays_.push_back(ray);
        }
    }
}

// Helper for worker-side ray generation
RayState RayBundle::generate_ray_at_index(int index,
                                         double observer_r,
                                         double impact_min,
                                         double impact_max,
                                         int num_theta,
                                         int num_phi,
                                         int num_impact,
                                         bool use_spherical,
                                         int num_rays_2d) const {
    
    if (use_spherical) {
        // Index decomposition
        int idx = index;
        int i_impact = idx % num_impact;
        idx /= num_impact;
        int i_phi = idx % num_phi;
        int i_theta = idx / num_phi;
        
        double t_theta = (num_theta > 1) ? double(i_theta) / (num_theta - 1) : 0.5;
        double theta = 0.1 + t_theta * (M_PI - 0.2);
        
        double t_phi = (num_phi > 1) ? double(i_phi) / num_phi : 0.0;
        double phi = t_phi * 2.0 * M_PI;
        
        double t_impact = (num_impact > 1) ? double(i_impact) / (num_impact - 1) : 0.5;
        double impact = impact_min + t_impact * (impact_max - impact_min);
        
        return initializer_->initialize_ray(observer_r, theta, phi, impact, 0.0);
    } else {
        // 2D Equatorial
        double t = (num_rays_2d > 1) ? double(index) / (num_rays_2d - 1) : 0.5;
        double impact = impact_min + t * (impact_max - impact_min);
        
        return initializer_->initialize_ray(observer_r, M_PI / 2.0, 0.0, impact, 0.0);
    }
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
        
        if (!std::isnan(ray.p[Physics::R])) {
            rays_.push_back(ray);
        }
    }
}

void RayBundle::print_summary() const {
    std::cout << "RayBundle: " << rays_.size() << " rays." << std::endl;
}

bool RayBundle::verify_all_null(double tolerance) const {
    for (const auto& ray : rays_) {
        if (!initializer_->verify_null_condition(ray, tolerance)) {
            return false;
        }
    }
    return true;
}

} // namespace Rays
