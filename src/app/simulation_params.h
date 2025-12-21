#ifndef SIMULATION_PARAMS_H
#define SIMULATION_PARAMS_H

namespace App {

// Runtime-adjustable simulation parameters
struct SimulationParams {
    double observer_r = 20.0;       // Observer distance from black hole
    double impact_min = 1.5;        // Minimum impact parameter (closer to critical)
    double impact_max = 8.0;        // Maximum impact parameter
    double lambda_step = 0.05;      // Integration step size
    double lambda_max = 100.0;      // Maximum affine parameter
    
    // Spherical bundle parameters (for 3D visualization)
    int num_theta = 4;              // Rings of latitude
    int num_phi = 6;               // Rays per ring  
    int num_impact = 3;             // Impact samples per angle
    // Total rays = num_theta * num_phi * num_impact = 72 rays
    
    bool use_spherical = true;      // Use spherical (3D) or uniform (2D) bundle
    int num_rays_2d = 20;           // Number of rays for 2D mode
    
    // Adjust observer radius
    void adjust_observer(double delta) {
        observer_r += delta;
        if (observer_r < 5.0) observer_r = 5.0;
        if (observer_r > 100.0) observer_r = 100.0;
    }
    
    // Adjust number of rays (affects num_phi in spherical mode)
    void adjust_rays(int delta) {
        if (use_spherical) {
            num_phi += delta / 2;
            if (num_phi < 4) num_phi = 4;
            if (num_phi > 16) num_phi = 16;
        } else {
            num_rays_2d += delta;
            if (num_rays_2d < 5) num_rays_2d = 5;
            if (num_rays_2d > 100) num_rays_2d = 100;
        }
    }
    
    // Get total ray count
    int get_total_rays() const {
        if (use_spherical) {
            return num_theta * num_phi * num_impact;
        } else {
            return num_rays_2d;
        }
    }
    
    // Adjust impact parameter range
    void adjust_impact_min(double delta) {
        impact_min += delta;
        if (impact_min < 0.5) impact_min = 0.5;
        if (impact_min > impact_max - 0.5) impact_min = impact_max - 0.5;
    }
    
    void adjust_impact_max(double delta) {
        impact_max += delta;
        if (impact_max < impact_min + 0.5) impact_max = impact_min + 0.5;
        if (impact_max > 15.0) impact_max = 15.0;
    }
};

} // namespace App

#endif // SIMULATION_PARAMS_H
