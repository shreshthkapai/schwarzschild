#include <emscripten/bind.h>
#include <vector>
#include <cmath>
#include "physics/schwarzschild_metric.h"
#include "physics/hamiltonian.h"
#include "numerics/integrator.h"
#include "rays/ray_initializer.h"
#include "rays/ray_bundle.h"

using namespace emscripten;

// Compute a batch of rays and return serialized data
// Format: [ray_id, termination_code, num_points, x0, y0, z0, x1, y1, z1, ..., ray_id, ...]
std::vector<float> compute_geodesic_batch(
    int start_index,
    int count,
    double observer_r,
    double impact_min,
    double impact_max,
    int num_theta,
    int num_phi,
    int num_impact,
    bool use_spherical,
    int num_rays_2d,
    double lambda_step,
    double lambda_max
) {
    using namespace Physics;
    using namespace Rays;
    using namespace Numerics;

    SchwarzschildMetric metric;
    Hamiltonian ham(&metric);
    RayInitializer ray_init(&metric);
    RayBundle bundle(&ray_init);
    GeodesicIntegrator integrator(&ham);
    
    // Configure integrator
    integrator.set_store_interval(Physics::STORE_INTERVAL);
    
    std::vector<float> result_buffer;
    // Pre-allocate to avoid frequent reallocations (estimate: count * 100 points * 3 coords)
    result_buffer.reserve(count * 300); 

    for (int i = 0; i < count; ++i) {
        int ray_idx = start_index + i;
        
        RayState ray = bundle.generate_ray_at_index(
            ray_idx, observer_r, impact_min, impact_max, 
            num_theta, num_phi, num_impact, use_spherical, num_rays_2d
        );
        
        if (std::isnan(ray.p[Physics::R])) {
            // Invalid ray (e.g. inside horizon or error), skip
            continue;
        }
        
        Geodesic geo = integrator.integrate(ray.x, ray.p, lambda_step, lambda_max);
        
        // Serialize
        // Header: [ray_id, termination_reason, num_points]
        result_buffer.push_back((float)ray_idx);
        result_buffer.push_back((float)geo.termination);
        result_buffer.push_back((float)geo.points.size());
        
        // Body: [x, y, z] for each point
        for (const auto& pt : geo.points) {
            // Convert to Cartesian for rendering
            // x = r * sin(theta) * cos(phi)
            // y = r * sin(theta) * sin(phi)
            // z = r * cos(theta)
            
            // Apply observer rotation if needed? 
            // For now, raw Schwarzschild to Cartesian
            float r = (float)pt.x[Physics::R];
            float th = (float)pt.x[Physics::THETA];
            float ph = (float)pt.x[Physics::PHI];
            
            float sin_th = std::sin(th);
            float x = r * sin_th * std::cos(ph);
            float y = r * sin_th * std::sin(ph);
            float z = r * std::cos(th);
            
            result_buffer.push_back(x);
            result_buffer.push_back(y);
            result_buffer.push_back(z);
            
            // Add color/error info? 
            // We can infer color from termination in JS, or add explicit color here
            // Let's keep it minimal: position only. 
            // Color is decided by termination code in JS.
        }
    }
    
    return result_buffer;
}

EMSCRIPTEN_BINDINGS(worker_module) {
    register_vector<float>("VectorFloat");
    
    function("compute_geodesic_batch", &compute_geodesic_batch);
}
    function("compute_geodesic_batch", &compute_geodesic_batch);
}
