#include <emscripten/bind.h>
#include <vector>
#include <cmath>
#include "physics/schwarzschild_metric.h"
#include "physics/hamiltonian.h"
#include "numerics/integrator.h"
#include "rays/ray_initializer.h"
#include "rays/ray_bundle.h"

using namespace emscripten;

// Compute a batch of geodesics and return serialized data
// Result format: [ray_id, termination_reason, num_points, x0, y0, z0, ...]
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
    // Pre-allocate buffer
    std::vector<float> result_buffer;
    result_buffer.reserve(count * 300); 

    for (int i = 0; i < count; ++i) {
        int ray_idx = start_index + i;
        
        RayState ray = bundle.generate_ray_at_index(
            ray_idx, observer_r, impact_min, impact_max, 
            num_theta, num_phi, num_impact, use_spherical, num_rays_2d
        );
        
        if (std::isnan(ray.p[Physics::R])) {
            continue;
        }
        
        Geodesic geo = integrator.integrate(ray.x, ray.p, lambda_step, lambda_max);
        
        // Serialized header
        result_buffer.push_back((float)ray_idx);
        result_buffer.push_back((float)geo.termination);
        result_buffer.push_back((float)geo.points.size());
        
        // Serialized body (Cartesian coordinates)
        for (const auto& pt : geo.points) {
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
            
        }
    }
    
    return result_buffer;
}

EMSCRIPTEN_BINDINGS(worker_module) {
    register_vector<float>("VectorFloat");
    
    function("compute_geodesic_batch", &compute_geodesic_batch);
}
