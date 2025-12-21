#include "numerics/integrator.h"
#include "physics/constants.h"
#include <cmath>
#include <iostream>

namespace Numerics {

GeodesicIntegrator::GeodesicIntegrator(const Physics::Hamiltonian* ham)
    : ham_(ham),
      constraint_tol_(Physics::CONSTRAINT_TOLERANCE),
      escape_radius_(100.0),
      store_interval_(1) {}

Geodesic GeodesicIntegrator::integrate(const double x0[4], const double p0[4],
                                      double lambda_step, double lambda_max) {
    Geodesic geodesic;
    
    // Store initial constants
    geodesic.E_initial = ham_->compute_energy(x0, p0);
    geodesic.L_initial = ham_->compute_angular_momentum(x0, p0);
    
    // Current state
    double x[4], p[4];
    for (int i = 0; i < 4; ++i) {
        x[i] = x0[i];
        p[i] = p0[i];
    }
    
    double lambda = 0.0;
    int step_count = 0;
    
    // Store initial point
    GeodesicPoint initial_point;
    compute_diagnostics(lambda, x, p, geodesic.E_initial, geodesic.L_initial, initial_point);
    geodesic.points.push_back(initial_point);
    
    // Integration loop
    while (lambda < lambda_max) {
        // RHS function for RK4: combines position and momentum into single state vector
        auto rhs = [this](double lam, const double* state, double* dstate) {
            const double* x_ptr = state;
            const double* p_ptr = state + 4;
            double* dx_ptr = dstate;
            double* dp_ptr = dstate + 4;
            
            this->ham_->compute_rhs(x_ptr, p_ptr, dx_ptr, dp_ptr);
        };
        
        // Pack state: [x[0..3], p[0..3]]
        double state[8], state_new[8];
        for (int i = 0; i < 4; ++i) {
            state[i] = x[i];
            state[i+4] = p[i];
        }
        
        // Adaptive step size
        // Reduce step size near strong gravity (r close to 2M)
        const double r = x[Physics::R];
        double current_step = lambda_step;
        
        if (r < 12.0) {
            // Smoothly scale step size based on distance from horizon (r=2M)
            // Range: 1.0 at r=12, down to 0.1 at r=2
            double dist = std::max(0.0, r - Physics::R_SCHWARZSCHILD);
            double scale = dist / 10.0;
            scale = std::max(0.05, std::min(1.0, scale));
            current_step *= scale;
        }
        
        // RK4 step
        RK4::step(lambda, state, 8, current_step, rhs, state_new);
        
        // Unpack
        for (int i = 0; i < 4; ++i) {
            x[i] = state_new[i];
            p[i] = state_new[i+4];
        }
        
        lambda += current_step;
        step_count++;
        
        // Store point at intervals
        if (step_count % store_interval_ == 0) {
            GeodesicPoint point;
            compute_diagnostics(lambda, x, p, geodesic.E_initial, geodesic.L_initial, point);
            
            // Update max errors
            geodesic.max_H_error = std::max(geodesic.max_H_error, point.H_error);
            geodesic.max_E_drift = std::max(geodesic.max_E_drift, point.E_drift);
            geodesic.max_L_drift = std::max(geodesic.max_L_drift, point.L_drift);
            
            geodesic.points.push_back(point);
        }
        
        // Check termination
        double H = ham_->compute_hamiltonian(x, p);
        double H_error = std::abs(H);
        geodesic.termination = check_termination(x, p, H_error, lambda, lambda_max);
        
        if (geodesic.termination != TerminationReason::RUNNING) {
            break;
        }
    }
    
    return geodesic;
}

void GeodesicIntegrator::compute_diagnostics(double lambda, const double x[4], const double p[4],
                                            double E_initial, double L_initial,
                                            GeodesicPoint& point) const {
    point.lambda = lambda;
    
    for (int i = 0; i < 4; ++i) {
        point.x[i] = x[i];
        point.p[i] = p[i];
    }
    
    point.H = ham_->compute_hamiltonian(x, p);
    point.E = ham_->compute_energy(x, p);
    point.L = ham_->compute_angular_momentum(x, p);
    
    point.H_error = std::abs(point.H);
    point.E_drift = std::abs(point.E - E_initial);
    point.L_drift = std::abs(point.L - L_initial);
}

TerminationReason GeodesicIntegrator::check_termination(const double x[4], const double p[4],
                                                       double H_error, double lambda, 
                                                       double lambda_max) const {
    using namespace Physics;
    
    const double r = x[R];
    
    // Check horizon crossing
    if (r < SINGULARITY_THRESHOLD) {
        return TerminationReason::HORIZON_CROSSED;
    }
    
    // Check escape
    if (r > escape_radius_) {
        return TerminationReason::ESCAPED;
    }
    
    // Check constraint violation
    if (H_error > constraint_tol_) {
        return TerminationReason::CONSTRAINT_VIOLATED;
    }
    
    // Check instability (NaN detection)
    if (std::isnan(r) || std::isnan(p[R])) {
        return TerminationReason::INSTABILITY;
    }
    
    // Check max lambda
    if (lambda >= lambda_max) {
        return TerminationReason::MAX_LAMBDA;
    }
    
    return TerminationReason::RUNNING;
}

} // namespace Numerics