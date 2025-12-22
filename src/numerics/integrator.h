#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "physics/hamiltonian.h"
#include "numerics/rk4.h"
#include <vector>

namespace Numerics {

// Termination reasons
enum class TerminationReason {
    RUNNING,           // Still integrating
    MAX_LAMBDA,        // Reached maximum affine parameter
    HORIZON_CROSSED,   // r < r_horizon
    ESCAPED,           // r > escape_radius
    CONSTRAINT_VIOLATED, // |H| > tolerance
    INSTABILITY        // Numerical instability detected
};

// Single point along geodesic with diagnostics
struct GeodesicPoint {
    double lambda;          // Affine parameter
    double x[4];            // Position (t, r, θ, φ)
    double p[4];            // Momentum (p_t, p_r, p_θ, p_φ)
    
    // Diagnostics
    double H;               // Hamiltonian constraint value
    double E;               // Energy
    double L;               // Angular momentum
    double Q;               // Carter constant (Task 21)
    double H_error;         // |H - 0|
    double E_drift;         // |E - E_0|
    double L_drift;         // |L - L_0|
};

// Complete geodesic trajectory with diagnostics
struct Geodesic {
    std::vector<GeodesicPoint> points;
    TerminationReason termination;
    
    // Initial values (for drift computation)
    double E_initial;
    double L_initial;
    
    // Statistics
    double max_H_error;
    double max_E_drift;
    double max_L_drift;
    
    // Validation (Task 18)
    double E_drift_pct;        // Final % drift
    double L_drift_pct;        // Final % drift
    std::vector<std::string> warning_flags; // Validation warnings
    
    Geodesic() : termination(TerminationReason::RUNNING),
                 E_initial(0), L_initial(0),
                 max_H_error(0), max_E_drift(0), max_L_drift(0),
                 E_drift_pct(0), L_drift_pct(0) {}
};

// Geodesic integrator with diagnostics
class GeodesicIntegrator {
public:
    GeodesicIntegrator(const Physics::Hamiltonian* ham);
    
    // Integrate a single geodesic
    // x0, p0: initial conditions
    // lambda_step: affine parameter step size
    // lambda_max: maximum affine parameter
    Geodesic integrate(const double x0[4], const double p0[4],
                      double lambda_step, double lambda_max);
    
    // Configuration
    void set_constraint_tolerance(double tol) { constraint_tol_ = tol; }
    void set_escape_radius(double r) { escape_radius_ = r; }
    void set_store_interval(int interval) { store_interval_ = interval; }
    
private:
    const Physics::Hamiltonian* ham_;
    
    // Termination criteria
    double constraint_tol_;
    double escape_radius_;
    int store_interval_;  // Store every N steps
    
    // Check termination conditions
    TerminationReason check_termination(const double x[4], const double p[4],
                                       double H_error, double lambda, double lambda_max) const;
    
    // Compute diagnostics for a point
    void compute_diagnostics(double lambda, const double x[4], const double p[4],
                           double E_initial, double L_initial,
                           GeodesicPoint& point) const;
                           
    // Task 19: Check effective potential
    bool verify_initial_conditions(const double x[4], double E, double L, std::string& error_msg) const;
    
    // Task 24: Constraint stabilization
    void stabilize_constraints(double x[4], double p[4], double E_target) const;
};

} // namespace Numerics

#endif // INTEGRATOR_H