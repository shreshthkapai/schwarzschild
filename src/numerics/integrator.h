#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "physics/hamiltonian.h"
#include "numerics/rk4.h"
#include <vector>

namespace Numerics {

// Termination reasons
enum class TerminationReason {
    RUNNING,             // Integration active
    MAX_LAMBDA,          // Max affine parameter reached
    HORIZON_CROSSED,     // Event horizon crossed
    ESCAPED,             // Escape radius reached
    CONSTRAINT_VIOLATED, // Hamiltonian constraint violated
    INSTABILITY          // Numerical instability
};

// Trajectory point with diagnostics
struct GeodesicPoint {
    double lambda;          // Affine parameter
    double x[4];            // Position (t, r, θ, φ)
    double p[4];            // Momentum (p_t, p_r, p_θ, p_φ)
    
    double H;               // Hamiltonian
    double E;               // Energy
    double L;               // Angular momentum
    double Q;               // Carter constant
    double H_error;         // Constraint error
    double E_drift;         // Energy drift
    double L_drift;         // Angular momentum drift
};

// Geodesic trajectory
struct Geodesic {
    std::vector<GeodesicPoint> points;
    TerminationReason termination;
    
    double E_initial;
    double L_initial;
    
    double max_H_error;
    double max_E_drift;
    double max_L_drift;
    
    double E_drift_pct;        // Final % drift
    double L_drift_pct;        // Final % drift
    std::vector<std::string> warning_flags;
    
    Geodesic() : termination(TerminationReason::RUNNING),
                 E_initial(0), L_initial(0),
                 max_H_error(0), max_E_drift(0), max_L_drift(0),
                 E_drift_pct(0), L_drift_pct(0) {}
};

// Geodesic integrator
class GeodesicIntegrator {
public:
    GeodesicIntegrator(const Physics::Hamiltonian* ham);
    
    // Integrate a single geodesic
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
    int store_interval_;
    
    // Check termination conditions
    TerminationReason check_termination(const double x[4], const double p[4],
                                       double H_error, double lambda, double lambda_max,
                                       double initial_r) const;
    
    // Compute diagnostics for a point
    void compute_diagnostics(double lambda, const double x[4], const double p[4],
                           double E_initial, double L_initial,
                           GeodesicPoint& point) const;
                           
    // Initial condition verification
    bool verify_initial_conditions(const double x[4], double E, double L, std::string& error_msg) const;
    
    // Constraint stabilization
    void stabilize_constraints(double x[4], double p[4], double E_target) const;
};

} // namespace Numerics

#endif // INTEGRATOR_H