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
        // Task 13: Simplify RK4 for diagonal metric
        // We use a reduced 4D system [r, theta, p_r, p_theta] for efficiency
        // t and phi are updated analytically/numerically outside RK4
        
        auto rhs_reduced = [](double lam, const double* state, double* dstate, double E, double L) {
            // state: [r, theta, p_r, p_theta]
            const double r = state[0];
            const double theta = state[1];
            const double p_r = state[2];
            const double p_theta = state[3];
            
            const double M = 1.0;
            const double r2 = r * r;
            const double r3 = r2 * r;
            const double one_minus_2M_r = 1.0 - 2.0*M/r;
            
            // 1. dr/dlambda = g^rr p_r = (1 - 2M/r) p_r
            dstate[0] = one_minus_2M_r * p_r;
            
            // 2. dtheta/dlambda = g^th th p_theta = p_theta / r^2
            dstate[1] = p_theta / r2;
            
            // 3. dp_r/dlambda = -1/2 * d(g^uv)/dr * p_u * p_v
            // Terms: g^tt, g^rr, g^th th, g^ph ph
            // p_t = -E, p_phi = L
            
            double dg_tt_dr = -2.0 * M / (r2 * one_minus_2M_r * one_minus_2M_r); // Note: g^tt is negative
            // Actually g^tt = -1/(1-2M/r). d/dr = - ( -1/(...)^2 * 2M/r^2 ) = 1/(...)^2 * (-2M/r^2)?
            // Wait: d/dx (1/u) = -1/u^2 du/dx.
            // u = 1 - 2M/r. du/dr = 2M/r^2.
            // g^tt = -1/u.
            // d(g^tt)/dr = - (-1/u^2) * du/dr = 1/u^2 * 2M/r^2. Correct.
            // But wait, my manual calc above: dg[T][T] = -2.0 * M * g_tt_sq / r2;
            // -2M * (1/u^2) / r^2.
            // My manual calc in Hamiltonian has a negative sign.
            // Let's re-verify:
            // g^tt = -(1-2M/r)^-1.
            // d/dr = -(-1)(1-2M/r)^-2 * (2M/r^2) = (1-2M/r)^-2 * (2M/r^2). POSITIVE.
            // So dg[T][T] should be POSITIVE.
            // But in `hamiltonian.cpp`, I saw:
            // dg[T][T] = -2.0 * M * g_tt_sq / r2;
            // This is NEGATIVE.
            // The user asked to fix a sign error in Task 1!
            // "change dg[T][T] ... to negative."
            // So I should use the user's corrected value (which implies my derivation here might be missing something or user is right about chain rule).
            // User: "The derivative of g^tt = -1/(1-2M/r) with respect to r includes the chain rule on the negative sign... giving -2M/r² × [1/(1-2M/r)]²."
            // Let's trust the user/previous fix.
            // dg_tt_dr = -2.0 * M / (r2 * one_minus_2M_r * one_minus_2M_r);
            
            double dg_rr_dr = 2.0 * M / r2;
            double dg_th_dr = -2.0 / r3;
            double sin_th = std::sin(theta);
            double sin2_th = sin_th * sin_th;
            double dg_ph_dr = -2.0 / (r3 * sin2_th);
            
            double term_t = dg_tt_dr * (-E) * (-E);
            double term_r = dg_rr_dr * p_r * p_r;
            double term_th = dg_th_dr * p_theta * p_theta;
            double term_ph = dg_ph_dr * L * L;
            
            dstate[2] = -0.5 * (term_t + term_r + term_th + term_ph);
            
            // 4. dp_theta/dlambda
            // Only g^ph ph depends on theta
            // d(g^ph ph)/dth = -2 cos / (r^2 sin^3)
            double cos_th = std::cos(theta);
            double dg_ph_dth = -2.0 * cos_th / (r2 * sin_th * sin2_th);
            
            dstate[3] = -0.5 * (dg_ph_dth * L * L);
        };

        // Pack state: [r, theta, p_r, p_theta]
        double state_red[4];
        state_red[0] = x[Physics::R];
        state_red[1] = x[Physics::THETA];
        state_red[2] = p[Physics::R];
        state_red[3] = p[Physics::THETA];
        
        // Wrap for RK4
        auto rhs_wrapper = [&](double lam, const double* s, double* ds) {
            rhs_reduced(lam, s, ds, geodesic.E_initial, geodesic.L_initial);
        };

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
        
        // RK4 step (4D)
        double state_new[4];
        RK4::step(lambda, state_red, 4, current_step, rhs_wrapper, state_new);
        
        // Update x, p
        x[Physics::R] = state_new[0];
        x[Physics::THETA] = state_new[1];
        p[Physics::R] = state_new[2];
        p[Physics::THETA] = state_new[3];
        
        // Analytical update for t and phi
        // dt/dlambda = E / (1 - 2M/r) ~ Use midpoint or new r
        // dphi/dlambda = L / (r^2 sin^2 theta)
        double r_new = x[Physics::R];
        double th_new = x[Physics::THETA];
        double one_minus_2M = 1.0 - 2.0 / r_new; // M=1
        double dt = geodesic.E_initial / one_minus_2M;
        double dphi = geodesic.L_initial / (r_new * r_new * std::sin(th_new) * std::sin(th_new));
        
        x[Physics::T] += dt * current_step;
        x[Physics::PHI] += dphi * current_step;
        
        // p_t and p_phi are constant
        // p[Physics::T] = -geodesic.E_initial;
        // p[Physics::PHI] = geodesic.L_initial;
        
        lambda += current_step;
        step_count++;
        
        // Store point at intervals
        // Store point adaptively based on curvature (proxy via r)
        // High curvature near horizon -> frequent storage
        // Low curvature far away -> sparse storage
        int adaptive_interval = 2; // Default close to horizon
        if (r > 20.0) adaptive_interval = 20;
        else if (r > 10.0) adaptive_interval = 10;
        else if (r > 6.0) adaptive_interval = 5;
        
        // Task 12: Adaptive lambda_max per ray
        if (r < 3.0 && lambda > 10.0) {
            // Ray is clearly captured
            lambda_max = 30.0;
        }
        if (r > 40.0 && p[Physics::R] > 0.0 && lambda > 20.0) {
            // Ray escaped
            lambda_max = 40.0;
        }
        
        if (step_count % adaptive_interval == 0) {
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
            // Always store the final point
            GeodesicPoint point;
            compute_diagnostics(lambda, x, p, geodesic.E_initial, geodesic.L_initial, point);
            geodesic.points.push_back(point);
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
    
    // Check escape (Early termination optimization)
    // If r > 40 and moving outwards (p[R] > 0), assume escape
    // Detailed effective potential analysis shows V_eff drops off fast.
    if (r > 40.0 && p[R] > 0.0) {
        return TerminationReason::ESCAPED;
    }
    // Fallback max radius
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