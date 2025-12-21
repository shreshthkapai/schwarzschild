#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "physics/schwarzschild_metric.h"

namespace Physics {

// Hamiltonian geodesic system for null geodesics in Schwarzschild spacetime
// Phase space: (x^μ, p_μ) where μ ∈ {t, r, θ, φ}
// Hamiltonian constraint: H = (1/2) g^μν p_μ p_ν = 0 for null geodesics

class Hamiltonian {
public:
    Hamiltonian(const SchwarzschildMetric* metric);
    
    // Compute Hamiltonian value H = (1/2) g^μν p_μ p_ν
    // For null geodesics, this should be ~0 (constraint)
    double compute_H(const double x[4], const double p[4]) const;
    
    // Hamilton's equations: dx^μ/dλ = ∂H/∂p_μ
    void compute_position_derivatives(const double x[4], const double p[4], double dx_dlambda[4]) const;
    
    // Hamilton's equations: dp_μ/dλ = -∂H/∂x^μ
    void compute_momentum_derivatives(const double x[4], const double p[4], double dp_dlambda[4]) const;
    
    // Combined RHS for integration: computes both dx/dλ and dp/dλ
    void compute_rhs(const double x[4], const double p[4], 
                     double dx_dlambda[4], double dp_dlambda[4]) const;
    
    // Constants of motion (for validation)
    double energy(const double x[4], const double p[4]) const;           // E = -p_t
    double angular_momentum(const double x[4], const double p[4]) const; // L = p_φ
    
private:
    const SchwarzschildMetric* metric_;
    
    // Helper: compute ∂g^μν/∂x^ρ for momentum derivatives
    void compute_metric_derivative(const double x[4], int rho, double dg[4][4]) const;
};

} // namespace Physics

#endif // HAMILTONIAN_H