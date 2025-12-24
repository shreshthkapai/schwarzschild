#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "physics/schwarzschild_metric.h"

namespace Physics {

// Null geodesic system
// Hamiltonian H = (1/2) g^μν p_μ p_ν

class Hamiltonian {
public:
    Hamiltonian(const SchwarzschildMetric* metric);
    
    // Hamiltonian value
    double compute_hamiltonian(const double x[4], const double p[4]) const;
    
    // dx^μ/dλ = ∂H/∂p_μ
    void compute_position_derivatives(const double x[4], const double p[4], double dx_dlambda[4]) const;
    
    // dp_μ/dλ = -∂H/∂x^μ
    void compute_momentum_derivatives(const double x[4], const double p[4], double dp_dlambda[4]) const;
    
    // RHS for integration
    void compute_rhs(const double x[4], const double p[4], 
                     double dx_dlambda[4], double dp_dlambda[4]) const;
    
    // Conserved quantities
    double compute_energy(const double x[4], const double p[4]) const;
    double compute_angular_momentum(const double x[4], const double p[4]) const;
    
private:
    const SchwarzschildMetric* metric_;
    
    // ∂g^μν/∂x^ρ
    void compute_metric_derivative(const double x[4], int rho, double dg[4][4]) const;
};

} // namespace Physics

#endif // HAMILTONIAN_H