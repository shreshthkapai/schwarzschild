#include "physics/hamiltonian.h"
#include <cmath>

namespace Physics {

Hamiltonian::Hamiltonian(const SchwarzschildMetric* metric) 
    : metric_(metric) {}

double Hamiltonian::compute_H(const double x[4], const double p[4]) const {
    // H = (1/2) g^μν p_μ p_ν
    
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    double H = 0.0;
    
    // For Schwarzschild (diagonal metric), only diagonal terms contribute
    for (int mu = 0; mu < 4; ++mu) {
        H += 0.5 * g_up[mu][mu] * p[mu] * p[mu];
    }
    
    return H;
}

void Hamiltonian::compute_position_derivatives(const double x[4], const double p[4], 
                                                double dx_dlambda[4]) const {
    // dx^μ/dλ = ∂H/∂p_μ = g^μν p_ν
    
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    // For diagonal metric: dx^μ/dλ = g^μμ p_μ (no sum)
    for (int mu = 0; mu < 4; ++mu) {
        dx_dlambda[mu] = g_up[mu][mu] * p[mu];
    }
}

void Hamiltonian::compute_momentum_derivatives(const double x[4], const double p[4], 
                                                double dp_dlambda[4]) const {
    // dp_μ/dλ = -∂H/∂x^μ = -(1/2) (∂g^αβ/∂x^μ) p_α p_β
    
    for (int mu = 0; mu < 4; ++mu) {
        dp_dlambda[mu] = 0.0;
        
        // Compute ∂g^αβ/∂x^μ
        double dg[4][4];
        compute_metric_derivative(x, mu, dg);
        
        // dp_μ/dλ = -(1/2) dg^αβ/dx^μ p_α p_β
        for (int alpha = 0; alpha < 4; ++alpha) {
            dp_dlambda[mu] -= 0.5 * dg[alpha][alpha] * p[alpha] * p[alpha];
        }
    }
}

void Hamiltonian::compute_rhs(const double x[4], const double p[4], 
                               double dx_dlambda[4], double dp_dlambda[4]) const {
    compute_position_derivatives(x, p, dx_dlambda);
    compute_momentum_derivatives(x, p, dp_dlambda);
}

double Hamiltonian::energy(const double x[4], const double p[4]) const {
    // Energy E = -p_t (conserved due to time translation symmetry)
    return -p[T];
}

double Hamiltonian::angular_momentum(const double x[4], const double p[4]) const {
    // Angular momentum L = p_φ (conserved due to axial symmetry)
    return p[PHI];
}

void Hamiltonian::compute_metric_derivative(const double x[4], int rho, double dg[4][4]) const {
    // Compute ∂g^μν/∂x^ρ numerically using finite differences
    // For Schwarzschild, we can do this analytically, but finite diff is more general
    
    const double h = 1e-8;  // Small step for numerical derivative
    
    double x_plus[4], x_minus[4];
    for (int i = 0; i < 4; ++i) {
        x_plus[i] = x[i];
        x_minus[i] = x[i];
    }
    
    x_plus[rho] += h;
    x_minus[rho] -= h;
    
    double g_plus[4][4], g_minus[4][4];
    metric_->compute_metric_contravariant(x_plus, g_plus);
    metric_->compute_metric_contravariant(x_minus, g_minus);
    
    // Central difference: dg/dx = (g(x+h) - g(x-h)) / (2h)
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            dg[mu][nu] = (g_plus[mu][nu] - g_minus[mu][nu]) / (2.0 * h);
        }
    }
}

} // namespace Physics