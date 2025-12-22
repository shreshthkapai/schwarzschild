#include "physics/hamiltonian.h"
#include <cmath>
#include "physics/constants.h"

namespace Physics {

Hamiltonian::Hamiltonian(const SchwarzschildMetric* metric) 
    : metric_(metric) {}

double Hamiltonian::compute_hamiltonian(const double x[4], const double p[4]) const {
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
        
        // dp_μ/dλ = -(1/2) (∂g^αβ/∂x^μ) p_α p_β
        // NOTE: This implementation assumes a diagonal metric (Schwarzschild) where g^αβ = 0 for α != β.
        // For a non-diagonal metric (e.g., Kerr), the inner loop must be a full double sum over α and β.
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

double Hamiltonian::compute_energy(const double x[4], const double p[4]) const {
    // Energy E = -p_t (conserved due to time translation symmetry)
    return -p[T];
}

double Hamiltonian::compute_angular_momentum(const double x[4], const double p[4]) const {
    // Angular momentum L = p_φ (conserved due to axial symmetry)
    return p[PHI];
}

void Hamiltonian::compute_metric_derivative(const double x[4], int rho, double dg[4][4]) const {
    // Analytical derivatives for Schwarzschild metric
    // Optimization: Removes 8+ metric evaluations per step compared to finite differences
    
    // Initialize to zero
    for(int i=0; i<4; ++i)
        for(int j=0; j<4; ++j)
            dg[i][j] = 0.0;
            
    double r = x[R];
    double theta = x[THETA];
    double M = Physics::M; // 1.0
    
    double r2 = r * r;
    double r3 = r2 * r;
    double one_minus_2M_r = 1.0 - 2.0*M/r;
    
    // Derivatives with respect to r
    if (rho == R) {
        // d(g^tt)/dr = 2M / (r^2 * (1-2M/r)^2)
        // g^tt = -1/(1-2M/r)
        // Note: The derivative is positive because g^tt increases from -inf to -1
        double g_tt_sq = 1.0 / (one_minus_2M_r * one_minus_2M_r);
        dg[T][T] = 2.0 * M * g_tt_sq / r2;
        
        // d(g^rr)/dr = 2M / r^2
        dg[R][R] = 2.0 * M / r2;
        
        // d(g^th th)/dr = -2 / r^3
        dg[THETA][THETA] = -2.0 / r3;
        
        // d(g^ph ph)/dr = -2 / (r^3 * sin^2(th))
        double sin_th = std::sin(theta);
        double sin2_th = sin_th * sin_th;
        dg[PHI][PHI] = -2.0 / (r3 * sin2_th);
    }
    // Derivatives with respect to theta
    else if (rho == THETA) {
        // Only g^ph ph depends on theta
        // d(g^ph ph)/dth = -2 * cot(th) * g^ph ph
        double sin_th = std::sin(theta);
        double cos_th = std::cos(theta);
        double sin_th_3 = sin_th * sin_th * sin_th;
        
        // g^ph ph = 1/(r^2 sin^2 th)
        // d/dth = 1/r^2 * (-2 sin^-3 * cos) = -2 cos / (r^2 sin^3)
        dg[PHI][PHI] = -2.0 * cos_th / (r2 * sin_th_3);
    }
}

} // namespace Physics
