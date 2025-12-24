#include "physics/hamiltonian.h"
#include <cmath>
#include "physics/constants.h"

namespace Physics {

Hamiltonian::Hamiltonian(const SchwarzschildMetric* metric) 
    : metric_(metric) {}

double Hamiltonian::compute_hamiltonian(const double x[4], const double p[4]) const {
    
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    double H = 0.0;
    
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
             H += 0.5 * g_up[mu][nu] * p[mu] * p[nu];
        }
    }
    
    return H;
}

void Hamiltonian::compute_position_derivatives(const double x[4], const double p[4], 
                                                double dx_dlambda[4]) const {
    
    double g_up[4][4];
    metric_->compute_metric_contravariant(x, g_up);
    
    for (int mu = 0; mu < 4; ++mu) {
        dx_dlambda[mu] = g_up[mu][mu] * p[mu];
    }
}

void Hamiltonian::compute_momentum_derivatives(const double x[4], const double p[4], 
                                                double dp_dlambda[4]) const {
    for (int mu = 0; mu < 4; ++mu) {
        dp_dlambda[mu] = 0.0;
        
        double dg[4][4];
        compute_metric_derivative(x, mu, dg);
        
        for (int alpha = 0; alpha < 4; ++alpha) {
            for (int beta = 0; beta < 4; ++beta) {
                dp_dlambda[mu] -= 0.5 * dg[alpha][beta] * p[alpha] * p[beta];
            }
        }
    }
}

void Hamiltonian::compute_rhs(const double x[4], const double p[4], 
                               double dx_dlambda[4], double dp_dlambda[4]) const {
    compute_position_derivatives(x, p, dx_dlambda);
    compute_momentum_derivatives(x, p, dp_dlambda);
}

double Hamiltonian::compute_energy(const double x[4], const double p[4]) const {
    // E = -p_t
    return -p[T];
}

double Hamiltonian::compute_angular_momentum(const double x[4], const double p[4]) const {
    // L = p_φ
    return p[PHI];
}

void Hamiltonian::compute_metric_derivative(const double x[4], int rho, double dg[4][4]) const {
    // Analytical derivatives
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
        double sin_th = std::sin(theta);
        double cos_th = std::cos(theta);
        double sin_th_3 = sin_th * sin_th * sin_th;
        
        // d/dth = 1/r^2 * (-2 sin^-3 * cos) = -2 cos / (r^2 sin^3)
        dg[PHI][PHI] = -2.0 * cos_th / (r2 * sin_th_3);
    }
}

} // namespace Physics
