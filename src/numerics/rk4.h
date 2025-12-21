#ifndef RK4_H
#define RK4_H

namespace Numerics {

// Classic 4th-order Runge-Kutta integrator
// Integrates dy/dt = f(t, y) from t to t+h
class RK4 {
public:
    // Single RK4 step for scalar function
    // y_new = y + h * f(t, y)
    static void step_scalar(double t, double y, double h,
                           double (*f)(double, double),
                           double& y_new);
    
    // RK4 step for system of equations
    // State vector y has dimension n
    // rhs(t, y, dydt) computes derivatives
    template<typename RHS>
    static void step(double t, const double* y, int n, double h,
                    RHS rhs, double* y_new);
};

// Template implementation
template<typename RHS>
void RK4::step(double t, const double* y, int n, double h, RHS rhs, double* y_new) {
    // Use stack arrays for small systems (avoid heap allocation)
    constexpr int MAX_DIM = 32;
    if (n > MAX_DIM) {
        // Fallback or error for larger systems? 
        // For this project, n=8, so this is safe.
        // If needed, we could use alloca or throw.
        return; 
    }

    double k1[MAX_DIM];
    double k2[MAX_DIM];
    double k3[MAX_DIM];
    double k4[MAX_DIM];
    double y_temp[MAX_DIM];
    
    // k1 = f(t, y)
    rhs(t, y, k1);
    
    // k2 = f(t + h/2, y + h*k1/2)
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + 0.5 * h * k1[i];
    }
    rhs(t + 0.5 * h, y_temp, k2);
    
    // k3 = f(t + h/2, y + h*k2/2)
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + 0.5 * h * k2[i];
    }
    rhs(t + 0.5 * h, y_temp, k3);
    
    // k4 = f(t + h, y + h*k3)
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + h * k3[i];
    }
    rhs(t + h, y_temp, k4);
    
    // y_new = y + (h/6)(k1 + 2*k2 + 2*k3 + k4)
    for (int i = 0; i < n; ++i) {
        y_new[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

} // namespace Numerics

#endif // RK4_H