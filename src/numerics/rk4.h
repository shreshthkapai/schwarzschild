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
    // Allocate temporary arrays
    double* k1 = new double[n];
    double* k2 = new double[n];
    double* k3 = new double[n];
    double* k4 = new double[n];
    double* y_temp = new double[n];
    
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
    
    // Clean up
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] y_temp;
}

} // namespace Numerics

#endif // RK4_H