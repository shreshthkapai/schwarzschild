#ifndef RK4_H
#define RK4_H

namespace Numerics {

// 4th-order Runge-Kutta integrator (dy/dt = f(t, y))
class RK4 {
public:
    // Scalar step
    static void step_scalar(double t, double y, double h,
                           double (*f)(double, double),
                           double& y_new);
    
    // Vector step
    template<typename RHS>
    static void step(double t, const double* y, int n, double h,
                    RHS rhs, double* y_new);
};

template<typename RHS>
void RK4::step(double t, const double* y, int n, double h, RHS rhs, double* y_new) {
    // Stack-allocated buffers
    constexpr int MAX_DIM = 32;
    if (n > MAX_DIM) return;

    double k1[MAX_DIM];
    double k2[MAX_DIM];
    double k3[MAX_DIM];
    double k4[MAX_DIM];
    double y_temp[MAX_DIM];
    
    // Stage 1
    rhs(t, y, k1);
    
    // Stage 2
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + 0.5 * h * k1[i];
    }
    rhs(t + 0.5 * h, y_temp, k2);
    
    // Stage 3
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + 0.5 * h * k2[i];
    }
    rhs(t + 0.5 * h, y_temp, k3);
    
    // Stage 4
    for (int i = 0; i < n; ++i) {
        y_temp[i] = y[i] + h * k3[i];
    }
    rhs(t + h, y_temp, k4);
    
    // Accumulate result
    for (int i = 0; i < n; ++i) {
        y_new[i] = y[i] + (h / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

} // namespace Numerics

#endif // RK4_H