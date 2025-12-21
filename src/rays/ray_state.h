#ifndef RAY_STATE_H
#define RAY_STATE_H

namespace Rays {

// Single ray state (position + momentum in phase space)
struct RayState {
    double x[4];  // Position: (t, r, θ, φ)
    double p[4];  // Momentum: (p_t, p_r, p_θ, p_φ)
    
    RayState() {
        for (int i = 0; i < 4; ++i) {
            x[i] = 0.0;
            p[i] = 0.0;
        }
    }
};

} // namespace Rays

#endif // RAY_STATE_H