#ifndef SCHWARZSCHILD_METRIC_H
#define SCHWARZSCHILD_METRIC_H

#include "physics/metric.h"
#include "physics/constants.h"

namespace Physics {

class SchwarzschildMetric : public Metric {
public:
    SchwarzschildMetric() = default;
    
    // Metric tensor implementations
    void compute_metric_covariant(const double x[4], double g[4][4]) const override;
    void compute_metric_contravariant(const double x[4], double g[4][4]) const override;
    bool is_valid_position(const double x[4]) const override;
    
    // Validation functions
    bool validate_symmetry(const double g[4][4]) const;
    bool validate_signature(const double g[4][4]) const;
    
private:
    // Helper: compute f(r) = 1 - 2M/r (appears everywhere in Schwarzschild)
    inline double f(double r) const {
        return 1.0 - 2.0 * M / r;
    }
};

} // namespace Physics

#endif // SCHWARZSCHILD_METRIC_H