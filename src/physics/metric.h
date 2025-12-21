#ifndef METRIC_H
#define METRIC_H

namespace Physics {

// Abstract base class for spacetime metrics
class Metric {
public:
    virtual ~Metric() = default;
    
    // Compute covariant metric components g_μν at position x
    virtual void compute_metric_covariant(const double x[4], double g[4][4]) const = 0;
    
    // Compute contravariant metric components g^μν at position x
    virtual void compute_metric_contravariant(const double x[4], double g[4][4]) const = 0;
    
    // Check if position is in valid region (outside singularity)
    virtual bool is_valid_position(const double x[4]) const = 0;
};

} // namespace Physics

#endif // METRIC_H