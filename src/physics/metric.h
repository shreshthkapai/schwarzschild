#ifndef METRIC_H
#define METRIC_H

namespace Physics {

// Spacetime metric base
class Metric {
public:
    virtual ~Metric() = default;
    
    // Covariant metric g_μν
    virtual void compute_metric_covariant(const double x[4], double g[4][4]) const = 0;
    
    // Contravariant metric g^μν
    virtual void compute_metric_contravariant(const double x[4], double g[4][4]) const = 0;
    
    // Position validity check
    virtual bool is_valid_position(const double x[4]) const = 0;
};

} // namespace Physics

#endif // METRIC_H