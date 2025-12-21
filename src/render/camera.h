#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>

namespace Render {

class Camera {
public:
    Camera();
    
    // Spherical coordinates
    float theta;    // Horizontal angle
    float phi;      // Vertical angle  
    float distance; // Distance from origin
    
    // Cartesian position (computed)
    float x, y, z;
    
    // Update Cartesian position from spherical
    void update();
    
    // Rotate camera by delta angles
    void rotate(float d_theta, float d_phi);
    
    // Zoom camera by factor
    void zoom(float factor);
    
    // Compute view matrix (look at origin)
    void compute_view_matrix(float view[16]) const;
    
    // Compute projection matrix
    void compute_proj_matrix(float proj[16], float aspect, float fov_degrees = 60.0f) const;
};

} // namespace Render

#endif // CAMERA_H
