#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>

namespace Render {

class Camera {
public:
    Camera();
    
    // Spherical coordinates
    float theta;    // Azimuthal angle
    float phi;      // Polar angle  
    float distance; // Radial distance
    
    // Cartesian coordinates
    float x, y, z;
    
    // Update position
    void update();
    
    // Rotate camera
    void rotate(float d_theta, float d_phi);
    
    // Zoom camera
    void zoom(float factor);
    
    // View matrix
    void compute_view_matrix(float view[16]) const;
    
    // Projection matrix
    void compute_proj_matrix(float proj[16], float aspect, float fov_degrees = 60.0f) const;
};

} // namespace Render

#endif // CAMERA_H
