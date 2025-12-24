#include "render/camera.h"
#include <cmath>

namespace Render {

Camera::Camera() 
    : theta(M_PI / 4.0f), phi(M_PI / 4.0f), distance(30.0f),
      x(0), y(0), z(0) {
    update();
}

void Camera::update() {
    x = distance * sin(phi) * cos(theta);
    y = distance * cos(phi);
    z = distance * sin(phi) * sin(theta);
}

void Camera::rotate(float d_theta, float d_phi) {
    theta -= d_theta;
    phi -= d_phi;
    
    // Clamp polar angle
    if (phi < 0.1f) phi = 0.1f;
    if (phi > M_PI - 0.1f) phi = M_PI - 0.1f;
    
    update();
}

void Camera::zoom(float factor) {
    distance *= factor;
    if (distance < 5.0f) distance = 5.0f;
    if (distance > 100.0f) distance = 100.0f;
    update();
}

void Camera::compute_view_matrix(float view[16]) const {
    // Look-at origin
    float eye[3] = {x, y, z};
    float center[3] = {0, 0, 0};
    float up[3] = {0, 1, 0};
    
    // Forward vector
    float f[3] = {center[0]-eye[0], center[1]-eye[1], center[2]-eye[2]};
    float len = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    f[0]/=len; f[1]/=len; f[2]/=len;
    
    // Right vector (f x up)
    float s[3] = {f[1]*up[2]-f[2]*up[1], f[2]*up[0]-f[0]*up[2], f[0]*up[1]-f[1]*up[0]};
    len = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
    s[0]/=len; s[1]/=len; s[2]/=len;
    
    // True up vector
    float u[3] = {s[1]*f[2]-s[2]*f[1], s[2]*f[0]-s[0]*f[2], s[0]*f[1]-s[1]*f[0]};
    
    // View matrix (column-major)
    view[0]=s[0]; view[1]=u[0]; view[2]=-f[0]; view[3]=0;
    view[4]=s[1]; view[5]=u[1]; view[6]=-f[1]; view[7]=0;
    view[8]=s[2]; view[9]=u[2]; view[10]=-f[2]; view[11]=0;
    view[12]=-s[0]*eye[0]-s[1]*eye[1]-s[2]*eye[2];
    view[13]=-u[0]*eye[0]-u[1]*eye[1]-u[2]*eye[2];
    view[14]=f[0]*eye[0]+f[1]*eye[1]+f[2]*eye[2];
    view[15]=1;
}

void Camera::compute_proj_matrix(float proj[16], float aspect, float fov_degrees) const {
    float fov = fov_degrees * M_PI / 180.0f;
    float f = 1.0f / tan(fov / 2.0f);
    float near = 0.1f;
    float far = 1000.0f;
    
    proj[0]=f/aspect; proj[1]=0; proj[2]=0; proj[3]=0;
    proj[4]=0; proj[5]=f; proj[6]=0; proj[7]=0;
    proj[8]=0; proj[9]=0; proj[10]=(far+near)/(near-far); proj[11]=-1;
    proj[12]=0; proj[13]=0; proj[14]=(2*far*near)/(near-far); proj[15]=0;
}

} // namespace Render
