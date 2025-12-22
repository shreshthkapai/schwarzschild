#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

namespace Render {

// Simple vertex structure for rendering
struct Vertex {
    float x, y, z;
    float r, g, b, a;  // Color
    
    Vertex(float x_, float y_, float z_, float r_=1.0f, float g_=1.0f, float b_=1.0f, float a_=1.0f)
        : x(x_), y(y_), z(z_), r(r_), g(g_), b(b_), a(a_) {}
};

// Generate wireframe sphere (for photon sphere)
std::vector<Vertex> generate_sphere(float radius, int segments, 
                                    float r, float g, float b, float a);

// Generate solid sphere triangles (for event horizon - filled black)
std::vector<Vertex> generate_solid_sphere(float radius, int segments,
                                          float r, float g, float b, float a);

// Generate unit sphere mesh (radius=1, centered at origin) for instancing
std::vector<Vertex> generate_unit_sphere_mesh(int segments);

// Generate accretion disk (glowing ring)
std::vector<Vertex> generate_accretion_disk(float inner_radius, float outer_radius,
                                            int segments, float thickness);

// Generate starfield background points
std::vector<Vertex> generate_starfield(int num_stars, float radius);

// Generate line strip from geodesic points
std::vector<Vertex> generate_geodesic_line(const std::vector<float>& points_x,
                                          const std::vector<float>& points_y,
                                          const std::vector<float>& points_z,
                                          const std::vector<float>& colors_r,
                                          const std::vector<float>& colors_g,
                                          const std::vector<float>& colors_b);

} // namespace Render

#endif // GEOMETRY_H