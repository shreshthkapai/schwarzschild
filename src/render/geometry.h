#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

namespace Render {

// Render vertex
struct Vertex {
    float x, y, z;
    float r, g, b, a;  // Color
    
    Vertex(float x_, float y_, float z_, float r_=1.0f, float g_=1.0f, float b_=1.0f, float a_=1.0f)
        : x(x_), y(y_), z(z_), r(r_), g(g_), b(b_), a(a_) {}
};

// Wireframe sphere
std::vector<Vertex> generate_sphere(float radius, int segments, 
                                    float r, float g, float b, float a);

// Solid sphere
std::vector<Vertex> generate_solid_sphere(float radius, int segments,
                                          float r, float g, float b, float a);

// Unit sphere (radius=1)
std::vector<Vertex> generate_unit_sphere_mesh(int segments);

// Accretion disk
std::vector<Vertex> generate_accretion_disk(float inner_radius, float outer_radius,
                                            int segments, float thickness);

// Starfield
std::vector<Vertex> generate_starfield(int num_stars, float radius);

// Geodesic line
std::vector<Vertex> generate_geodesic_line(const std::vector<float>& points_x,
                                          const std::vector<float>& points_y,
                                          const std::vector<float>& points_z,
                                          const std::vector<float>& colors_r,
                                          const std::vector<float>& colors_g,
                                          const std::vector<float>& colors_b);

} // namespace Render

#endif // GEOMETRY_H