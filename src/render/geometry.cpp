#include "render/geometry.h"
#include <cmath>
#include <cstdlib>

namespace Render {

std::vector<Vertex> generate_sphere(float radius, int segments,
                                    float r, float g, float b, float a) {
    std::vector<Vertex> vertices;
    const float pi = 3.14159265359f;
    
    // Latitude lines (for wireframe)
    for (int lat = 0; lat < segments; ++lat) {
        float theta1 = pi * float(lat) / float(segments);
        
        for (int lon = 0; lon < segments; ++lon) {
            float phi1 = 2.0f * pi * float(lon) / float(segments);
            // Just generate points for lines
            float x = radius * std::sin(theta1) * std::cos(phi1);
            float y = radius * std::cos(theta1);
            float z = radius * std::sin(theta1) * std::sin(phi1);
            vertices.push_back(Vertex(x, y, z, r, g, b, a));
            
            x = radius * std::sin(theta1) * std::cos(phi1 + 0.1f); 
            y = radius * std::cos(theta1);
            z = radius * std::sin(theta1) * std::sin(phi1 + 0.1f);
            vertices.push_back(Vertex(x, y, z, r, g, b, a));
        }
    }
    return vertices;
}

std::vector<Vertex> generate_solid_sphere(float radius, int segments,
                                          float r, float g, float b, float a) {
    std::vector<Vertex> vertices;
    const float pi = 3.14159265359f;
    
    for (int lat = 0; lat < segments; ++lat) {
        float theta1 = pi * float(lat) / float(segments);
        float theta2 = pi * float(lat+1) / float(segments);
        
        for (int lon = 0; lon < segments; ++lon) {
            float phi1 = 2.0f * pi * float(lon) / float(segments);
            float phi2 = 2.0f * pi * float(lon+1) / float(segments);
            
            // 4 corners of a quad
            float x1 = radius * std::sin(theta1) * std::cos(phi1);
            float y1 = radius * std::cos(theta1);
            float z1 = radius * std::sin(theta1) * std::sin(phi1);
            
            float x2 = radius * std::sin(theta1) * std::cos(phi2);
            float y2 = radius * std::cos(theta1);
            float z2 = radius * std::sin(theta1) * std::sin(phi2);
            
            float x3 = radius * std::sin(theta2) * std::cos(phi2);
            float y3 = radius * std::cos(theta2);
            float z3 = radius * std::sin(theta2) * std::sin(phi2);
            
            float x4 = radius * std::sin(theta2) * std::cos(phi1);
            float y4 = radius * std::cos(theta2);
            float z4 = radius * std::sin(theta2) * std::sin(phi1);
            
            // 2 triangles
            vertices.push_back(Vertex(x1, y1, z1, r, g, b, a));
            vertices.push_back(Vertex(x2, y2, z2, r, g, b, a));
            vertices.push_back(Vertex(x3, y3, z3, r, g, b, a));
            
            vertices.push_back(Vertex(x1, y1, z1, r, g, b, a));
            vertices.push_back(Vertex(x3, y3, z3, r, g, b, a));
            vertices.push_back(Vertex(x4, y4, z4, r, g, b, a));
        }
    }
    return vertices;
}

std::vector<Vertex> generate_unit_sphere_mesh(int segments) {
    // Return wireframe grid for instancing
    return generate_sphere(1.0f, segments, 1.0f, 1.0f, 1.0f, 1.0f);
}

std::vector<Vertex> generate_accretion_disk(float inner_radius, float outer_radius,
                                            int segments, float thickness) {
    std::vector<Vertex> vertices;
    const float pi = 3.14159265359f;
    
    for (int i = 0; i < segments; ++i) {
        float phi1 = 2.0f * pi * float(i) / float(segments);
        float phi2 = 2.0f * pi * float(i+1) / float(segments);
        
        float r_in = inner_radius;
        float r_out = outer_radius;
        
        // Inner points
        float x1 = r_in * std::cos(phi1);
        float z1 = r_in * std::sin(phi1);
        
        float x2 = r_in * std::cos(phi2);
        float z2 = r_in * std::sin(phi2);
        
        // Outer points
        float x3 = r_out * std::cos(phi2);
        float z3 = r_out * std::sin(phi2);
        
        float x4 = r_out * std::cos(phi1);
        float z4 = r_out * std::sin(phi1);
        
        // Colors (Hot orange to red gradient)
        // Inner: Hot/White
        float r_in_c = 1.0f, g_in_c = 0.9f, b_in_c = 0.5f;
        // Outer: Red/Dark
        float r_out_c = 0.8f, g_out_c = 0.2f, b_out_c = 0.0f;
        
        float a = 0.6f;
        
        vertices.push_back(Vertex(x1, 0, z1, r_in_c, g_in_c, b_in_c, a));
        vertices.push_back(Vertex(x2, 0, z2, r_in_c, g_in_c, b_in_c, a));
        vertices.push_back(Vertex(x3, 0, z3, r_out_c, g_out_c, b_out_c, a));
        
        vertices.push_back(Vertex(x1, 0, z1, r_in_c, g_in_c, b_in_c, a));
        vertices.push_back(Vertex(x3, 0, z3, r_out_c, g_out_c, b_out_c, a));
        vertices.push_back(Vertex(x4, 0, z4, r_out_c, g_out_c, b_out_c, a));
    }
    return vertices;
}

std::vector<Vertex> generate_starfield(int num_stars, float radius) {
    std::vector<Vertex> vertices;
    
    for (int i = 0; i < num_stars; ++i) {
        // Random point on sphere
        float theta = std::acos(1.0f - 2.0f * float(rand()) / RAND_MAX);
        float phi = 2.0f * 3.14159f * float(rand()) / RAND_MAX;
        
        float x = radius * std::sin(theta) * std::cos(phi);
        float y = radius * std::cos(theta);
        float z = radius * std::sin(theta) * std::sin(phi);
        
        // Random brightness
        float brightness = 0.5f + 0.5f * float(rand()) / RAND_MAX;
        
        vertices.push_back(Vertex(x, y, z, brightness, brightness, brightness, 1.0f));
    }
    return vertices;
}

std::vector<Vertex> generate_geodesic_line(const std::vector<float>& points_x,
                                          const std::vector<float>& points_y,
                                          const std::vector<float>& points_z,
                                          const std::vector<float>& colors_r,
                                          const std::vector<float>& colors_g,
                                          const std::vector<float>& colors_b) {
    std::vector<Vertex> vertices;
    
    size_t n = points_x.size();
    vertices.reserve(n);
    
    for (size_t i = 0; i < n; ++i) {
        float r = colors_r.empty() ? 1.0f : colors_r[i];
        float g = colors_g.empty() ? 1.0f : colors_g[i];
        float b = colors_b.empty() ? 1.0f : colors_b[i];
        
        vertices.push_back(Vertex(points_x[i], points_y[i], points_z[i], r, g, b, 1.0f));
    }
    
    return vertices;
}

} // namespace Render
