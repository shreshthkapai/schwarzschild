#ifndef RENDERER_H
#define RENDERER_H

#include "render/geometry.h"
#include "numerics/integrator.h"
#include "rays/ray_bundle.h"
#include <GLES3/gl3.h>
#include <vector>

namespace Render {

// Color mapping modes
enum class ColorMode {
    SOLID,           // Single color
    BY_ERROR,        // Color by Hamiltonian error
    BY_TERMINATION   // Color by capture/escape
};

class Renderer {
public:
    Renderer();
    ~Renderer();
    
    // Initialize GL resources
    bool initialize();
    
    // Render scene
    void render(int width, int height, const float view_matrix[16], const float proj_matrix[16]);
    
    // Set geodesics to render
    void set_geodesics(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Configuration
    void set_show_horizon(bool show) { show_horizon_ = show; }
    void set_show_photon_sphere(bool show) { show_photon_sphere_ = show; }
    void set_show_accretion_disk(bool show) { show_accretion_disk_ = show; }
    void set_show_starfield(bool show) { show_starfield_ = show; }
    void set_color_mode(ColorMode mode) { color_mode_ = mode; }
    
    // Getters for controls
    bool get_show_horizon() const { return show_horizon_; }
    bool get_show_photon_sphere() const { return show_photon_sphere_; }
    bool get_show_accretion_disk() const { return show_accretion_disk_; }
    bool get_show_starfield() const { return show_starfield_; }
    ColorMode get_color_mode() const { return color_mode_; }
    
private:
    // Shader programs
    GLuint shader_program_;
    GLint u_mvp_;
    
    // Vertex buffers
    GLuint vao_;
    GLuint vbo_horizon_;       // Solid sphere (triangles)
    GLuint vbo_photon_sphere_; // Wireframe sphere (lines)
    GLuint vbo_accretion_disk_;// Gradient disk (triangles)
    GLuint vbo_starfield_;     // Star points (points)
    GLuint vbo_geodesics_; // Single batched VBO
    
    // Geometry data
    std::vector<Vertex> horizon_vertices_;
    std::vector<Vertex> photon_sphere_vertices_;
    std::vector<Vertex> accretion_disk_vertices_;
    std::vector<Vertex> starfield_vertices_;
    std::vector<Vertex> geodesic_vertices_; // Flattened vectors
    
    // Settings
    bool show_horizon_;
    bool show_photon_sphere_;
    bool show_accretion_disk_;
    bool show_starfield_;
    ColorMode color_mode_;
    
    // Helper functions
    void compile_shaders();
    void update_horizon_geometry();
    void update_photon_sphere_geometry();
    void update_accretion_disk_geometry();
    void update_starfield_geometry();
    void update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Drawing helpers
    void draw_lines(const std::vector<Vertex>& vertices, GLuint vbo);
    void draw_triangles(const std::vector<Vertex>& vertices, GLuint vbo);
    void draw_points(const std::vector<Vertex>& vertices, GLuint vbo);
    
    // Convert Schwarzschild (r, θ, φ) to Cartesian (x, y, z)
    void to_cartesian(double r, double theta, double phi, float& x, float& y, float& z) const;
    
    // Color mapping functions
    void get_error_color(double error, float& r, float& g, float& b) const;
    void get_termination_color(Numerics::TerminationReason reason, float& r, float& g, float& b) const;
};

} // namespace Render

#endif // RENDERER_H