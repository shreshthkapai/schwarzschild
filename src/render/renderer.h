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
    BY_TERMINATION,  // Color by capture/escape
    DOPPLER,         // Doppler shift
    LENSING          // Lensing grid
};

class Renderer {
public:
    Renderer();
    ~Renderer();
    
    // Resource initialization
    bool initialize();
    
    // Main render pass
    void render(int width, int height, const float view_matrix[16], const float proj_matrix[16]);
    
    // Update geodesics
    void update_geodesics(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Interactive ray
    void set_interactive_ray(const Numerics::Geodesic& ray);
    
    // Buffer update
    void update_geodesics_from_buffer(const std::vector<float>& data);
    
    // Clear dynamic data
    void clear_geodesics();
    
    // Configuration
    void set_show_horizon(bool show) { show_horizon_ = show; }
    void set_show_photon_sphere(bool show) { show_photon_sphere_ = show; }
    void set_show_accretion_disk(bool show) { show_accretion_disk_ = show; }
    void set_show_starfield(bool show) { show_starfield_ = show; }
    void set_show_einstein_ring(bool show) { show_einstein_ring_ = show; }
    void set_color_mode(ColorMode mode) { color_mode_ = mode; }
    
    // Getters for controls
    bool get_show_horizon() const { return show_horizon_; }
    bool get_show_photon_sphere() const { return show_photon_sphere_; }
    bool get_show_accretion_disk() const { return show_accretion_disk_; }
    bool get_show_starfield() const { return show_starfield_; }
    bool get_show_einstein_ring() const { return show_einstein_ring_; }
    ColorMode get_color_mode() const { return color_mode_; }
    
private:
    // Shader programs
    GLuint shader_program_;
    GLint u_mvp_;
    GLint u_scale_; 
    
    // Vertex buffers
    GLuint vao_;
    GLuint vbo_static_;     // Static vbo
    
    // Sphere instancing
    GLuint vbo_sphere_mesh_;
    GLint count_sphere_mesh_;

    // Termination vbos
    GLuint vbo_captured_;
    GLuint vbo_escaped_;
    GLuint vbo_other_;
    
    // Interactive vbo
    GLuint vbo_interactive_;
    std::vector<Vertex> interactive_vertices_;
    
    GLint offset_accretion_disk_, count_accretion_disk_;
    GLint offset_starfield_, count_starfield_;
    
    // Geometry data
    std::vector<Vertex> static_vertices_;
    
    // Termination geometry data
    std::vector<Vertex> captured_vertices_;
    std::vector<Vertex> escaped_vertices_;
    std::vector<Vertex> other_vertices_;
    
    // Settings
    bool show_horizon_;
    bool show_photon_sphere_;
    bool show_accretion_disk_;
    bool show_starfield_;
    bool show_einstein_ring_;
    ColorMode color_mode_;
    
    // Internal helpers
    bool compile_shaders();
    void build_static_geometry();
    void update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Drawing dispatch
    void draw_layout(GLenum mode, GLint first, GLint count);
    void draw_einstein_ring(const float view[16], const float proj[16]);
    
    // Convert Schwarzschild (r, θ, φ) to Cartesian (x, y, z)
    void to_cartesian(double r, double theta, double phi, float& x, float& y, float& z) const;
    
    // Color mapping functions
    void get_error_color(double error, float& r, float& g, float& b) const;
    void get_termination_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
    void get_doppler_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
    void get_lensing_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
};

} // namespace Render

#endif // RENDERER_H