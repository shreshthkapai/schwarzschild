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
    void update_geodesics(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Update geodesics from serialized buffer (for worker integration)
    void update_geodesics_from_buffer(const std::vector<float>& data);
    
    // Clear dynamic geodesics
    void clear_geodesics();
    
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
    GLuint vbo_static_;     // Merged static geometry
    GLuint vbo_geodesics_;  // Single batched VBO for dynamic rays
    
    // Geometry offsets and counts in vbo_static_
    GLint offset_horizon_, count_horizon_;
    GLint offset_photon_sphere_, count_photon_sphere_;
    GLint offset_accretion_disk_, count_accretion_disk_;
    GLint offset_starfield_, count_starfield_;
    
    // Geometry data (temporary storage during init)
    std::vector<Vertex> static_vertices_;
    std::vector<Vertex> geodesic_vertices_; 
    
    // Settings
    bool show_horizon_;
    bool show_photon_sphere_;
    bool show_accretion_disk_;
    bool show_starfield_;
    ColorMode color_mode_;
    
    // Helper functions
    bool compile_shaders();
    void build_static_geometry(); // Merged generation
    void update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Drawing helpers
    void draw_layout(GLenum mode, GLint first, GLint count);
    
    // Convert Schwarzschild (r, θ, φ) to Cartesian (x, y, z)
    void to_cartesian(double r, double theta, double phi, float& x, float& y, float& z) const;
    
    // Color mapping functions
    void get_error_color(double error, float& r, float& g, float& b) const;
    void get_termination_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
};

} // namespace Render

#endif // RENDERER_H