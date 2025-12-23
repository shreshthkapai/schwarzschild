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
    DOPPLER,         // Task 27: Color by Doppler shift
    LENSING          // Task 28: Gravitational Lensing Grid
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
    
    // Task 30: Set interactive ray
    void set_interactive_ray(const Numerics::Geodesic& ray);
    
    // Update geodesics from serialized buffer (for worker integration)
    void update_geodesics_from_buffer(const std::vector<float>& data);
    
    // Clear dynamic geodesics
    void clear_geodesics();
    
    // Configuration
    void set_show_horizon(bool show) { show_horizon_ = show; }
    void set_show_photon_sphere(bool show) { show_photon_sphere_ = show; }
    void set_show_accretion_disk(bool show) { show_accretion_disk_ = show; }
    void set_show_starfield(bool show) { show_starfield_ = show; }
    void set_show_einstein_ring(bool show) { show_einstein_ring_ = show; } // Task 29
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
    GLint u_scale_; // Promoted to member
    
    // Vertex buffers
    GLuint vao_;
    GLuint vbo_static_;     // Merged static geometry (Stars, Disk)
    
    // Task 15: Shared sphere mesh for instancing
    GLuint vbo_sphere_mesh_;
    GLint count_sphere_mesh_;

    // Task 14: Separate VBOs for different termination types
    GLuint vbo_captured_;
    GLuint vbo_escaped_;
    GLuint vbo_other_;
    
    // Task 30: Interactive Ray VBO
    GLuint vbo_interactive_;
    std::vector<Vertex> interactive_vertices_;
    
    // Geometry offsets and counts in vbo_static_
    // Horizon and Photon Sphere are now drawn using vbo_sphere_mesh_
    GLint offset_accretion_disk_, count_accretion_disk_;
    GLint offset_starfield_, count_starfield_;
    
    // Geometry data (temporary storage during init)
    std::vector<Vertex> static_vertices_;
    
    // Task 14: Separate vertex arrays
    std::vector<Vertex> captured_vertices_;
    std::vector<Vertex> escaped_vertices_;
    std::vector<Vertex> other_vertices_;
    
    // Settings
    bool show_horizon_;
    bool show_photon_sphere_;
    bool show_accretion_disk_;
    bool show_starfield_;
    bool show_einstein_ring_; // Task 29
    ColorMode color_mode_;
    
    // Helper functions
    bool compile_shaders();
    void build_static_geometry(); // Merged generation
    void update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics);
    
    // Drawing helpers
    void draw_layout(GLenum mode, GLint first, GLint count);
    void draw_einstein_ring(const float view[16], const float proj[16]); // Task 29
    
    // Convert Schwarzschild (r, θ, φ) to Cartesian (x, y, z)
    void to_cartesian(double r, double theta, double phi, float& x, float& y, float& z) const;
    
    // Color mapping functions
    void get_error_color(double error, float& r, float& g, float& b) const;
    void get_termination_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
    // Task 27: Color by Doppler shift
    void get_doppler_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;
    
    // Task 28: Gravitational Lensing
    void get_lensing_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const;

    // Task 26: Bloom Post-processing Resources
    GLuint fbo_main_, tex_main_, rbo_depth_;
    GLuint fbo_bright_, tex_bright_;
    GLuint fbo_blur_[2], tex_blur_[2];
    
    GLuint shader_bloom_;     // Extract bright color
    GLuint shader_blur_;      // Gaussian blur
    GLuint shader_composite_; // Combine scene + bloom
    
    GLuint vao_quad_, vbo_quad_; // Full screen quad
    
    bool enable_bloom_; 
    
    // Bloom Helpers
    bool init_bloom_resources(int width, int height);
    void resize_bloom_resources(int width, int height);
    void render_quad();
    void render_bloom_pass(int width, int height);
    
public:
    void set_enable_bloom(bool enable) { enable_bloom_ = enable; }
    bool get_enable_bloom() const { return enable_bloom_; }
};

} // namespace Render

#endif // RENDERER_H