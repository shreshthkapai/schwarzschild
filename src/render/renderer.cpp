#include "render/renderer.h"
#include "physics/constants.h"
#include <cmath>
#include <iostream>

namespace Render {

// Simple vertex shader
const char* vertex_shader_src = 
"#version 300 es\n"
"precision mediump float;\n"
"\n"
"in vec3 a_position;\n"
"in vec4 a_color;\n"
"\n"
"uniform mat4 u_mvp;\n"
"\n"
"out vec4 v_color;\n"
"\n"
"void main() {\n"
"    gl_Position = u_mvp * vec4(a_position, 1.0);\n"
"    v_color = a_color;\n"
"}\n";

// Simple fragment shader
const char* fragment_shader_src = 
"#version 300 es\n"
"precision mediump float;\n"
"\n"
"in vec4 v_color;\n"
"out vec4 fragColor;\n"
"\n"
"void main() {\n"
"    fragColor = v_color;\n"
"}\n";

Renderer::Renderer()
    : shader_program_(0), u_mvp_(-1), vao_(0), 
      vbo_static_(0), vbo_geodesics_(0),
      offset_horizon_(0), count_horizon_(0),
      offset_photon_sphere_(0), count_photon_sphere_(0),
      offset_accretion_disk_(0), count_accretion_disk_(0),
      offset_starfield_(0), count_starfield_(0),
      show_horizon_(true), show_photon_sphere_(true),
      show_accretion_disk_(true), show_starfield_(true),
      color_mode_(ColorMode::BY_TERMINATION) {}

Renderer::~Renderer() {
    if (vbo_static_) glDeleteBuffers(1, &vbo_static_);
    if (vbo_geodesics_) glDeleteBuffers(1, &vbo_geodesics_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (shader_program_) glDeleteProgram(shader_program_);
}

bool Renderer::initialize() {
    if (!compile_shaders()) return false;
    
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);
    
    build_static_geometry();
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(2.0f);
    
    return true;
}

bool Renderer::compile_shaders() {
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertex_shader_src, nullptr);
    glCompileShader(vs);
    
    GLint success;
    glGetShaderiv(vs, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(vs, 512, nullptr, log);
        std::cerr << "Vertex shader compilation failed:\n" << log << std::endl;
        return false;
    }
    
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragment_shader_src, nullptr);
    glCompileShader(fs);
    
    glGetShaderiv(fs, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(fs, 512, nullptr, log);
        std::cerr << "Fragment shader compilation failed:\n" << log << std::endl;
        return false;
    }
    
    shader_program_ = glCreateProgram();
    glAttachShader(shader_program_, vs);
    glAttachShader(shader_program_, fs);
    glLinkProgram(shader_program_);
    
    glGetProgramiv(shader_program_, GL_LINK_STATUS, &success);
    if (!success) {
        char log[512];
        glGetProgramInfoLog(shader_program_, 512, nullptr, log);
        std::cerr << "Shader program linking failed:\n" << log << std::endl;
        return false;
    }
    
    glUseProgram(shader_program_);
    u_mvp_ = glGetUniformLocation(shader_program_, "u_mvp");
    
    return true;
}

void Renderer::build_static_geometry() {
    using namespace Geometry;
    
    static_vertices_.clear();
    
    // 1. Horizon (Black Sphere)
    std::vector<Vertex> horizon_verts;
    generate_sphere(Physics::R_SCHWARZSCHILD, 32, 32, {0.0f, 0.0f, 0.0f, 1.0f}, horizon_verts);
    offset_horizon_ = 0;
    count_horizon_ = horizon_verts.size();
    static_vertices_.insert(static_vertices_.end(), horizon_verts.begin(), horizon_verts.end());
    
    // 2. Photon Sphere (Wireframe/Transparent)
    std::vector<Vertex> photon_verts;
    // Using points for now, or lines if I can generate lines
    // For simplicity, let's use a sphere and draw as points or lines
    generate_sphere(Physics::R_PHOTON_SPHERE, 32, 32, {1.0f, 1.0f, 0.0f, 0.3f}, photon_verts);
    offset_photon_sphere_ = static_vertices_.size();
    count_photon_sphere_ = photon_verts.size();
    static_vertices_.insert(static_vertices_.end(), photon_verts.begin(), photon_verts.end());
    
    // 3. Accretion Disk (Flat Ring)
    std::vector<Vertex> disk_verts;
    generate_disk(Physics::R_ISCO, 20.0, 64, {0.8f, 0.4f, 0.1f, 0.4f}, {0.4f, 0.1f, 0.0f, 0.0f}, disk_verts);
    offset_accretion_disk_ = static_vertices_.size();
    count_accretion_disk_ = disk_verts.size();
    static_vertices_.insert(static_vertices_.end(), disk_verts.begin(), disk_verts.end());
    
    // 4. Starfield (Points at infinity)
    std::vector<Vertex> star_verts;
    // Random stars
    for(int i=0; i<1000; ++i) {
        float theta = (float(rand()) / RAND_MAX) * M_PI;
        float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
        float r = 500.0f;
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        star_verts.push_back({{x,y,z}, {1.0f, 1.0f, 1.0f, 1.0f}});
    }
    offset_starfield_ = static_vertices_.size();
    count_starfield_ = star_verts.size();
    static_vertices_.insert(static_vertices_.end(), star_verts.begin(), star_verts.end());
    
    // Upload Static VBO
    glGenBuffers(1, &vbo_static_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
    glBufferData(GL_ARRAY_BUFFER, static_vertices_.size() * sizeof(Vertex), static_vertices_.data(), GL_STATIC_DRAW);
    
    // Clear temporary data
    static_vertices_.clear();
}

void Renderer::update_geodesics(const std::vector<Numerics::Geodesic>& geodesics) {
    update_geodesic_geometry(geodesics);
}

void Renderer::update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics) {
    geodesic_vertices_.clear();
    
    for (const auto& geo : geodesics) {
        if (geo.points.size() < 2) continue;
        
        float r, g, b;
        get_termination_color(geo, r, g, b);
        Vertex::Color color = {r, g, b, 0.8f}; // High alpha for visibility
        
        // Generate line segments
        for (size_t i = 0; i < geo.points.size() - 1; ++i) {
            const auto& p1 = geo.points[i];
            const auto& p2 = geo.points[i+1];
            
            float x1, y1, z1;
            float x2, y2, z2;
            
            to_cartesian(p1.x[Physics::R], p1.x[Physics::THETA], p1.x[Physics::PHI], x1, y1, z1);
            to_cartesian(p2.x[Physics::R], p2.x[Physics::THETA], p2.x[Physics::PHI], x2, y2, z2);
            
            geodesic_vertices_.push_back({{x1, y1, z1}, color});
            geodesic_vertices_.push_back({{x2, y2, z2}, color});
        }
    }
    
    if (vbo_geodesics_ == 0) {
        glGenBuffers(1, &vbo_geodesics_);
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
    // Use DYNAMIC_DRAW for frequent updates
    glBufferData(GL_ARRAY_BUFFER, geodesic_vertices_.size() * sizeof(Vertex), geodesic_vertices_.data(), GL_DYNAMIC_DRAW);
}

void Renderer::update_geodesics_from_buffer(const std::vector<float>& data) {
    if (data.empty()) return;
    
    // Parse buffer
    size_t idx = 0;
    while (idx < data.size()) {
        // Header: [ray_id, termination, num_points]
        if (idx + 3 > data.size()) break;
        
        // float ray_id = data[idx++]; // Unused for rendering
        idx++; 
        float termination = data[idx++];
        int num_points = (int)data[idx++];
        
        if (num_points < 2) {
            idx += num_points * 3;
            continue;
        }
        
        // Color
        float r, g, b;
        // Simple color based on termination
        // 0: Horizon, 1: Escape, 2: Other
        int term_code = (int)termination;
        if (term_code == 0) { // Horizon
            r = 0.0f; g = 0.0f; b = 0.0f;
        } else if (term_code == 1) { // Escape
            r = 0.2f; g = 0.4f; b = 0.8f; // Blue-ish
        } else {
            r = 1.0f; g = 0.0f; b = 0.0f; // Error/Max steps
        }
        Vertex::Color color = {r, g, b, 0.8f};
        
        // Points: [x, y, z] * num_points
        // We need to generate line segments (p1, p2), (p2, p3)...
        
        // Read first point
        float x1 = data[idx++];
        float y1 = data[idx++];
        float z1 = data[idx++];
        
        for (int i = 0; i < num_points - 1; ++i) {
            float x2 = data[idx++];
            float y2 = data[idx++];
            float z2 = data[idx++];
            
            geodesic_vertices_.push_back({{x1, y1, z1}, color});
            geodesic_vertices_.push_back({{x2, y2, z2}, color});
            
            x1 = x2;
            y1 = y2;
            z1 = z2;
        }
    }
    
    // Re-upload VBO
    if (vbo_geodesics_ == 0) {
        glGenBuffers(1, &vbo_geodesics_);
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
    glBufferData(GL_ARRAY_BUFFER, geodesic_vertices_.size() * sizeof(Vertex), geodesic_vertices_.data(), GL_DYNAMIC_DRAW);
}

void Renderer::clear_geodesics() {
    geodesic_vertices_.clear();
    if (vbo_geodesics_) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    }
}

void Renderer::render(int width, int height, const float view_matrix[16], const float proj_matrix[16]) {
    glViewport(0, 0, width, height);
    
    // Compute MVP = Proj * View
    // Note: This ignores Model matrix (Identity)
    float mvp[16];
    
    // Matrix multiplication: MVP = P * V
    // OpenGL is Column-Major
    for (int i = 0; i < 4; ++i) { // Row
        for (int j = 0; j < 4; ++j) { // Col
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += proj_matrix[k*4 + i] * view_matrix[j*4 + k]; 
            }
            mvp[i + 4*j] = sum;
        }
    }
    
    glUniformMatrix4fv(u_mvp_, 1, GL_FALSE, mvp);
    
    // --- Draw Static Geometry (Single VBO bind) ---
    glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
    
    // 1. Starfield
    if (show_starfield_) {
        glDrawArrays(GL_POINTS, offset_starfield_, count_starfield_);
    }
    
    // 2. Accretion Disk (blended)
    if (show_accretion_disk_) {
        glDepthMask(GL_FALSE);
        glDrawArrays(GL_TRIANGLES, offset_accretion_disk_, count_accretion_disk_);
        glDepthMask(GL_TRUE);
    }
    
    // 3. Horizon
    if (show_horizon_) {
        glDrawArrays(GL_TRIANGLES, offset_horizon_, count_horizon_);
    }
    
    // 4. Photon Sphere
    if (show_photon_sphere_) {
        // Use GL_POINTS or GL_LINES depending on generation
        // For sphere points:
        glDrawArrays(GL_POINTS, offset_photon_sphere_, count_photon_sphere_);
    }
    
    // --- Draw Dynamic Geometry ---
    if (!geodesic_vertices_.empty()) {
        // Additive blending for glow effect
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        glDrawArrays(GL_LINES, 0, geodesic_vertices_.size());
        
        // Restore standard blending
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
}

void Renderer::draw_layout(GLenum mode, GLint first, GLint count) {
    glDrawArrays(mode, first, count);
}

void Renderer::to_cartesian(double r, double theta, double phi, float& x, float& y, float& z) const {
    x = r * std::sin(theta) * std::cos(phi);
    y = r * std::cos(theta);
    z = r * std::sin(theta) * std::sin(phi);
}

void Renderer::get_error_color(double error, float& r, float& g, float& b) const {
    double log_error = std::log10(std::max(error, 1e-16));
    double t = (log_error + 16.0) / 10.0;
    t = std::max(0.0, std::min(1.0, t));
    r = t; g = 1.0 - 0.5*t; b = 0.0;
}

void Renderer::get_termination_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const {
    using namespace Numerics;
    switch(geo.termination) {
        case TerminationReason::HORIZON_CROSSED:
            r = 0.0f; g = 0.0f; b = 0.0f; // Black (fell in)
            break;
        case TerminationReason::ESCAPED:
            {
                // Procedural Sky / Einstein Ring effect
                if (!geo.points.empty()) {
                    const auto& last_p = geo.points.back().p;
                    float p_th = float(last_p[Physics::THETA]);
                    float p_ph = float(last_p[Physics::PHI]);
                    
                    // Checkerboard pattern
                    float u = p_th * 10.0f;
                    float v = p_ph * 10.0f;
                    bool check = (int(std::abs(u)) + int(std::abs(v))) % 2 == 0;
                    
                    if (check) {
                        r = 0.1f; g = 0.2f; b = 0.4f; // Dark Blue
                    } else {
                         r = 0.05f; g = 0.1f; b = 0.2f; // Darker
                    }
                    
                    // Add a "Milky Way" band
                    if (std::abs(p_th) < 0.2f) {
                        r += 0.4f; g += 0.3f; b += 0.4f;
                    }
                    
                } else {
                    r = 0.0f; g = 0.1f; b = 0.3f; // Background
                }
            }
            break;
        case TerminationReason::MAX_STEPS_EXCEEDED:
        default:
            r = 0.5f; g = 0.0f; b = 0.0f; // Red debug
            break;
    }
}

} // namespace Render
