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
"layout(location = 0) in vec3 a_position;\n"
"layout(location = 1) in vec4 a_color;\n"
"\n"
"uniform mat4 u_mvp;\n"
"uniform float u_scale;\n" // Added scale uniform
"\n"
"out vec4 v_color;\n"
"\n"
"void main() {\n"
"    gl_Position = u_mvp * vec4(a_position * u_scale, 1.0);\n"
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
      vbo_static_(0), 
      vbo_sphere_mesh_(0), count_sphere_mesh_(0),
      vbo_captured_(0), vbo_escaped_(0), vbo_other_(0),
      offset_accretion_disk_(0), count_accretion_disk_(0),
      offset_starfield_(0), count_starfield_(0),
      show_horizon_(true), show_photon_sphere_(true),
      show_accretion_disk_(true), show_starfield_(true),
      color_mode_(ColorMode::BY_TERMINATION) {}

Renderer::~Renderer() {
    if (vbo_static_) glDeleteBuffers(1, &vbo_static_);
    if (vbo_sphere_mesh_) glDeleteBuffers(1, &vbo_sphere_mesh_);
    if (vbo_captured_) glDeleteBuffers(1, &vbo_captured_);
    if (vbo_escaped_) glDeleteBuffers(1, &vbo_escaped_);
    if (vbo_other_) glDeleteBuffers(1, &vbo_other_);
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
    GLint u_scale = glGetUniformLocation(shader_program_, "u_scale");
    glUniform1f(u_scale, 1.0f); // Default scale
    
    return true;
}

void Renderer::build_static_geometry() {
    // using namespace Geometry; // Geometry namespace does not exist
    
    static_vertices_.clear();
    
    // 1. Sphere Mesh (Unit Sphere) - Task 15
    // Shared mesh for Horizon and Photon Sphere
    std::vector<Vertex> sphere_verts = generate_unit_sphere_mesh(32);
    count_sphere_mesh_ = sphere_verts.size();
    
    // Upload Sphere VBO
    glGenBuffers(1, &vbo_sphere_mesh_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_mesh_);
    glBufferData(GL_ARRAY_BUFFER, sphere_verts.size() * sizeof(Vertex), sphere_verts.data(), GL_STATIC_DRAW);
    
    // 2. Accretion Disk (Flat Ring)
    // Note: Colors are hardcoded in generate_accretion_disk for gradient
    std::vector<Vertex> disk_verts = generate_accretion_disk(Physics::R_ISCO, 20.0f, 64, 0.05f);
    offset_accretion_disk_ = static_vertices_.size();
    count_accretion_disk_ = disk_verts.size();
    static_vertices_.insert(static_vertices_.end(), disk_verts.begin(), disk_verts.end());
    
    // 3. Starfield (Points at infinity)
    std::vector<Vertex> star_verts;
    // Random stars
    for(int i=0; i<1000; ++i) {
        float theta = (float(rand()) / RAND_MAX) * M_PI;
        float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
        float r = 500.0f;
        float x = r * sin(theta) * cos(phi);
        float y = r * cos(theta);
        float z = r * sin(theta) * sin(phi);
        star_verts.emplace_back(x, y, z, 1.0f, 1.0f, 1.0f, 1.0f);
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
    captured_vertices_.clear();
    escaped_vertices_.clear();
    other_vertices_.clear();
    
    for (const auto& geo : geodesics) {
        if (geo.points.size() < 2) continue;
        
        float r, g, b;
        get_termination_color(geo, r, g, b);
        float color_a = 0.8f; 
        
        // Select target vector based on termination
        std::vector<Vertex>* target_verts = &other_vertices_;
        if (geo.termination == Numerics::TerminationReason::HORIZON_CROSSED) {
            target_verts = &captured_vertices_;
        } else if (geo.termination == Numerics::TerminationReason::ESCAPED) {
            target_verts = &escaped_vertices_;
        }

        // Determine stride
        const size_t MAX_VERTICES = 50000;
        int stride = 1;
        size_t current_count = captured_vertices_.size() + escaped_vertices_.size() + other_vertices_.size();
        if (current_count > MAX_VERTICES) {
            stride = 4;
        } else if (current_count > MAX_VERTICES * 0.7) {
            stride = 2;
        }

        // Generate line segments with stride
        for (size_t i = 0; i < geo.points.size() - 1; i += stride) {
            size_t next_idx = i + stride;
            if (next_idx >= geo.points.size()) next_idx = geo.points.size() - 1;
            
            if (i >= next_idx) break;
            
            const auto& p1 = geo.points[i];
            const auto& p2 = geo.points[next_idx];
            
            float x1, y1, z1;
            float x2, y2, z2;
            
            to_cartesian(p1.x[Physics::R], p1.x[Physics::THETA], p1.x[Physics::PHI], x1, y1, z1);
            to_cartesian(p2.x[Physics::R], p2.x[Physics::THETA], p2.x[Physics::PHI], x2, y2, z2);
            
            target_verts->emplace_back(x1, y1, z1, r, g, b, color_a);
            target_verts->emplace_back(x2, y2, z2, r, g, b, color_a);
            
            if (next_idx == geo.points.size() - 1) break;
        }
    }
    
    // Re-upload VBOs
    if (vbo_captured_ == 0) glGenBuffers(1, &vbo_captured_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_captured_);
    glBufferData(GL_ARRAY_BUFFER, captured_vertices_.size() * sizeof(Vertex), captured_vertices_.data(), GL_DYNAMIC_DRAW);
    
    if (vbo_escaped_ == 0) glGenBuffers(1, &vbo_escaped_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_escaped_);
    glBufferData(GL_ARRAY_BUFFER, escaped_vertices_.size() * sizeof(Vertex), escaped_vertices_.data(), GL_DYNAMIC_DRAW);
    
    if (vbo_other_ == 0) glGenBuffers(1, &vbo_other_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_other_);
    glBufferData(GL_ARRAY_BUFFER, other_vertices_.size() * sizeof(Vertex), other_vertices_.data(), GL_DYNAMIC_DRAW);
}

void Renderer::update_geodesics_from_buffer(const std::vector<float>& data) {
    if (data.empty()) return;
    
    // Task 17: Vertex budget
    const size_t MAX_VERTICES = 50000;
    
    // Clear previous
    captured_vertices_.clear();
    escaped_vertices_.clear();
    other_vertices_.clear();
    
    // Parse buffer
    size_t idx = 0;
    while (idx < data.size()) {
        // Header: [ray_id, termination, num_points]
        if (idx + 3 > data.size()) break;
        
        idx++; // ray_id
        float termination = data[idx++];
        int num_points = (int)data[idx++];
        
        if (num_points < 2) {
            idx += num_points * 3;
            continue;
        }
        
        // Select target vector
        int term_code = (int)termination;
        std::vector<Vertex>* target_verts = &other_vertices_;
        if (term_code == 0) target_verts = &captured_vertices_; // Horizon
        else if (term_code == 1) target_verts = &escaped_vertices_; // Escaped
        
        // Determine stride
        int stride = 1;
        size_t current_count = captured_vertices_.size() + escaped_vertices_.size() + other_vertices_.size();
        if (current_count > MAX_VERTICES) {
            stride = 4;
        } else if (current_count > MAX_VERTICES * 0.7) {
            stride = 2;
        }
        
        // Color
        float r, g, b;
        if (term_code == 0) { // Horizon
            r = 0.0f; g = 0.0f; b = 0.0f;
        } else if (term_code == 1) { // Escape
            r = 0.2f; g = 0.4f; b = 0.8f;
        } else {
            r = 1.0f; g = 0.0f; b = 0.0f;
        }
        float color_a = 0.8f;
        
        // Read first point
        float x1 = data[idx++];
        float y1 = data[idx++];
        float z1 = data[idx++];
        
        for (int i = 1; i < num_points; ++i) {
            float x2 = data[idx++];
            float y2 = data[idx++];
            float z2 = data[idx++];
            
            if (i % stride == 0 || i == num_points - 1) {
                target_verts->emplace_back(x1, y1, z1, r, g, b, color_a);
                target_verts->emplace_back(x2, y2, z2, r, g, b, color_a);
                
                x1 = x2;
                y1 = y2;
                z1 = z2;
            }
        }
    }
    
    // Re-upload VBOs
    if (vbo_captured_ == 0) glGenBuffers(1, &vbo_captured_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_captured_);
    glBufferData(GL_ARRAY_BUFFER, captured_vertices_.size() * sizeof(Vertex), captured_vertices_.data(), GL_DYNAMIC_DRAW);
    
    if (vbo_escaped_ == 0) glGenBuffers(1, &vbo_escaped_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_escaped_);
    glBufferData(GL_ARRAY_BUFFER, escaped_vertices_.size() * sizeof(Vertex), escaped_vertices_.data(), GL_DYNAMIC_DRAW);
    
    if (vbo_other_ == 0) glGenBuffers(1, &vbo_other_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_other_);
    glBufferData(GL_ARRAY_BUFFER, other_vertices_.size() * sizeof(Vertex), other_vertices_.data(), GL_DYNAMIC_DRAW);
}

void Renderer::clear_geodesics() {
    captured_vertices_.clear();
    escaped_vertices_.clear();
    other_vertices_.clear();
    
    if (vbo_captured_) { glBindBuffer(GL_ARRAY_BUFFER, vbo_captured_); glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW); }
    if (vbo_escaped_) { glBindBuffer(GL_ARRAY_BUFFER, vbo_escaped_); glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW); }
    if (vbo_other_) { glBindBuffer(GL_ARRAY_BUFFER, vbo_other_); glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW); }
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
    
    // Get uniforms
    GLint u_scale = glGetUniformLocation(shader_program_, "u_scale");
    GLint u_color_tint = glGetUniformLocation(shader_program_, "u_color_tint");
    
    // Reset tint
    glUniform4f(u_color_tint, 1.0f, 1.0f, 1.0f, 1.0f);
    
    // --- Draw Static Geometry ---
    // 1. Starfield (u_scale = 1.0)
    if (show_starfield_) {
        glUniform1f(u_scale, 1.0f);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        
        glDrawArrays(GL_POINTS, offset_starfield_, count_starfield_);
    }
    
    // 2. Accretion Disk (blended)
    if (show_accretion_disk_) {
        glUniform1f(u_scale, 1.0f);
        if (!show_starfield_) glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
        
        glDepthMask(GL_FALSE);
        glDrawArrays(GL_TRIANGLES, offset_accretion_disk_, count_accretion_disk_);
        glDepthMask(GL_TRUE);
    }
    
    // 3. Spheres (Horizon and Photon Sphere) - Task 15 (Reusing Sphere Mesh)
    if (show_horizon_ || show_photon_sphere_) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_mesh_);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        
        // Horizon (Scale = R_SCHWARZSCHILD = 2.0)
        if (show_horizon_) {
            glUniform1f(u_scale, Physics::R_SCHWARZSCHILD);
            // Horizon should be black. Mesh is white (1,1,1,1).
            // Tint: 0,0,0,1
            glUniform4f(u_color_tint, 0.0f, 0.0f, 0.0f, 1.0f);
            glDrawArrays(GL_LINES, 0, count_sphere_mesh_);
        }
        
        // Photon Sphere (Scale = R_PHOTON_SPHERE = 3.0)
        if (show_photon_sphere_) {
            glUniform1f(u_scale, Physics::R_PHOTON_SPHERE);
            // Photon Sphere should be Yellow/Transparent.
            // Tint: 1,1,0,0.3
            glUniform4f(u_color_tint, 1.0f, 1.0f, 0.0f, 0.3f);
            
            glDepthMask(GL_FALSE); // Transparent
            glDrawArrays(GL_LINES, 0, count_sphere_mesh_);
            glDepthMask(GL_TRUE);
        }
        
        // Reset tint
        glUniform4f(u_color_tint, 1.0f, 1.0f, 1.0f, 1.0f);
    }
    
    // --- Draw Dynamic Geometry (Task 14: Separate VBOs) ---
    glUniform1f(u_scale, 1.0f); // Reset scale
    
    if (!captured_vertices_.empty() || !escaped_vertices_.empty() || !other_vertices_.empty()) {
        // Standard alpha blending for visibility
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Draw Escaped (Blue-ish)
        if (!escaped_vertices_.empty()) {
            glBindBuffer(GL_ARRAY_BUFFER, vbo_escaped_);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
             glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
             glDrawArrays(GL_LINES, 0, escaped_vertices_.size());
        }
        
        // Draw Captured (Black/Invisible?)
        if (!captured_vertices_.empty()) {
            glBindBuffer(GL_ARRAY_BUFFER, vbo_captured_);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
            glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
            glDrawArrays(GL_LINES, 0, captured_vertices_.size());
        }
        
        // Draw Other (Yellow)
        if (!other_vertices_.empty()) {
            glBindBuffer(GL_ARRAY_BUFFER, vbo_other_);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
            glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
            glDrawArrays(GL_LINES, 0, other_vertices_.size());
        }
        
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
            r = 1.0f; g = 0.4f; b = 0.1f; // Orange-red (visible, fell in)
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
                        r = 0.3f; g = 0.7f; b = 1.0f; // Bright Cyan
                    } else {
                         r = 0.2f; g = 0.5f; b = 0.9f; // Lighter Blue
                    }
                    
                    // Add a "Milky Way" band
                    if (std::abs(p_th) < 0.2f) {
                        r += 0.4f; g += 0.3f; b += 0.4f;
                    }
                    
                } else {
                    r = 0.3f; g = 0.6f; b = 1.0f; // Bright blue fallback
                }
            }
            break;
        case Numerics::TerminationReason::MAX_LAMBDA:
            r = 1.0f; g = 1.0f; b = 0.0f; // Yellow
            break;
    }
}

} // namespace Render
