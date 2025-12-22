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
"uniform float u_scale;\n"
"uniform vec4 u_color_tint;\n" // Added color tint
"\n"
"out vec4 v_color;\n"
"\n"
"void main() {\n"
"    gl_Position = u_mvp * vec4(a_position * u_scale, 1.0);\n"
"    v_color = a_color * u_color_tint;\n" // Apply tint
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

// --- BLOOM SHADERS ---

// 1. Bloom Extraction (Threshold)
const char* bloom_shader_src = 
"#version 300 es\n"
"precision mediump float;\n"
"out vec4 fragColor;\n"
"in vec2 v_texCoord;\n"
"uniform sampler2D u_image;\n"
"\n"
"void main() {\n"
"    vec4 color = texture(u_image, v_texCoord);\n"
"    // Calculate brightness (simple gray scale)\n"
"    float brightness = dot(color.rgb, vec3(0.2126, 0.7152, 0.0722));\n"
"    if(brightness > Physics::BLOOM_THRESHOLD) {\n"
"        fragColor = vec4(color.rgb, 1.0);\n"
"    } else {\n"
"        fragColor = vec4(0.0, 0.0, 0.0, 1.0);\n"
"    }\n"
"}\n";

// 2. Gaussian Blur (9-tap)
const char* blur_shader_src = 
"#version 300 es\n"
"precision mediump float;\n"
"out vec4 fragColor;\n"
"in vec2 v_texCoord;\n"
"uniform sampler2D u_image;\n"
"uniform bool u_horizontal;\n"
"uniform float u_weight[5];\n"
"\n"
"void main() {\n"
"    vec2 tex_offset = 1.0 / vec2(textureSize(u_image, 0));\n" // Gets size of single texel
"    vec3 result = texture(u_image, v_texCoord).rgb * u_weight[0];\n"
"    if(u_horizontal) {\n"
"        for(int i = 1; i < 5; ++i) {\n"
"            result += texture(u_image, v_texCoord + vec2(tex_offset.x * float(i), 0.0)).rgb * u_weight[i];\n"
"            result += texture(u_image, v_texCoord - vec2(tex_offset.x * float(i), 0.0)).rgb * u_weight[i];\n"
"        }\n"
"    } else {\n"
"        for(int i = 1; i < 5; ++i) {\n"
"            result += texture(u_image, v_texCoord + vec2(0.0, tex_offset.y * float(i))).rgb * u_weight[i];\n"
"            result += texture(u_image, v_texCoord - vec2(0.0, tex_offset.y * float(i))).rgb * u_weight[i];\n"
"        }\n"
"    }\n"
"    fragColor = vec4(result, 1.0);\n"
"}\n";

// 3. Composite (Scene + Bloom)
const char* composite_shader_src = 
"#version 300 es\n"
"precision mediump float;\n"
"out vec4 fragColor;\n"
"in vec2 v_texCoord;\n"
"uniform sampler2D u_scene;\n"
"uniform sampler2D u_bloom;\n"
"uniform float u_exposure;\n"
"\n"
"void main() {\n"
"    const float gamma = 2.2;\n"
"    vec3 hdrColor = texture(u_scene, v_texCoord).rgb;\n"
"    vec3 bloomColor = texture(u_bloom, v_texCoord).rgb;\n"
"    \n"
"    // Additive blending\n"
"    hdrColor += bloomColor;\n"
"    \n"
"    // Tone mapping (exposure)\n"
"    vec3 result = vec3(1.0) - exp(-hdrColor * u_exposure);\n"
"    \n"
"    // Gamma correction\n"
"    result = pow(result, vec3(1.0 / gamma));\n"
"    \n"
"    fragColor = vec4(result, 1.0);\n"
"}\n";

// Pass-through Vertex Shader for Quads
const char* quad_vertex_shader_src = 
"#version 300 es\n"
"layout (location = 0) in vec2 aPos;\n"
"layout (location = 1) in vec2 aTexCoords;\n"
"out vec2 v_texCoord;\n"
"\n"
"void main() {\n"
"    v_texCoord = aTexCoords;\n"
"    gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);\n"
"}\n";

Renderer::Renderer()
    : shader_program_(0), u_mvp_(-1), u_scale_(-1), vao_(0), 
      vbo_static_(0), 
      vbo_sphere_mesh_(0), count_sphere_mesh_(0),
      vbo_captured_(0), vbo_escaped_(0), vbo_other_(0),
      offset_accretion_disk_(0), count_accretion_disk_(0),
      offset_starfield_(0), count_starfield_(0),
      show_horizon_(true), show_photon_sphere_(true),
      show_accretion_disk_(true),
      show_starfield_(true),
      show_einstein_ring_(true),
      color_mode_(ColorMode::BY_TERMINATION),
      fbo_main_(0), tex_main_(0), rbo_depth_(0),
      fbo_bright_(0), tex_bright_(0),
      fbo_blur_{0, 0}, tex_blur_{0, 0},
      shader_bloom_(0), shader_blur_(0), shader_composite_(0),
      vao_quad_(0), vbo_quad_(0),
      enable_bloom_(false) {}

Renderer::~Renderer() {
    if (vbo_static_) glDeleteBuffers(1, &vbo_static_);
    if (vbo_sphere_mesh_) glDeleteBuffers(1, &vbo_sphere_mesh_);
    if (vbo_captured_) glDeleteBuffers(1, &vbo_captured_);
    if (vbo_escaped_) glDeleteBuffers(1, &vbo_escaped_);
    if (vbo_other_) glDeleteBuffers(1, &vbo_other_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (shader_program_) glDeleteProgram(shader_program_);
    
    if (enable_bloom_) {
        if (fbo_main_) glDeleteFramebuffers(1, &fbo_main_);
        if (tex_main_) glDeleteTextures(1, &tex_main_);
        if (rbo_depth_) glDeleteRenderbuffers(1, &rbo_depth_);
        if (fbo_bright_) glDeleteFramebuffers(1, &fbo_bright_);
        if (tex_bright_) glDeleteTextures(1, &tex_bright_);
        if (fbo_blur_[0]) glDeleteFramebuffers(2, fbo_blur_);
        if (tex_blur_[0]) glDeleteTextures(2, tex_blur_);
        if (shader_bloom_) glDeleteProgram(shader_bloom_);
        if (shader_blur_) glDeleteProgram(shader_blur_);
        if (shader_composite_) glDeleteProgram(shader_composite_);
        if (vao_quad_) glDeleteVertexArrays(1, &vao_quad_);
        if (vbo_quad_) glDeleteBuffers(1, &vbo_quad_);
    }
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
    u_scale_ = glGetUniformLocation(shader_program_, "u_scale");
    glUniform1f(u_scale_, 1.0f); // Default scale
    
    GLint u_color_tint = glGetUniformLocation(shader_program_, "u_color_tint");
    glUniform4f(u_color_tint, 1.0f, 1.0f, 1.0f, 1.0f); // Default tint (no change)
    
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
        // Use color based on current mode
        if (color_mode_ == ColorMode::BY_TERMINATION) {
            get_termination_color(geo, r, g, b);
        } else if (color_mode_ == ColorMode::BY_ERROR) {
            get_error_color(geo.max_H_error, r, g, b);
        } else if (color_mode_ == ColorMode::DOPPLER) {
            get_doppler_color(geo, r, g, b);
        } else if (color_mode_ == ColorMode::LENSING) {
            get_lensing_color(geo, r, g, b);
        } else { // SOLID
            r = 1.0f; g = 1.0f; b = 1.0f;
        }
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
    if (enable_bloom_) {
        if (tex_main_ == 0) {
            init_bloom_resources(width, height);
        } else {
            resize_bloom_resources(width, height);
        }
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_main_);
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    glViewport(0, 0, width, height);
    
    // Clear whatever is bound
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Setup for Scene Rendering
    glUseProgram(shader_program_);
    
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
    // u_scale_ is member
    GLint u_color_tint = glGetUniformLocation(shader_program_, "u_color_tint");
    
    // Reset tint
    glUniform4f(u_color_tint, 1.0f, 1.0f, 1.0f, 1.0f);
    
    // --- Draw Static Geometry ---
    // 1. Starfield (u_scale = 1.0)
    if (show_starfield_) {
        glUniform1f(u_scale_, 1.0f);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        
        glDrawArrays(GL_POINTS, offset_starfield_, count_starfield_);
    }
    
    // 2. Accretion Disk (blended)
    if (show_accretion_disk_) {
        glUniform1f(u_scale_, 1.0f);
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
            glUniform1f(u_scale_, Physics::R_SCHWARZSCHILD);
            // Horizon should be black. Mesh is white (1,1,1,1).
            // Tint: 0,0,0,1
            glUniform4f(u_color_tint, 0.0f, 0.0f, 0.0f, 1.0f);
            glDrawArrays(GL_LINES, 0, count_sphere_mesh_);
        }
        
        // Photon Sphere (Scale = R_PHOTON_SPHERE = 3.0)
        if (show_photon_sphere_) {
            glUniform1f(u_scale_, Physics::R_PHOTON_SPHERE);
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
    
    // Draw Dynamic Geometry
    glUniform1f(u_scale_, 1.0f); // Reset scale
    
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
        
        // Draw Interactive Ray (Task 30)
        if (!interactive_vertices_.empty()) {
            if (vbo_interactive_ == 0) glGenBuffers(1, &vbo_interactive_);
            glBindBuffer(GL_ARRAY_BUFFER, vbo_interactive_);
            glBufferData(GL_ARRAY_BUFFER, interactive_vertices_.size() * sizeof(Vertex), interactive_vertices_.data(), GL_DYNAMIC_DRAW);
            
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
            glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
            
            // Draw thick line
            // glLineWidth(3.0); // Not supported in all browsers but good to try
            glDrawArrays(GL_LINE_STRIP, 0, interactive_vertices_.size());
            // glLineWidth(1.0);
        }
        
        // Draw Einstein Ring Overlay
        if (show_einstein_ring_) {
            draw_einstein_ring(view_matrix, proj_matrix);
        }

        // Restore standard blending
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    
    // 2. Bloom Post-processing
    if (enable_bloom_) {
        render_bloom_pass(width, height);
    }
}

void Renderer::set_interactive_ray(const Numerics::Geodesic& ray) {
    interactive_vertices_.clear();
    if (ray.points.size() < 2) return;
    
    std::vector<float> points_x, points_y, points_z;
    points_x.reserve(ray.points.size());
    points_y.reserve(ray.points.size());
    points_z.reserve(ray.points.size());
    
    for (const auto& p : ray.points) {
        float x, y, z;
        to_cartesian(p.x[1], p.x[2], p.x[3], x, y, z);
        points_x.push_back(x);
        points_y.push_back(y);
        points_z.push_back(z);
    }
    
    // Bright Green for interactive ray
    std::vector<float> r(ray.points.size(), 0.2f);
    std::vector<float> g(ray.points.size(), 1.0f);
    std::vector<float> b(ray.points.size(), 0.2f);
    
    interactive_vertices_ = generate_geodesic_line(points_x, points_y, points_z, r, g, b);
}

void Renderer::draw_einstein_ring(const float view[16], const float proj[16]) {
    // Critical radius b = sqrt(27) approx 5.196
    const float radius = 5.196f;
    
    // Use the Sphere Mesh but scale it
    glUniform1f(u_scale_, radius); 
    
    // Set Tint to Gold to ensure visibility over white starfield/disk
    // We need to look up the uniform again as it's not a member, or we could add it.
    // Optimization: Store u_color_tint_ as member. But looking up is fine for now.
    GLint u_color_tint = glGetUniformLocation(shader_program_, "u_color_tint");
    glUniform4f(u_color_tint, 1.0f, 0.8f, 0.0f, 1.0f); // Opaque Gold
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_sphere_mesh_);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
    
    // Draw wireframe
    glDrawArrays(GL_LINES, 0, count_sphere_mesh_); 
    
    // Reset Tint
    glUniform4f(u_color_tint, 1.0f, 1.0f, 1.0f, 1.0f);
}

void Renderer::get_doppler_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const {
    if (geo.points.empty()) { r=1; g=1; b=1; return; }
    
    // Rigorous Relativistic Doppler
    // v_r (ray) < 0 => Photon arriving (v > 0)
    
    const auto& p = geo.points[0];
    double radius = p.x[1];
    double one_minus_2M = 1.0 - 2.0 * Physics::M / radius;
    double p_r_cov = p.p[1];
    double E = geo.E_initial;
    
    if (std::abs(E) < 1e-9) { r=1; g=1; b=1; return; }

    // Calculate radial velocity v = dr/dt = (p^r / p^t)
    double p_r_con = one_minus_2M * p_r_cov;
    double p_t_con = E / one_minus_2M;
    double v_r_ray = p_r_con / p_t_con;

    // Photon velocity relative to observer: v_photon = -v_ray (Ray back-traced)
    double v_photon = -v_r_ray;
    
    // Clamp
    if (v_photon > 0.999) v_photon = 0.999;
    if (v_photon < -0.999) v_photon = -0.999;
    
    // Doppler Factor D = sqrt((1-v)/(1+v))
    double D = std::sqrt((1.0 - v_photon) / (1.0 + v_photon));
    
    // Map D to Color
    // D < 1 (Redshift/Climbing out) -> Red
    // D > 1 (Blueshift/Falling in) -> Blue
    
    // Base Luminance
    float intensity = 1.0f;
    
    if (D < 1.0) {
        // Redshift (D goes 1->0)
        float t = std::pow(float(D), 0.5f); // Curve it
        r = 1.0f * intensity;
        g = t * intensity;
        b = t * intensity; // Fade to red
    } else {
        // Blueshift (D goes 1->Inf)
        // Cap at D=3
        double D_clamped = D > 3.0 ? 3.0 : D;
        float t = float(D_clamped - 1.0) / 2.0f; // 0 to 1
        r = (1.0f - t) * intensity;
        g = (1.0f - t) * intensity;
        b = 1.0f * intensity; // Fade to blue
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
            r = Physics::COLOR_RED_R; g = Physics::COLOR_RED_G; b = Physics::COLOR_RED_B;
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

// Gravitational Lensing Grid
void Renderer::get_lensing_color(const Numerics::Geodesic& geo, float& r, float& g, float& b) const {
    if (geo.termination != Numerics::TerminationReason::ESCAPED || geo.points.empty()) {
        r = 0.0f; g = 0.0f; b = 0.0f;
        return;
    }
    
    const auto& last_p = geo.points.back();
    double theta = last_p.x[Physics::THETA];
    double phi = last_p.x[Physics::PHI];
    
    double density = 10.0; 
    double u = theta / 3.14159;
    double v = phi / (2.0 * 3.14159);
    
    int check = (int(u * density) + int(v * density * 2.0)) % 2;
    
    if (check == 0) {
        r = 1.0f; g = 1.0f; b = 1.0f;
    } else {
        r = 0.2f; g = 0.2f; b = 0.2f;
    }
}

// Bloom Helpers (Internal linkage)
GLuint compile_shader_helper(const char* source, GLenum type) {
    GLuint shader = glCreateShader(type);
    if (!shader) return 0;
    
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        std::cerr << "Shader compile error: " << log << std::endl;
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

GLuint link_program_helper(GLuint vs, GLuint fs) {
    if (!vs || !fs) return 0;
    
    GLuint prog = glCreateProgram();
    if (!prog) return 0;
    
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint success;
    glGetProgramiv(prog, GL_LINK_STATUS, &success);
    if (!success) {
        char log[512];
        glGetProgramInfoLog(prog, 512, nullptr, log);
        std::cerr << "Program link error: " << log << std::endl;
        glDeleteProgram(prog);
        return 0;
    }
    return prog;
}

bool Renderer::init_bloom_resources(int width, int height) {
    GLuint vs_quad = compile_shader_helper(quad_vertex_shader_src, GL_VERTEX_SHADER);
    if (!vs_quad) return false;
    
    GLuint fs_bloom = compile_shader_helper(bloom_shader_src, GL_FRAGMENT_SHADER);
    if (!fs_bloom) {
        glDeleteShader(vs_quad);
        return false;
    }
    
    shader_bloom_ = link_program_helper(vs_quad, fs_bloom);
    
    GLuint fs_blur = compile_shader_helper(blur_shader_src, GL_FRAGMENT_SHADER);
    if (!fs_blur) {
        glDeleteShader(vs_quad);
        glDeleteShader(fs_bloom);
        if (shader_bloom_) glDeleteProgram(shader_bloom_);
        return false;
    }
    
    shader_blur_ = link_program_helper(vs_quad, fs_blur);
    
    GLuint fs_comp = compile_shader_helper(composite_shader_src, GL_FRAGMENT_SHADER);
    if (!fs_comp) {
        glDeleteShader(vs_quad);
        glDeleteShader(fs_bloom);
        glDeleteShader(fs_blur);
        if (shader_bloom_) glDeleteProgram(shader_bloom_);
        if (shader_blur_) glDeleteProgram(shader_blur_);
        return false;
    }
    
    shader_composite_ = link_program_helper(vs_quad, fs_comp);
    
    if (!shader_bloom_ || !shader_blur_ || !shader_composite_) {
        glDeleteShader(vs_quad);
        glDeleteShader(fs_bloom);
        glDeleteShader(fs_blur);
        glDeleteShader(fs_comp);
        if (shader_bloom_) glDeleteProgram(shader_bloom_);
        if (shader_blur_) glDeleteProgram(shader_blur_);
        if (shader_composite_) glDeleteProgram(shader_composite_);
        return false;
    }
    
    glDeleteShader(vs_quad);
    glDeleteShader(fs_bloom);
    glDeleteShader(fs_blur);
    glDeleteShader(fs_comp);
    
    glGenFramebuffers(1, &fbo_main_);
    glGenTextures(1, &tex_main_);
    glGenRenderbuffers(1, &rbo_depth_);
    
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_main_);
    glBindTexture(GL_TEXTURE_2D, tex_main_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_main_, 0);
    
    glBindRenderbuffer(GL_RENDERBUFFER, rbo_depth_);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo_depth_);
    
    glGenFramebuffers(2, fbo_blur_);
    glGenTextures(2, tex_blur_);
    for (int i = 0; i < 2; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur_[i]);
        glBindTexture(GL_TEXTURE_2D, tex_blur_[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_blur_[i], 0);
    }
    
    glGenVertexArrays(1, &vao_quad_);
    glGenBuffers(1, &vbo_quad_);
    glBindVertexArray(vao_quad_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_quad_);
    float quadVertices[] = {
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f,  0.0f, 1.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return true;
}

void Renderer::resize_bloom_resources(int width, int height) {
    if (!enable_bloom_ || tex_main_ == 0) return;
    
    glBindTexture(GL_TEXTURE_2D, tex_main_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    
    glBindRenderbuffer(GL_RENDERBUFFER, rbo_depth_);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    
    for (int i = 0; i < 2; i++) {
        glBindTexture(GL_TEXTURE_2D, tex_blur_[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
    }
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Renderer::render_quad() {
    glBindVertexArray(vao_quad_);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
}

void Renderer::render_bloom_pass(int width, int height) {
    if (!enable_bloom_) return;
    
    bool horizontal = true;
    int amount = Physics::BLOOM_BLUR_PASSES;
    
    glUseProgram(shader_blur_);
    float weights[5] = {0.227027f, 0.1945946f, 0.1216216f, 0.054054f, 0.016216f};
    glUniform1fv(glGetUniformLocation(shader_blur_, "u_weight"), 5, weights);
    
    for (int i = 0; i < amount; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur_[horizontal]); 
        glUniform1i(glGetUniformLocation(shader_blur_, "u_horizontal"), horizontal);
        
        if (i == 0) {
           glUseProgram(shader_bloom_);
           glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur_[0]);
           glClear(GL_COLOR_BUFFER_BIT);
           glActiveTexture(GL_TEXTURE0);
           glBindTexture(GL_TEXTURE_2D, tex_main_); 
           render_quad();
           
           glUseProgram(shader_blur_); 
           continue; 
        }
        
        glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur_[horizontal]); 
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex_blur_[!horizontal]); 
        render_quad();
        horizontal = !horizontal;
    }
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glUseProgram(shader_composite_);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex_main_);
    glUniform1i(glGetUniformLocation(shader_composite_, "u_scene"), 0);
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, tex_blur_[!horizontal]);
    glUniform1i(glGetUniformLocation(shader_composite_, "u_bloom"), 1);
    glUniform1f(glGetUniformLocation(shader_composite_, "u_exposure"), 1.0f);
    render_quad();
}

} // namespace Render
