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
      vbo_horizon_(0), vbo_photon_sphere_(0),
      vbo_accretion_disk_(0), vbo_starfield_(0), vbo_geodesics_(0),
      show_horizon_(true), show_photon_sphere_(true),
      show_accretion_disk_(true), show_starfield_(true),
      color_mode_(ColorMode::BY_TERMINATION) {}

Renderer::~Renderer() {
    if (vbo_horizon_) glDeleteBuffers(1, &vbo_horizon_);
    if (vbo_photon_sphere_) glDeleteBuffers(1, &vbo_photon_sphere_);
    if (vbo_accretion_disk_) glDeleteBuffers(1, &vbo_accretion_disk_);
    if (vbo_starfield_) glDeleteBuffers(1, &vbo_starfield_);
    if (vbo_geodesics_) glDeleteBuffers(1, &vbo_geodesics_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (shader_program_) glDeleteProgram(shader_program_);
}

bool Renderer::initialize() {
    compile_shaders();
    
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);
    
    update_horizon_geometry();
    update_photon_sphere_geometry();
    update_accretion_disk_geometry();
    update_starfield_geometry();
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(2.0f);
    
    return true;
}

void Renderer::compile_shaders() {
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertex_shader_src, nullptr);
    glCompileShader(vs);
    
    GLint success;
    glGetShaderiv(vs, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(vs, 512, nullptr, log);
        std::cerr << "Vertex shader compilation failed:\n" << log << std::endl;
    }
    
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragment_shader_src, nullptr);
    glCompileShader(fs);
    
    glGetShaderiv(fs, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[512];
        glGetShaderInfoLog(fs, 512, nullptr, log);
        std::cerr << "Fragment shader compilation failed:\n" << log << std::endl;
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
    }
    
    u_mvp_ = glGetUniformLocation(shader_program_, "u_mvp");
    
    glDeleteShader(vs);
    glDeleteShader(fs);
}

void Renderer::update_horizon_geometry() {
    // Solid black sphere for event horizon
    horizon_vertices_ = generate_solid_sphere(Physics::R_SCHWARZSCHILD, 40, 0.0f, 0.0f, 0.0f, 1.0f);
    
    if (!vbo_horizon_) glGenBuffers(1, &vbo_horizon_);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_horizon_);
    glBufferData(GL_ARRAY_BUFFER, 
                 horizon_vertices_.size() * sizeof(Vertex),
                 horizon_vertices_.data(), 
                 GL_STATIC_DRAW);
}

void Renderer::update_photon_sphere_geometry() {
    // Wireframe yellow sphere for photon sphere (unstable orbit)
    photon_sphere_vertices_ = generate_sphere(Physics::R_PHOTON_SPHERE, 30, 1.0f, 0.8f, 0.0f, 0.3f);
    
    if (!vbo_photon_sphere_) glGenBuffers(1, &vbo_photon_sphere_);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_photon_sphere_);
    glBufferData(GL_ARRAY_BUFFER,
                 photon_sphere_vertices_.size() * sizeof(Vertex),
                 photon_sphere_vertices_.data(),
                 GL_STATIC_DRAW);
}

void Renderer::update_accretion_disk_geometry() {
    // Glowing accretion disk: inner=4M, outer=12M
    accretion_disk_vertices_ = generate_accretion_disk(4.0f, 12.0f, 60, 0.0f);
    
    if (!vbo_accretion_disk_) glGenBuffers(1, &vbo_accretion_disk_);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_accretion_disk_);
    glBufferData(GL_ARRAY_BUFFER,
                 accretion_disk_vertices_.size() * sizeof(Vertex),
                 accretion_disk_vertices_.data(),
                 GL_STATIC_DRAW);
}

void Renderer::update_starfield_geometry() {
    // Background stars at distance 800M
    // Increased count to 2000 to drastically improve "silhouette" effect of the black hole
    starfield_vertices_ = generate_starfield(500, 800.0f);
    
    if (!vbo_starfield_) glGenBuffers(1, &vbo_starfield_);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_starfield_);
    glBufferData(GL_ARRAY_BUFFER,
                 starfield_vertices_.size() * sizeof(Vertex),
                 starfield_vertices_.data(),
                 GL_STATIC_DRAW);
}

void Renderer::set_geodesics(const std::vector<Numerics::Geodesic>& geodesics) {
    update_geodesic_geometry(geodesics);
}

void Renderer::update_geodesic_geometry(const std::vector<Numerics::Geodesic>& geodesics) {
    geodesic_vertices_.clear();
    
    // Reserve estimated memory to prevent reallocations
    size_t est_vertices = geodesics.size() * 100; // Conservative estimate
    geodesic_vertices_.reserve(est_vertices);

    for (const auto& geo : geodesics) {
        if (geo.points.size() < 2) continue;

        // Determine color based on mode
        float r = 1.0f, g = 1.0f, b = 1.0f;
        if (color_mode_ == ColorMode::BY_TERMINATION) {
            get_termination_color(geo.termination, r, g, b);
        }

        // Generate line segments for this ray
        for (size_t i = 0; i < geo.points.size() - 1; ++i) {
            const auto& p1 = geo.points[i];
            const auto& p2 = geo.points[i+1];

            if (color_mode_ == ColorMode::BY_ERROR) {
                get_error_color(p1.H_error, r, g, b);
            }

            float x1, y1, z1;
            to_cartesian(p1.x[Physics::R], p1.x[Physics::THETA], p1.x[Physics::PHI], x1, y1, z1);
            Vertex v1 = {x1, y1, z1, r, g, b, 1.0f};

            if (color_mode_ == ColorMode::BY_ERROR) {
                get_error_color(p2.H_error, r, g, b);
            }

            float x2, y2, z2;
            to_cartesian(p2.x[Physics::R], p2.x[Physics::THETA], p2.x[Physics::PHI], x2, y2, z2);
            Vertex v2 = {x2, y2, z2, r, g, b, 1.0f};

            geodesic_vertices_.push_back(v1);
            geodesic_vertices_.push_back(v2);
        }
    }
    
    if (!vbo_geodesics_) glGenBuffers(1, &vbo_geodesics_);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
    glBufferData(GL_ARRAY_BUFFER,
                 geodesic_vertices_.size() * sizeof(Vertex),
                 geodesic_vertices_.data(),
                 GL_STATIC_DRAW);
}

void Renderer::render(int width, int height, const float view_matrix[16], const float proj_matrix[16]) {
    // Dark deep purple/blue background - lighter than before to see black hole
    glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glUseProgram(shader_program_);
    glBindVertexArray(vao_);
    
    // Compute MVP
    float mvp[16];
    for (int i = 0; i < 16; ++i) mvp[i] = 0.0f;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                mvp[i*4 + j] += proj_matrix[i*4 + k] * view_matrix[k*4 + j];
            }
        }
    }
    
    glUniformMatrix4fv(u_mvp_, 1, GL_FALSE, mvp);
    
    // 1. Draw Starfield (Background)
    if (show_starfield_) {
        draw_points(starfield_vertices_, vbo_starfield_);
    }
    
    // 2. Draw Accretion Disk (blended)
    if (show_accretion_disk_) {
        glDepthMask(GL_FALSE); // Don't write depth, so disk is transparent
        draw_triangles(accretion_disk_vertices_, vbo_accretion_disk_);
        glDepthMask(GL_TRUE);
    }
    
    // 3. Draw Event Horizon (solid black sphere)
    if (show_horizon_) {
        draw_triangles(horizon_vertices_, vbo_horizon_);
    }
    
    // 4. Draw Geodesic Rays (Batched)
    if (!geodesic_vertices_.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        glDrawArrays(GL_LINES, 0, geodesic_vertices_.size());
    }
    
    // 5. Draw Photon Sphere (wireframe overlay)
    if (show_photon_sphere_) {
        draw_lines(photon_sphere_vertices_, vbo_photon_sphere_);
    }
}

void Renderer::draw_lines(const std::vector<Vertex>& vertices, GLuint vbo) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
    glDrawArrays(GL_LINE_STRIP, 0, vertices.size());
}

void Renderer::draw_triangles(const std::vector<Vertex>& vertices, GLuint vbo) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
    glDrawArrays(GL_TRIANGLES, 0, vertices.size());
}

void Renderer::draw_points(const std::vector<Vertex>& vertices, GLuint vbo) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
    glDrawArrays(GL_POINTS, 0, vertices.size());
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

void Renderer::get_termination_color(Numerics::TerminationReason reason, float& r, float& g, float& b) const {
    using namespace Numerics;
    switch(reason) {
        case TerminationReason::HORIZON_CROSSED:
            r = 1.0f; g = 0.1f; b = 0.1f; // Red glow
            break;
        case TerminationReason::ESCAPED:
            r = 0.4f; g = 0.8f; b = 1.0f; // Cyan/Blue glow
            break;
        default:
            r = 0.5f; g = 0.5f; b = 0.5f;
            break;
    }
}

} // namespace Render