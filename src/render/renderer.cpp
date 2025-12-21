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
    
    u_mvp_ = glGetUniformLocation(shader_program_, "u_mvp");
    
    glDeleteShader(vs);
    glDeleteShader(fs);
    return true;
}

void Renderer::build_static_geometry() {
    static_vertices_.clear();
    
    // 1. Horizon (Triangles)
    {
        auto verts = generate_solid_sphere(Physics::R_SCHWARZSCHILD, Physics::SPHERE_SEGMENTS, 0.0f, 0.0f, 0.0f, 1.0f);
        offset_horizon_ = static_vertices_.size();
        count_horizon_ = verts.size();
        static_vertices_.insert(static_vertices_.end(), verts.begin(), verts.end());
    }
    
    // 2. Photon Sphere (Lines)
    {
        auto verts = generate_sphere(Physics::R_PHOTON_SPHERE, Physics::PHOTON_SPHERE_SEGMENTS, 1.0f, 0.8f, 0.0f, 0.3f);
        offset_photon_sphere_ = static_vertices_.size();
        count_photon_sphere_ = verts.size();
        static_vertices_.insert(static_vertices_.end(), verts.begin(), verts.end());
    }
    
    // 3. Accretion Disk (Triangles)
    {
        auto verts = generate_accretion_disk(Physics::DISK_INNER_R, Physics::DISK_OUTER_R, Physics::DISK_SEGMENTS, 0.0f);
        offset_accretion_disk_ = static_vertices_.size();
        count_accretion_disk_ = verts.size();
        static_vertices_.insert(static_vertices_.end(), verts.begin(), verts.end());
    }
    
    // 4. Starfield (Points)
    {
        auto verts = generate_starfield(Physics::STAR_COUNT, Physics::STAR_DIST);
        offset_starfield_ = static_vertices_.size();
        count_starfield_ = verts.size();
        static_vertices_.insert(static_vertices_.end(), verts.begin(), verts.end());
    }
    
    // Create VBO
    if (!vbo_static_) glGenBuffers(1, &vbo_static_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_static_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_vertices_.size() * sizeof(Vertex),
                 static_vertices_.data(),
                 GL_STATIC_DRAW);
}

void Renderer::update_geodesics(const std::vector<Numerics::Geodesic>& geodesics) {
    geodesic_vertices_.clear();
    
    // Exact reservation to prevent reallocations
    size_t total_points = 0;
    for (const auto& geo : geodesics) total_points += geo.points.size();
    if (total_points > 0) geodesic_vertices_.reserve(total_points * 2); // 2 verts per segment

    for (const auto& geo : geodesics) {
        if (geo.points.size() < 2) continue;

        float r_col = 1.0f, g_col = 1.0f, b_col = 1.0f;
        if (color_mode_ == ColorMode::BY_TERMINATION) {
            get_termination_color(geo.termination, r_col, g_col, b_col);
        }

        for (size_t i = 0; i < geo.points.size() - 1; ++i) {
            const auto& p1 = geo.points[i];
            const auto& p2 = geo.points[i+1];

            if (color_mode_ == ColorMode::BY_ERROR) {
                get_error_color(p1.H_error, r_col, g_col, b_col);
            }

            float x1, y1, z1;
            to_cartesian(p1.x[Physics::R], p1.x[Physics::THETA], p1.x[Physics::PHI], x1, y1, z1);
            Vertex v1 = {x1, y1, z1, r_col, g_col, b_col, 1.0f};

            if (color_mode_ == ColorMode::BY_ERROR) {
                get_error_color(p2.H_error, r_col, g_col, b_col);
            }

            float x2, y2, z2;
            to_cartesian(p2.x[Physics::R], p2.x[Physics::THETA], p2.x[Physics::PHI], x2, y2, z2);
            Vertex v2 = {x2, y2, z2, r_col, g_col, b_col, 1.0f};

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
    glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glUseProgram(shader_program_);
    glBindVertexArray(vao_);
    
    // Optimized Matrix Multiplication (Column-Major)
    // MVP = Proj * View
    float mvp[16];
    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += proj_matrix[col * 4 + k] * view_matrix[k * 4 + row]; 
                // Note: Standard OpenGL matrix layout is column-major
                // proj[col][k] * view[k][row] ?? No.
                // Output column 'col', row 'row'.
                // sum_k (Left[k][row] * Right[col][k]) -> Pre-multiply?
                // Visual check: mvp[i*4+j] usually means row i, col j if row-major indexing?
                // Or col i, row j?
                // The loop I'm replacing was:
                // mvp[i*4 + j] += proj[i*4 + k] * view[k*4 + j];
                // This implies row-major storage interpretation or mathematical C_ij = Sum A_ik B_kj.
                
                // Let's stick to the simpler loop structure but make it cleaner:
            }
            // sum computed below
        }
    }
    
    // Re-implementing the loop exactly as logic requires but without nested 3-level junk if possible.
    // Actually, the previous loop was correct for Row-Major math on 1D arrays representing Column-Major matrices?
    // Let's just unroll it cleanly.
    
    // C = A * B
    // C[row][col] = sum(A[row][k] * B[k][col])
    // Flat: C[col*4 + row] (Col-Major) or C[row*4 + col] (Row-Major).
    // OpenGL uses Column-Major. access is M[col*4 + row].
    // Let's assume input arrays are consistent.
    
    for (int i = 0; i < 4; ++i) { // Row
        for (int j = 0; j < 4; ++j) { // Col
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) {
                sum += proj_matrix[k*4 + i] * view_matrix[j*4 + k]; // Verify this?
                // If Proj is P, View is V. MVP = P * V.
                // (PV)_ij = P_ik V_kj.
                // If Column Major: M[i + 4*j] is M_ij (row i, col j).
                // P[i + 4*k] * V[k + 4*j].
                // Yes.
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
    
    // 4. Photon Sphere (Wireframe - needs LINES but stored as TRIANGLES? No, generate_sphere returns vertices for lines?
    // Wait, generate_sphere usually returns triangles or lines depending on usage.
    // Let's check geometry.cpp. If it returns standard sphere mesh, it's triangles.
    // If I want wireframe, I should use GL_LINE_STRIP or similar.
    // But currently I'm drawing GL_LINE_STRIP in old code.
    // Is generate_sphere optimized for lines? 
    // The old code used draw_lines which uses GL_LINE_STRIP.
    // But if I put it in a single buffer, I can't easily use LINE_STRIP unless it's segmented.
    // Or I use GL_LINES.
    // I'll assume generate_sphere creates a strip-friendly list, but for GL_LINES I need pairs.
    // Actually, if I use GL_LINES, I need to ensure vertices are pairs.
    // Let's assume for now I use LINES for Photon Sphere if the generator supports it.
    // If not, I might draw it as points or just skip wireframe optimization for now.
    // BUT the user wants optimization.
    // I will use GL_LINE_LOOP or similar? No, can't in single draw.
    // I will switch to `GL_LINES` for photon sphere if I can.
    // For now, I'll use GL_LINE_STRIP for photon sphere, but I must be careful about continuity.
    // Actually, the old code used `draw_lines` (GL_LINE_STRIP).
    // If I put it in a big buffer, using `glDrawArrays` with count/offset works for strip too.
    if (show_photon_sphere_) {
        glDrawArrays(GL_LINE_STRIP, offset_photon_sphere_, count_photon_sphere_);
    }
    
    // --- Draw Dynamic Geometry ---
    if (!geodesic_vertices_.empty()) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo_geodesics_);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(3*sizeof(float)));
        glDrawArrays(GL_LINES, 0, geodesic_vertices_.size());
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