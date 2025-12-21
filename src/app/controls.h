#ifndef CONTROLS_H
#define CONTROLS_H

#include "app/simulation_params.h"
#include "render/renderer.h"
#include "render/camera.h"
#include <functional>

namespace App {

// Keyboard/mouse input handler
class Controls {
public:
    Controls();
    
    // Set references to controllable objects
    void set_camera(Render::Camera* cam) { camera_ = cam; }
    void set_renderer(Render::Renderer* rend) { renderer_ = rend; }
    void set_params(SimulationParams* params) { params_ = params; }
    
    // Set callback for refiring rays
    void set_refire_callback(std::function<void()> callback) { refire_callback_ = callback; }
    
    // Handle keyboard input (returns true if handled)
    bool on_key_down(const char* key);
    
    // Handle mouse input
    void on_mouse_down(double x, double y);
    void on_mouse_up();
    void on_mouse_move(double x, double y);
    void on_wheel(double delta_y);
    
    // Get current states
    bool is_dragging() const { return is_dragging_; }
    
    // Print current parameters to console
    void print_params() const;
    
private:
    Render::Camera* camera_;
    Render::Renderer* renderer_;
    SimulationParams* params_;
    std::function<void()> refire_callback_;
    
    // Mouse state
    bool is_dragging_;
    double last_mouse_x_;
    double last_mouse_y_;
    
    // Toggle helpers
    void toggle_horizon();
    void toggle_photon_sphere();
    void cycle_color_mode();
};

} // namespace App

#endif // CONTROLS_H
