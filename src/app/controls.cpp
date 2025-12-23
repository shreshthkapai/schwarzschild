#include "app/controls.h"
#include <iostream>
#include <cstring>
#include <emscripten/html5.h> // For canvas size

namespace App {

Controls::Controls()
    : camera_(nullptr), renderer_(nullptr), params_(nullptr),
      is_dragging_(false), last_mouse_x_(0), last_mouse_y_(0) {}

bool Controls::on_key_down(const char* key) {
    if (!key) return false;
    
    // Debug: print received key
    std::cout << "[Controls] Key pressed: '" << key << "'" << std::endl;
    
    // Toggles
    if (strcmp(key, "h") == 0 || strcmp(key, "H") == 0) {
        toggle_horizon();
        return true;
    }
    if (strcmp(key, "p") == 0 || strcmp(key, "P") == 0) {
        toggle_photon_sphere();
        return true;
    }
    if (strcmp(key, "c") == 0 || strcmp(key, "C") == 0) {
        cycle_color_mode();
        return true;
    }
    
    return false;
}

void Controls::on_mouse_down(double x, double y) {
    is_dragging_ = true;
    last_mouse_x_ = x;
    last_mouse_y_ = y;
}

void Controls::on_mouse_up() {
    is_dragging_ = false;
}

void Controls::on_mouse_move(double x, double y) {
    if (is_dragging_ && camera_) {
        double dx = x - last_mouse_x_;
        double dy = y - last_mouse_y_;
        camera_->rotate(dx * 0.005f, dy * 0.005f);
        last_mouse_x_ = x;
        last_mouse_y_ = y;
    }
}

void Controls::on_wheel(double delta_y) {
    if (camera_) {
        camera_->zoom(1.0 + delta_y * 0.001);
    }
}

void Controls::toggle_horizon() {
    if (renderer_) {
        bool current = renderer_->get_show_horizon();
        renderer_->set_show_horizon(!current);
        std::cout << "[Controls] Horizon: " << (!current ? "ON" : "OFF") << std::endl;
    }
}

void Controls::toggle_photon_sphere() {
    if (renderer_) {
        bool current = renderer_->get_show_photon_sphere();
        renderer_->set_show_photon_sphere(!current);
        std::cout << "[Controls] Photon sphere: " << (!current ? "ON" : "OFF") << std::endl;
    }
}

void Controls::cycle_color_mode() {
    if (renderer_) {
        Render::ColorMode current = renderer_->get_color_mode();
        Render::ColorMode next;
        const char* mode_name;
        
        switch (current) {
            case Render::ColorMode::BY_TERMINATION:
                next = Render::ColorMode::BY_ERROR;
                mode_name = "BY_ERROR";
                break;
            case Render::ColorMode::BY_ERROR:
                next = Render::ColorMode::SOLID;
                mode_name = "SOLID";
                break;
            case Render::ColorMode::SOLID:
                next = Render::ColorMode::DOPPLER; 
                mode_name = "DOPPLER";
                break;
            case Render::ColorMode::DOPPLER:
                next = Render::ColorMode::LENSING; // Task 28
                mode_name = "LENSING";
                break;
            default:
                next = Render::ColorMode::BY_TERMINATION;
                mode_name = "BY_TERMINATION";
                break;
        }
        
        renderer_->set_color_mode(next);
        std::cout << "[Controls] Color mode: " << mode_name << std::endl;
        
        if (refresh_callback_) {
            refresh_callback_();
        }
    }
}

} // namespace App
