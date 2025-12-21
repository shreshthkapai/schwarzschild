#include <emscripten.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>
#include <cmath>
#include <iostream>

#include "physics/constants.h"
#include "physics/schwarzschild_metric.h"
#include "physics/hamiltonian.h"
#include "numerics/integrator.h"
#include "rays/ray_initializer.h"
#include "rays/ray_bundle.h"
#include "render/renderer.h"
#include "render/camera.h"
#include "app/simulation_params.h"
#include "app/controls.h"

// Global state
Render::Camera camera;
Render::Renderer* renderer = nullptr;
App::SimulationParams params;
App::Controls controls;
std::vector<Numerics::Geodesic> geodesics;

// Forward declarations
void setup_geodesics();
void refire_rays();

// Mouse callbacks
EM_BOOL mouse_callback(int eventType, const EmscriptenMouseEvent* e, void* userData) {
    if (eventType == EMSCRIPTEN_EVENT_MOUSEDOWN) {
        controls.on_mouse_down(e->clientX, e->clientY);
    } else if (eventType == EMSCRIPTEN_EVENT_MOUSEUP) {
        controls.on_mouse_up();
    } else if (eventType == EMSCRIPTEN_EVENT_MOUSEMOVE) {
        controls.on_mouse_move(e->clientX, e->clientY);
    }
    return true;
}

EM_BOOL wheel_callback(int eventType, const EmscriptenWheelEvent* e, void* userData) {
    controls.on_wheel(e->deltaY);
    return true;
}

EM_BOOL key_callback(int eventType, const EmscriptenKeyboardEvent* e, void* userData) {
    if (eventType == EMSCRIPTEN_EVENT_KEYDOWN) {
        controls.on_key_down(e->key);
    }
    return true;
}

void render_frame() {
    int width, height;
    emscripten_get_canvas_element_size("#canvas", &width, &height);
    
    glViewport(0, 0, width, height);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    
    float view[16];
    float proj[16];
    float aspect = float(width) / float(height);
    
    camera.compute_view_matrix(view);
    camera.compute_proj_matrix(proj, aspect);
    
    if (renderer) {
        renderer->render(width, height, view, proj);
    }
}

void setup_geodesics() {
    using namespace Physics;
    using namespace Rays;
    using namespace Numerics;
    
    std::cout << "\n=== Setting up geodesics ===" << std::endl;
    std::cout << "Observer r: " << params.observer_r << std::endl;
    std::cout << "Impact range: [" << params.impact_min << ", " << params.impact_max << "]" << std::endl;
    std::cout << "Mode: " << (params.use_spherical ? "3D Spherical" : "2D Equatorial") << std::endl;
    
    SchwarzschildMetric metric;
    Hamiltonian ham(&metric);
    RayInitializer ray_init(&metric);
    GeodesicIntegrator integrator(&ham);
    
    RayBundle bundle(&ray_init);
    
    if (params.use_spherical) {
        // 3D spherical ray bundle - rays from all angles
        bundle.generate_spherical_bundle(params.observer_r, 
                                         params.impact_min, 
                                         params.impact_max,
                                         params.num_theta,
                                         params.num_phi,
                                         params.num_impact);
    } else {
        // 2D equatorial bundle (original behavior)
        bundle.generate_uniform_bundle(params.observer_r, 
                                       params.impact_min, 
                                       params.impact_max, 
                                       params.num_rays_2d);
    }
    
    std::cout << "Integrating " << bundle.size() << " geodesics..." << std::endl;
    
    integrator.set_store_interval(15);  // Store more points for smoother curves
    geodesics.clear();
    geodesics.reserve(bundle.size());
    
    for (const auto& ray : bundle.get_rays()) {
        Geodesic geo = integrator.integrate(ray.x, ray.p, 
                                            params.lambda_step, params.lambda_max);
        geodesics.push_back(geo);
    }
    
    // Statistics
    int captured = 0, escaped = 0, other = 0;
    for (const auto& geo : geodesics) {
        if (geo.termination == TerminationReason::HORIZON_CROSSED) captured++;
        else if (geo.termination == TerminationReason::ESCAPED) escaped++;
        else other++;
    }
    
    std::cout << "Results: " << captured << " captured, " 
              << escaped << " escaped, " << other << " other" << std::endl;
    std::cout << "===========================" << std::endl;
}

void refire_rays() {
    setup_geodesics();
    if (renderer) {
        renderer->set_geodesics(geodesics);
    }
}

int main() {
    std::cout << "=== Schwarzschild Geodesics Visualization ===" << std::endl;
    std::cout << "Event horizon: r = " << Physics::R_SCHWARZSCHILD << std::endl;
    std::cout << "Photon sphere: r = " << Physics::R_PHOTON_SPHERE << std::endl;
    std::cout << "\n--- Controls ---" << std::endl;
    std::cout << "Mouse drag: Rotate camera" << std::endl;
    std::cout << "Scroll: Zoom" << std::endl;
    std::cout << "Arrow Up/Down: Adjust observer distance" << std::endl;
    std::cout << "Arrow Left/Right: Adjust number of rays" << std::endl;
    std::cout << "[ / ]: Adjust impact parameter range" << std::endl;
    std::cout << "R: Refire rays" << std::endl;
    std::cout << "H: Toggle horizon" << std::endl;
    std::cout << "P: Toggle photon sphere" << std::endl;
    std::cout << "C: Cycle color mode" << std::endl;
    std::cout << "I: Print current parameters" << std::endl;
    std::cout << "----------------\n" << std::endl;
    
    // Initialize WebGL
    EmscriptenWebGLContextAttributes attrs;
    emscripten_webgl_init_context_attributes(&attrs);
    attrs.majorVersion = 2;
    attrs.minorVersion = 0;
    
    EMSCRIPTEN_WEBGL_CONTEXT_HANDLE ctx = emscripten_webgl_create_context("#canvas", &attrs);
    emscripten_webgl_make_context_current(ctx);
    
    std::cout << "WebGL context created" << std::endl;
    
    // Initialize renderer
    renderer = new Render::Renderer();
    renderer->initialize();
    std::cout << "Renderer initialized" << std::endl;
    
    // Setup controls
    controls.set_camera(&camera);
    controls.set_renderer(renderer);
    controls.set_params(&params);
    controls.set_refire_callback(refire_rays);
    
    // Setup initial geodesics
    setup_geodesics();
    renderer->set_geodesics(geodesics);
    
    // Input callbacks
    emscripten_set_mousedown_callback("#canvas", nullptr, true, mouse_callback);
    emscripten_set_mouseup_callback("#canvas", nullptr, true, mouse_callback);
    emscripten_set_mousemove_callback("#canvas", nullptr, true, mouse_callback);
    emscripten_set_wheel_callback("#canvas", nullptr, true, wheel_callback);
    emscripten_set_keydown_callback(EMSCRIPTEN_EVENT_TARGET_DOCUMENT, nullptr, true, key_callback);
    
    std::cout << "Ready! Use keyboard controls listed above." << std::endl;
    
    // Start render loop
    emscripten_set_main_loop(render_frame, 0, true);
    
    return 0;
}

#include <emscripten/bind.h>

using namespace emscripten;

// Web interface wrappers
void web_update_params(double observer_r, double impact_min, double impact_max, 
                      int num_theta, int num_phi, bool use_spherical) {
    params.observer_r = observer_r;
    params.impact_min = impact_min;
    params.impact_max = impact_max;
    params.num_theta = num_theta;
    params.num_phi = num_phi;
    params.use_spherical = use_spherical;
    
    // Auto-refire when params change
    refire_rays();
}

void web_set_toggles(bool horizon, bool photon, bool disk, bool stars) {
    if (renderer) {
        renderer->set_show_horizon(horizon);
        renderer->set_show_photon_sphere(photon);
        renderer->set_show_accretion_disk(disk);
        renderer->set_show_starfield(stars);
        // Note: Color mode cycle is handled separately or can be added here
    }
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("updateParams", &web_update_params);
    function("setToggles", &web_set_toggles);
}