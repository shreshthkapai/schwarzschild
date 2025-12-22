#include <emscripten.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <list>
#include <map>
#include <sstream>
#include <iomanip>

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

// Encapsulated application state
struct AppState {
    Render::Camera camera;
    std::unique_ptr<Render::Renderer> renderer;
    App::SimulationParams params;
    App::Controls controls;
    std::vector<Numerics::Geodesic> geodesics;
    
    AppState() : renderer(nullptr) {}
};

// Singleton instance managed by unique_ptr
static std::unique_ptr<AppState> g_app;

// Task 9: LRU Cache for Geodesic Library
class GeodesicCache {
public:
    static const size_t MAX_SIZE = 1000; // Limit to prevent memory bloat

    struct CacheEntry {
        std::string key;
        Numerics::Geodesic geodesic;
    };

    bool get(const std::string& key, Numerics::Geodesic& out_geo) {
        auto it = map_.find(key);
        if (it == map_.end()) return false;
        
        // Move to front (MRU)
        list_.splice(list_.begin(), list_, it->second);
        out_geo = it->second->geodesic;
        return true;
    }

    void put(const std::string& key, const Numerics::Geodesic& geo) {
        // Check if exists
        auto it = map_.find(key);
        if (it != map_.end()) {
            list_.splice(list_.begin(), list_, it->second);
            it->second->geodesic = geo;
            return;
        }

        // Evict if full
        if (list_.size() >= MAX_SIZE) {
            auto last = list_.end();
            last--;
            map_.erase(last->key);
            list_.pop_back();
        }

        // Insert
        list_.push_front({key, geo});
        map_[key] = list_.begin();
    }
    
    // Key generation from RayState
    static std::string make_key(const Rays::RayState& ray) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2);
        // Key by initial position and momentum
        // x: r, theta, phi (t is ignored)
        // p: p_r, p_theta, p_phi (p_t is determined by conservation)
        ss << "r:" << ray.x[Physics::R] 
           << "_th:" << ray.x[Physics::THETA] 
           << "_ph:" << ray.x[Physics::PHI]
           << "_pr:" << std::setprecision(3) << ray.p[Physics::R]
           << "_pth:" << ray.p[Physics::THETA]
           << "_pph:" << ray.p[Physics::PHI];
        return ss.str();
    }
    
private:
    std::list<CacheEntry> list_;
    std::map<std::string, std::list<CacheEntry>::iterator> map_;
};

static GeodesicCache g_cache;

// Forward declarations
void setup_geodesics();
void refire_rays();

// Mouse callbacks
EM_BOOL mouse_callback(int eventType, const EmscriptenMouseEvent* e, void* userData) {
    if (!g_app) return false;
    
    if (eventType == EMSCRIPTEN_EVENT_MOUSEDOWN) {
        g_app->controls.on_mouse_down(e->clientX, e->clientY);
    } else if (eventType == EMSCRIPTEN_EVENT_MOUSEUP) {
        g_app->controls.on_mouse_up();
    } else if (eventType == EMSCRIPTEN_EVENT_MOUSEMOVE) {
        g_app->controls.on_mouse_move(e->clientX, e->clientY);
    }
    return true;
}

EM_BOOL wheel_callback(int eventType, const EmscriptenWheelEvent* e, void* userData) {
    if (!g_app) return false;
    g_app->controls.on_wheel(e->deltaY);
    return true;
}

EM_BOOL key_callback(int eventType, const EmscriptenKeyboardEvent* e, void* userData) {
    if (!g_app) return false;
    if (eventType == EMSCRIPTEN_EVENT_KEYDOWN) {
        g_app->controls.on_key_down(e->key);
    }
    return true;
}

void render_frame() {
    if (!g_app || !g_app->renderer) return;

    int width, height;
    emscripten_get_canvas_element_size("#canvas", &width, &height);
    
    glViewport(0, 0, width, height);
    glClearColor(0.0f, 0.0f, 0.05f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    float view[16];
    float proj[16];
    float aspect = float(width) / float(height);
    
    g_app->camera.compute_view_matrix(view);
    g_app->camera.compute_proj_matrix(proj, aspect);
    
    g_app->renderer->render(width, height, view, proj);
}

void setup_geodesics() {
    using namespace Physics;
    using namespace Rays;
    using namespace Numerics;
    
    if (!g_app) return;
    
    std::cout << "\n=== Setting up geodesics ===" << std::endl;
    std::cout << "Observer r: " << g_app->params.observer_r << std::endl;
    std::cout << "Impact range: [" << g_app->params.impact_min << ", " << g_app->params.impact_max << "]" << std::endl;
    std::cout << "Mode: " << (g_app->params.use_spherical ? "3D Spherical" : "2D Equatorial") << std::endl;
    
    SchwarzschildMetric metric;
    Hamiltonian ham(&metric);
    RayInitializer ray_init(&metric);
    GeodesicIntegrator integrator(&ham);
    
    RayBundle bundle(&ray_init);
    
    if (g_app->params.use_spherical) {
        // 3D spherical ray bundle - rays from all angles
        bundle.generate_spherical_bundle(g_app->params.observer_r, 
                                         g_app->params.impact_min, 
                                         g_app->params.impact_max,
                                         g_app->params.num_theta,
                                         g_app->params.num_phi,
                                         g_app->params.num_impact);
    } else {
        // 2D equatorial bundle (original behavior)
        bundle.generate_uniform_bundle(g_app->params.observer_r, 
                                       g_app->params.impact_min, 
                                       g_app->params.impact_max, 
                                       g_app->params.num_rays_2d);
    }
    
    std::cout << "Integrating " << bundle.size() << " geodesics..." << std::endl;
    
    integrator.set_store_interval(Physics::STORE_INTERVAL);  // Store more points for smoother curves
    g_app->geodesics.clear();
    g_app->geodesics.reserve(bundle.size());
    
    int cache_hits = 0;
    
    for (const auto& ray : bundle.get_rays()) {
        std::string key = GeodesicCache::make_key(ray);
        Numerics::Geodesic geo;
        
        if (g_cache.get(key, geo)) {
            cache_hits++;
        } else {
            geo = integrator.integrate(ray.x, ray.p, 
                                       g_app->params.lambda_step, g_app->params.lambda_max);
            g_cache.put(key, geo);
        }
        
        g_app->geodesics.push_back(geo);
    }
    
    if (cache_hits > 0) {
        std::cout << "Cache hits: " << cache_hits << " / " << bundle.size() << std::endl;
    }
    
    // Statistics
    int captured = 0, escaped = 0, other = 0;
    for (const auto& geo : g_app->geodesics) {
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
    if (g_app && g_app->renderer) {
        g_app->renderer->update_geodesics(g_app->geodesics);
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
    
    // Initialize AppState
    g_app = std::make_unique<AppState>();
    
    // Initialize WebGL
    EmscriptenWebGLContextAttributes attrs;
    emscripten_webgl_init_context_attributes(&attrs);
    attrs.majorVersion = 2;
    attrs.minorVersion = 0;
    
    EMSCRIPTEN_WEBGL_CONTEXT_HANDLE ctx = emscripten_webgl_create_context("#canvas", &attrs);
    if (ctx <= 0) {
        std::cerr << "WebGL context creation failed!" << std::endl;
        emscripten_run_script("alert('WebGL context creation failed! Your browser may not support WebGL 2.0.')");
        return 1;
    }
    emscripten_webgl_make_context_current(ctx);
    
    std::cout << "WebGL context created" << std::endl;
    
    // Initialize renderer
    g_app->renderer = std::make_unique<Render::Renderer>();
    if (!g_app->renderer->initialize()) {
        std::cerr << "Renderer initialization failed!" << std::endl;
        emscripten_run_script("alert('Renderer initialization failed! Check console for shader errors.')");
        return 1;
    }
    std::cout << "Renderer initialized" << std::endl;
    
    // Setup controls
    g_app->controls.set_camera(&g_app->camera);
    g_app->controls.set_renderer(g_app->renderer.get());
    g_app->controls.set_params(&g_app->params);
    g_app->controls.set_refire_callback(refire_rays);
    
    // Setup initial geodesics
    setup_geodesics();
    g_app->renderer->update_geodesics(g_app->geodesics);
    
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
    if (!g_app) return;
    
    g_app->params.observer_r = observer_r;
    g_app->params.impact_min = impact_min;
    g_app->params.impact_max = impact_max;
    g_app->params.num_theta = num_theta;
    g_app->params.num_phi = num_phi;
    g_app->params.use_spherical = use_spherical;
    
    // Auto-refire when params change
    refire_rays();
}

void web_set_toggles(bool horizon, bool photon, bool disk, bool stars) {
    if (g_app && g_app->renderer) {
        g_app->renderer->set_show_horizon(horizon);
        g_app->renderer->set_show_photon_sphere(photon);
        g_app->renderer->set_show_accretion_disk(disk);
        g_app->renderer->set_show_starfield(stars);
        // Note: Color mode cycle is handled separately or can be added here
    }
}

void web_update_geodesics_from_buffer(const std::vector<float>& data) {
    if (g_app && g_app->renderer) {
        g_app->renderer->update_geodesics_from_buffer(data);
    }
}

void web_clear_geodesics() {
    if (g_app && g_app->renderer) {
        g_app->renderer->clear_geodesics();
    }
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("updateParams", &web_update_params);
    function("setToggles", &web_set_toggles);
    function("updateGeodesicsFromBuffer", &web_update_geodesics_from_buffer);
    function("clearGeodesics", &web_clear_geodesics);
}