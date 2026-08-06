#include "model.hpp"
#include "renderer.hpp"
#include "utils.hpp"
#include <chrono>
#include <thread>

const double SMOOTHING_LENGTH = 1000000 * 0.3;
const double DENSITY_IN_REST  = 1e20 * 5;
const double TEMPERATURE_K    = 15.0;
const double VISCOSITY_K      = 1e19;

const double SHEAR_MODULUS = 1e28 * 2;
const double BULK_MODULUS = 1e29 * 3;

const double PARTICLE_RADIUS  = 10000;
const double BARNES_HUT_THETA = 0.5;

const double TIME_DELTA = 3;

const double SCALE = MOON_EARTH_DIST / 20;

// const double WATER_LAYER_THICKNESS  = 1000000; // м
const int    WATER_PARTICLES_NUMBER = 500;

int main() {
    LiquidsSettings l_settings {
        SMOOTHING_LENGTH,
        DENSITY_IN_REST,
        TEMPERATURE_K * 1.5,
        VISCOSITY_K * 0.02
    };

    ElasticSolidsSettings e_settings {
        SMOOTHING_LENGTH * 200,
        SHEAR_MODULUS,
        BULK_MODULUS
    };

    GravitySettings g_settings {
        PARTICLE_RADIUS,
        BARNES_HUT_THETA,
        INF
    };

    ModelSettings settings {
        l_settings,
        e_settings,
        g_settings,
        TIME_DELTA,
        SCALE
    };

    Model model(settings); // TODO: scale

    Planet Earth {
        glm::dvec2 { 0, 0 },
        EARTH_RADIUS,
        EARTH_MASS,
        glm::dvec2 { 0, 0 },
    };

    Planet Moon {
        glm::dvec2 { MOON_EARTH_DIST / 20, 0 },
        MOON_RADIUS,
        MOON_MASS,
        glm::dvec2 { 0, calc_orbital_velocity(EARTH_MASS, MOON_EARTH_DIST / 15) },
        // glm::dvec2 { 0, 0 },
        // MOON_ANGULAR_VELOCITY * 200 
    };

    model.add_planet("Земля", Earth);
    // model.add_particle_planet(Moon, WATER_PARTICLES_NUMBER);
    model.add_particle_satellite("Земля", Moon, WATER_PARTICLES_NUMBER, calc_orbital_velocity(EARTH_MASS, MOON_EARTH_DIST / 20));

    Renderer renderer;

    renderer.init();

    while (true) {
        int step = 40;

        auto start = std::chrono::high_resolution_clock::now();
        auto end = start;

        int cnt = 0;
        while (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() < step) {
            model.next_step();
            end = std::chrono::high_resolution_clock::now();
            cnt++;
        }

        std::cout << "steps in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms: " << cnt << '\n';

        renderer.tick(model.get_frame());

        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
}