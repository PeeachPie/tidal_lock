#include "model.hpp"
#include "renderer.hpp"
#include "utils.hpp"
#include <chrono>
#include <thread>

const double PARTICLE_DIST    = 4500;  // в метрах. Выбор зависит от радиуса планеты и производительности системы
const double SMOOTHING_LENGTH = 15000; // в метрах. Следует выбирать так, чтобы кол-во соседей было около 30

const double DENSITY_IN_REST  = DESPINA_MASS / (PI * DESPINA_RADIUS * DESPINA_RADIUS);
const double TEMPERATURE_K    = 2000.0;   // коэф давления, необходим для больших планет, чтобы частицы не слипались под действием гравитации 
const double VISCOSITY_K      = 1e13 * 6; // коэф диссипативных сил - влияет на скорость приливного захвата

const double SHEAR_MODULUS = 50e9; // сопротивление сдвигам
const double BULK_MODULUS  = 90e9; // сопротивление измененям объема

const double PARTICLE_RADIUS  = 10000; // порог близости частиц для расчета сил гравитации
const double BARNES_HUT_THETA = 0.5;   // 0.0-1.0 Чем ближе к 0.0, тем дольше, но точнее расчет сил гравитации

const double TIME_DELTA = 1; // шаг симуляции в секундах

const double SCALE           = MOON_EARTH_DIST / 8; // начальное приближение камеры
const double SATELLITE_SCALE = 100;                 // увеличение спутника при отрисовке

const bool DEBUG = true; // логирование в консоль и debug_log.csv
const int STEP = 40;     // минимальный шаг отрисовки в мс

int main() {
    LiquidsSettings l_settings {
        SMOOTHING_LENGTH,
        DENSITY_IN_REST,
        TEMPERATURE_K,
        0.0
    };

    ElasticSolidsSettings e_settings {
        SMOOTHING_LENGTH * 3,
        SHEAR_MODULUS,
        BULK_MODULUS,
        VISCOSITY_K
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
        SCALE,
        SATELLITE_SCALE,
        DEBUG
    };

    Model model(settings); // TODO: scale

    Planet Neptune {
        glm::dvec2 { 0, 0 },
        NEPTUNE_RADIUS,
        NEPTUNE_MASS
    };

    Planet Despina {
        glm::dvec2 { NEPTUNE_DESPINA_DIST, 0 },
        DESPINA_RADIUS,
        DESPINA_MASS,
        glm::dvec2 { 0, calc_orbital_velocity(NEPTUNE_MASS, NEPTUNE_DESPINA_DIST) },
        glm::dvec2 { 0, 0 },
        -MOON_ANGULAR_VELOCITY * 100  
    };

    model.add_planet("Нептун", Neptune);

    // model.add_particle_planet(Despina, PARTICLE_DIST);
    model.add_particle_satellite(
        "Нептун", 
        Despina, 
        PARTICLE_DIST, 
        calc_orbital_velocity(NEPTUNE_MASS, NEPTUNE_DESPINA_DIST)
    );

    Renderer renderer;

    renderer.init();

    while (true) {
        auto start = std::chrono::high_resolution_clock::now();
        auto end = start;

        int cnt = 0;
        while (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() < STEP) {
            model.next_step();
            end = std::chrono::high_resolution_clock::now();
            cnt++;
        }

        std::cout << "steps in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms: " << cnt << '\n';

        renderer.tick(model.get_frame());

        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
}