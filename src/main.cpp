#include "model.hpp"
#include "renderer.hpp"
#include <chrono>
#include <thread>

int main() {
    Model model(SMOOTHING_LENGTH, TIME_DELTA, MOON_EARTH_DIST / 20); // TODO: scale

    Planet Earth {
        glm::dvec2 { 0, 0 },
        glm::dvec2 { 0, 0 },
        0,
        EARTH_RADIUS,
        EARTH_MASS,
        EARTH_ANGULAR_VELOCITY,
        0
    };

    Planet Moon {
        glm::dvec2 { MOON_EARTH_DIST / 15, 0 },
        glm::dvec2 { 0, 0 },
        MOON_EARTH_DIST / 15,
        MOON_RADIUS,
        EARTH_MASS * 4.593296553894743 * 1.0,
        MOON_ANGULAR_VELOCITY * 150,
        MOON_VELOCITY
    };

    // Planet Moon2 {
    //     glm::dvec2 { 0, MOON_EARTH_DIST / 10 },
    //     glm::dvec2 { 0, 0 },
    //     MOON_EARTH_DIST / 10,
    //     0,
    //     0,
    //     0,
    //     0
    // };

    model.add_planet("Земля", Earth);
    // model.add_planet("Луна", Moon2);
    model.add_planet("Луна", Moon);

    model.add_water_to_planet("Земля", WATER_PARTICLES_NUMBER, WATER_LAYER_THICKNESS);
    // model.add_water_to_planet("Луна2", WATER_PARTICLES_NUMBER / 3, WATER_LAYER_THICKNESS);
    // model.add_water_to_planet("Луна", WATER_PARTICLES_NUMBER / 2, MOON_RADIUS);

    Renderer renderer;

    renderer.init();

    std::cout << model.calc_energy() << '\n';
    int i = 0;

    while (true) {

        // auto start = std::chrono::high_resolution_clock::now();
        // for (int i = 0; i < 10; i++) {
            model.next_step();
        // }
        // auto end = std::chrono::high_resolution_clock::now();

        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "next_step: " << duration.count() << " ms\n";

        renderer.tick(model.get_frame());

        if (i % 100 == 0) {
            std::cout << i << ' ' << model.calc_energy() << '\n';
        }
        i++;
        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
}