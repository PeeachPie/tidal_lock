#include "model.hpp"
#include "renderer.hpp"
#include "barnes_hut.hpp"
#include <chrono>
#include <thread>

int main() {
    Model model(SMOOTHING_LENGTH, TIME_DELTA, MOON_EARTH_DIST / 20); // TODO: scale

    Planet Earth {
        glm::dvec2 { 0, 0 },
        EARTH_RADIUS,
        EARTH_MASS,
        glm::dvec2 { 0, -3222 },
    };

    Planet Moon {
        glm::dvec2 { MOON_EARTH_DIST / 20, 0 },
        MOON_RADIUS,
        EARTH_MASS,
        glm::dvec2 { 0, 3222 },
        // glm::dvec2 { 0, 0 },
        // MOON_ANGULAR_VELOCITY * 300
    };

    model.add_planet("Земля", Earth);
    model.add_planet("Луна", Moon);

    model.add_water_to_planet("Земля", WATER_PARTICLES_NUMBER, WATER_LAYER_THICKNESS);
    // model.add_water_to_planet("Луна", WATER_PARTICLES_NUMBER, MOON_RADIUS);

    Renderer renderer;

    renderer.init();

    while (true) {

        // auto start = std::chrono::high_resolution_clock::now();
        // for (int i = 0; i < 10; i++) {
            model.next_step();
        // }
        // auto end = std::chrono::high_resolution_clock::now();

        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "next_step: " << duration.count() << " ms\n";

        renderer.tick(model.get_frame());

        // std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    }

    // BarnesHutQuadTree tree(20);

    // Particle a({10, 10});
    // Particle b({5, 5});
    // Particle c({0, 0});
    // Particle d({-5, -5});
    // Particle e({-10, -10});
    // Particle f({-5, 5});

    // tree.insert(a, 1, 1);
    // tree.insert(b, 1, 2);
    // tree.insert(c, 1, 3);
    // tree.insert(d, 1, 4);
    // tree.insert(e, 1, 5);
    // tree.insert(f, 1, 6);

    // glm::dvec2 fs = tree.calc_force(c, 1, 3);

    // std::cout << fs.x << ' ' << fs.y << '\n';

    // tree.print();
}