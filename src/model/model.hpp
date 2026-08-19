#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include <fstream>
#include <cassert>

#include "utils.hpp"
#include "constants.hpp"
#include "particle.hpp"
#include "elastic_solids.hpp"
#include "liquid.hpp"
#include "gravity.hpp"

struct ModelSettings {
    LiquidsSettings liquids_settings;
    ElasticSolidsSettings elastic_solids_settings;
    GravitySettings gravity_settings;

    double time_step;

    double scale;
    double planet_scale;

    bool debug;
};

class Model {
private:
    LiquidsSolver liquids_solver;
    ElasticSolidsSolver elastic_solids_solver;
    GravitySolver gravity_solver;

    // TODO: CollisionsSolver

    double time_step_;

    double scale_;
    double planet_scale_;

    bool debug_;

    double current_time = 0;
    std::map<std::string, Planet> planets_;
    std::vector<Particle> p_;

    double p_r_ = 0.05 * 40;

    std::vector<Particle> initial_p_;
    glm::dvec2 initial_bc_;

    std::ofstream debug_csv_;
    long long step_counter_ = 0;
public:
    Model(ModelSettings &settings);

    void add_planet(std::string name, Planet planet);

    void add_water_to_planet(std::string planet_name, double m, double h, double d);

    void add_particle_planet(Planet planet, double d);

    void add_particle_satellite(std::string planet_name, const Planet &satellite, double d, double velocity);

    Frame get_frame();

    void next_step();

private:
    void _add_particle(Particle particle);
    void _add_water_to_planet(const Planet& planet, double m, double h, double d);

    void _planets_next_step(const std::map<std::string, glm::dvec2> &forces);
    void _particles_next_step(const std::vector<glm::dvec2> &forces);

    glm::dvec2 _get_barycenter();
    std::pair<double, double> _get_angular_velocity(glm::dvec2 satellite_bc, std::string planet_name);
};