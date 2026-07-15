#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <glm/glm.hpp>
#include "utils.hpp"
#include "constants.hpp"
#include "barnes_hut.hpp"
#include <iostream>
#include "particle.hpp"

class Model {
private:
    // double h_;
    double time_step_;
    double current_time = 0;

    std::map<std::string, Planet> planets_;
    std::vector<Particle> p_;

    // TODO: as param
    double p_mass_ = PARTICLE_MASS;
    double p_smoothing_lenght_ = SMOOTHING_LENGTH;
    double p_density_in_rest = DENSITY_IN_REST;
    double p_temperature_k = TEMPERATURE_K;
    double p_r_ = 0.05 * 40;

    double scale_;

public:
    Model(double h, double time_delta, double scale);

    void add_planet(std::string name, Planet planet);

    void add_water_to_planet(std::string planet_name, int p_count, double h);

    Frame get_frame();

    void next_step();

    double calc_energy() const;

private:
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> _get_particle_grid();
    std::vector<int> _get_particle_neighborhood(int pi, const std::unordered_map<int, std::unordered_map<int, std::vector<int>>> &particle_grid);
    std::vector<std::vector<int>> _get_neighborhoods();

    void _planets_next_step();
    void _particles_next_step();
    
    glm::dvec2 _particle_gravity_forces(int i, BarnesHutQuadTree &tree);
    glm::dvec2 _planet_gravity_forces(std::string planet_name);
    double _sph_calc_particle_density(int i, const std::vector<int> &neighborhood);
    std::vector<double> _sph_calc_density(const std::vector<std::vector<int>> &neighborhoods);
    std::vector<double> _sph_calc_pressure(std::vector<double> density, const std::vector<std::vector<int>> &neighborhoods);
    glm::dvec2 _particle_pressure_forces(int i, const std::vector<int> &neighborhood, const std::vector<double> &density, const std::vector<double> pressure);
    glm::dvec2 _particle_viscosity_forces(int i, const std::vector<int> &neighborhood, const std::vector<double> &density);
    
    double _calc_density_in_rest();

    double _sph_w_density(double dist);
    double _sph_w_pressure_grad(double dist);
    double _sph_w_viscosity_laplacian(double dist);
    
    void add_particle(Particle particle);
};