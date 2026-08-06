#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>

#include "utils.hpp"
#include "constants.hpp"
#include "particle.hpp"
#include "neighborhood_search.hpp"

struct LiquidsSettings {
    // double p_mass;
    double p_smoothing_length;
    double p_density_in_rest;
    double p_temperature_k;
    double p_viscousity_k;
};

class LiquidsSolver {
private:
    // double p_mass_;
    double p_smoothing_length_;
    double p_density_in_rest_;
    double p_temperature_k_;
    double p_viscousity_k_;

public:
    LiquidsSolver(LiquidsSettings &settings);

    std::vector<glm::dvec2> calc_forces(const std::vector<Particle> &p);

private:
    double _sph_calc_particle_density(
        int i, 
        const std::vector<Particle> &p, 
        const std::vector<int> &neighborhood
    );

    std::vector<double> _sph_calc_density(
        const std::vector<Particle> &p, 
        const std::vector<std::vector<int>> &neighborhoods
    );

    std::vector<double> _sph_calc_pressure(
        const std::vector<Particle> &p, 
        const std::vector<double> &density, 
        const std::vector<std::vector<int>> &neighborhoods
    );

    glm::dvec2 _sph_particle_pressure_force(
        int i, 
        const std::vector<Particle> &p, 
        const std::vector<int> &neighborhood, 
        const std::vector<double> &density, 
        const std::vector<double> &pressure
    );

    glm::dvec2 _sph_particle_viscosity_force(
        int i, 
        const std::vector<Particle> &p, 
        const std::vector<int> &neighborhood, 
        const std::vector<double> &density
    );

    double _sph_w_density(double dist);
    double _sph_w_pressure_grad(double dist);
    double _sph_w_viscosity_laplacian(double dist);
};