#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <glm/glm.hpp>
#include "utils.hpp"
#include "constants.hpp"
#include <iostream>

struct Particle {
    glm::dvec2 pos;
    glm::dvec2 v;
    glm::dvec2 a;

    Particle(): pos(0), v(0), a(0) {};
    Particle(glm::dvec2 pos, glm::dvec2 v=glm::dvec2(0, 0), glm::dvec2 a=glm::dvec2(0, 0)): 
        pos(pos), v(v), a(a) {};
};

struct ParticleFrame {
    glm::dvec2 pos;
};

struct PlanetFrame {
    glm::dvec2 pos;
    double r;
};

struct Planet {
    glm::dvec2 pos;
    glm::dvec2 rotation_focus;
    double rotation_r;
    double r;
    double m;
    double a_v;
    double v;
};

struct Frame {
    std::vector<PlanetFrame> planets;
    std::vector<ParticleFrame> particles;

    float particle_r;
};

class Model {
private:
    // double h_;
    double time_step_;
    double current_time = 0;

    std::map<std::string, Planet> planets_;
    std::vector<Particle> p_;

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
    glm::dvec2 _get_barycenter();
    glm::dvec2 _particle_gravity_forces(int i);
    glm::dvec2 _particle_centrifugal_forces(int i, glm::dvec2 barycenter);
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