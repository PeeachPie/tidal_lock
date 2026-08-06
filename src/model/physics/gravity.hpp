#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include "utils.hpp"

#include "constants.hpp"
#include "barnes_hut.hpp"
#include "particle.hpp"

struct GravitySettings {
    double p_radius;
    double theta;
    double bound;
};

class GravitySolver {
private:
    double p_radius_;
    double theta_;
    double bound_;
public:
    GravitySolver(const GravitySettings& settings);

    std::pair<std::vector<glm::dvec2>, std::map<std::string, glm::dvec2>> calc_forces(
        const std::vector<Particle> &p, 
        const std::map<std::string, Planet> &planet
    );
};