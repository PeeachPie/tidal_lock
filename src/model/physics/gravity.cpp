#include "gravity.hpp"

GravitySolver::GravitySolver(const GravitySettings &settings): 
    p_radius_(settings.p_radius), 
    theta_(settings.theta), 
    bound_(settings.bound)
{};

std::pair<std::vector<glm::dvec2>, std::map<std::string, glm::dvec2>> 
GravitySolver::calc_forces(
    const std::vector<Particle> &p, 
    const std::map<std::string, Planet> &planet
) {
    BarnesHutQuadTree tree(p_radius_, theta_, bound_);
    for (int i = 0; i < p.size(); i++) {
        tree.insert(p[i], i);
    }

    int id = p.size();
    for (auto &pr: planet) {
        const Planet& pl = pr.second;
        tree.insert({ pl.pos, pl.m }, id);
        id++;
    }

    std::vector<glm::dvec2> p_forces(p.size());
    std::map<std::string, glm::dvec2> pl_forces;

    #pragma omp parallel for
    for (int i = 0; i < p.size(); i++) {
        p_forces[i] = tree.calc_force(p[i], i);
    }

    id = p.size();
    for (auto &pr: planet) {
        const Planet& pl = pr.second;
        pl_forces[pr.first] = tree.calc_force({ pl.pos, pl.m }, id);
        id++;
    }

    return {
        p_forces,
        pl_forces
    };
}