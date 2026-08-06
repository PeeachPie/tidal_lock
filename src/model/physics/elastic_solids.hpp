#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include "utils.hpp"

#include "constants.hpp"
#include "neighborhood_search.hpp"
#include "particle.hpp"

struct ElasticSolidsSettings {
    double p_smoothing_length;
    double p_shear_modulus;
    double p_bulk_modulus;
};

class ElasticSolidsSolver {
private:
    double p_smoothing_length_;
    double p_shear_modulus_;
    double p_bulk_modulus_;

    bool initial_state_fixed_ = false;

    std::vector<double> rest_volumes_;
    std::vector<std::vector<int>> initial_neighborhoods_;
    std::vector<std::unordered_map<int, glm::dvec2>> precalc_grads_;
    std::vector<std::unordered_map<int, glm::dvec2>> initial_offsets_;

    double p_r_ = 0.05 * 40;

    double scale_;

public:
    ElasticSolidsSolver(ElasticSolidsSettings &settings);
    std::vector<glm::dvec2> calc_forces(const std::vector<Particle> &p);

private:
    void _fix_initial_state(const std::vector<Particle> &p);

    double _sph_w_density(double dist);
    double _sph_w_grad(double dist);

    double _sph_calc_particle_density(
        int i, 
        const std::vector<Particle> &p, 
        const std::vector<int> &neighborhood
    );
    std::vector<double> _sph_calc_density(
        const std::vector<Particle> &p, 
        const std::vector<std::vector<int>> 
        &neighborhoods
    );

    std::vector<double> _sph_calc_volumes(const std::vector<Particle> &p, const std::vector<double> &density);

    std::vector<std::unordered_map<int, glm::dvec2>> _calc_initial_offsets(
        const std::vector<Particle> &p,
        const std::vector<std::vector<int>> &neighborhoods
    );
    
    glm::dmat2 _sph_correction_matrix(int i, const std::vector<Particle> &p);

    std::vector<std::unordered_map<int, glm::dvec2>> _sph_calc_corrected_kernel_gradients(
        const std::vector<Particle> &p
    );

    glm::dmat2 _sph_calc_particle_deformation_grad(
        int i,
        const std::vector<Particle> &p
    );

    std::pair<glm::dmat2, glm::dmat2> _sph_calc_particle_corotated_deformation_grad(
        int i, 
        const std::vector<Particle> &p
    );

    glm::dmat2 _infinitesimal_strain_tensor(glm::dmat2 deformation_grad);
    glm::dmat2 _piola_kirchhoff_stress_tensor(glm::dmat2 infinitesimal_strain_tensor);

    std::pair<std::vector<glm::dmat2>, std::vector<glm::dmat2>> _sph_calc_stress_tensors(
        const std::vector<Particle> &p
    );

    glm::dvec2 _sph_particle_elastic_force(
        int i,
        const std::vector<glm::dmat2> &stress_tensors, 
        const std::vector<glm::dmat2> &rotation
    );
};