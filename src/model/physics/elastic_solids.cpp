#include "elastic_solids.hpp"

ElasticSolidsSolver::ElasticSolidsSolver(ElasticSolidsSettings &settings): 
    p_smoothing_length_(settings.p_smoothing_length), 
    p_shear_modulus_(settings.p_shear_modulus),
    p_bulk_modulus_(settings.p_bulk_modulus) {};

std::vector<std::unordered_map<int, glm::dvec2>> ElasticSolidsSolver::_calc_initial_offsets(
    const std::vector<Particle> &p,
    const std::vector<std::vector<int>> &neighborhoods
) {
    std::vector<std::unordered_map<int, glm::dvec2>> offesets(p.size());

    for (int i = 0; i < p.size(); i++) {
        for (int j: neighborhoods[i]) {
            offesets[i][j] = p[i].pos - p[j].pos;
        }
    }

    return offesets;
}

double ElasticSolidsSolver::_sph_w_density(double dist) {
    if (dist >= p_smoothing_length_) return 0;
    return 4 * pow((1 - (dist * dist) / (p_smoothing_length_ * p_smoothing_length_)), 3) / PI;
}

double ElasticSolidsSolver::_sph_calc_particle_density(
    int i, 
    const std::vector<Particle> &p, 
    const std::vector<int> &neighborhood
) {
    double density = 0.0;

    for (int j: neighborhood) {
        density += p[j].m * _sph_w_density(glm::distance(p[i].pos, p[j].pos));
    }

    return density;
}

std::vector<double> ElasticSolidsSolver::_sph_calc_density(
    const std::vector<Particle> &p, 
    const std::vector<std::vector<int>> &neighborhoods
) {
    std::vector<double> density(p.size());

    #pragma omp parallel for
    for (int i = 0; i < p.size(); i++) {
        density[i] = _sph_calc_particle_density(i, p, neighborhoods[i]);
    }

    return density;
}

void ElasticSolidsSolver::_fix_initial_state(const std::vector<Particle> &p) {
    initial_neighborhoods_ = get_neighborhoods(p, p_smoothing_length_);
    initial_offsets_ = _calc_initial_offsets(p, initial_neighborhoods_);
    rest_volumes_ = _sph_calc_volumes(p, _sph_calc_density(p, initial_neighborhoods_));
    precalc_grads_ = _sph_calc_corrected_kernel_gradients(p);
    initial_state_fixed_ = true;
}

double ElasticSolidsSolver::_sph_w_grad(double dist) {
    if (dist >= p_smoothing_length_) return 0;
    return -30 * pow(1 * (1 - dist / p_smoothing_length_), 2) / PI;
}

std::vector<double> ElasticSolidsSolver::_sph_calc_volumes(const std::vector<Particle> &p, const std::vector<double> &density) {
    std::vector<double> volume(density.size());

    #pragma omp parallel for
    for (int i = 0; i < volume.size(); i++) {
        volume[i] = p[i].m / density[i];
    }

    return volume;
}

glm::dmat2 ElasticSolidsSolver::_sph_correction_matrix(int i, const std::vector<Particle> &p) {
    glm::dmat2 correction_matrix(0.0, 0.0, 0.0, 0.0);

    for (int j: initial_neighborhoods_[i]) {
        if (i == j) continue;

        glm::dvec2 r = p[i].pos - p[j].pos;
        glm::dvec2 n = glm::normalize(r);

        correction_matrix += tensor_product(rest_volumes_[j] * (_sph_w_grad(glm::length(r)) * n), -r);
    }

    return moore_penrose_inverse(correction_matrix);
}

std::vector<std::unordered_map<int, glm::dvec2>> 
ElasticSolidsSolver::_sph_calc_corrected_kernel_gradients(const std::vector<Particle> &p) {
    std::vector<std::unordered_map<int, glm::dvec2>> corrected_kernel_gradients(p.size());

    // TODO: #pragma omp parallel for
    for (int i = 0; i < p.size(); i++) {
        glm::dmat2 L = _sph_correction_matrix(i, p);
        // TODO: #pragma omp parallel for
        for (int j: initial_neighborhoods_[i]) {
            if (i == j) corrected_kernel_gradients[i][j] = glm::dvec2(0.0, 0.0);
            else {
                // TODO: возможно другой градиент
                // TODO: n вычислять в функции градиента
                glm::dvec2 r = p[i].pos - p[j].pos;
                glm::dvec2 n = glm::normalize(r);

                corrected_kernel_gradients[i][j] = L * (_sph_w_grad(glm::length(r)) * n);
            }
        }
    }

    return corrected_kernel_gradients;
}

glm::dmat2 ElasticSolidsSolver::_sph_calc_particle_deformation_grad(int i, const std::vector<Particle> &p) {
    glm::dmat2 deformation_grad(0.0, 0.0, 0.0, 0.0);

    for (int j: initial_neighborhoods_[i]) {
        // TODO: i==j?
        glm::dvec2 r = p[j].pos - p[i].pos;
        deformation_grad += tensor_product(rest_volumes_[j] * r, precalc_grads_[i].at(j));
    }

    return deformation_grad;
}

std::pair<glm::dmat2, glm::dmat2> 
ElasticSolidsSolver::_sph_calc_particle_corotated_deformation_grad(int i, const std::vector<Particle> &p) {
    glm::dmat2 F = _sph_calc_particle_deformation_grad(i, p);
    glm::dmat2 R = extract_rotation(F);

    glm::dmat2 F_corotated(1.0);

    for (int j: initial_neighborhoods_[i]) {
        F_corotated += tensor_product(rest_volumes_[j] * ((p[j].pos - p[i].pos) - R * initial_offsets_[j].at(i)), R * precalc_grads_[i].at(j));
    }

    return { F_corotated, R };
}

glm::dmat2 ElasticSolidsSolver::_infinitesimal_strain_tensor(glm::dmat2 deformation_grad) {
    return 0.5 * (deformation_grad + glm::transpose(deformation_grad)) - glm::dmat2(1.0);
}

glm::dmat2 ElasticSolidsSolver::_piola_kirchhoff_stress_tensor(glm::dmat2 strain_tensor) {
    return 
        2.0 * p_shear_modulus_ * strain_tensor + 
        (p_bulk_modulus_ - (2.0/3.0) * p_shear_modulus_) * trace(strain_tensor) * glm::dmat2(1.0);
}

std::pair<std::vector<glm::dmat2>, std::vector<glm::dmat2>> 
ElasticSolidsSolver::_sph_calc_stress_tensors(const std::vector<Particle> &p) {
    std::vector<glm::dmat2> P(p.size());
    std::vector<glm::dmat2> R(p.size());

    #pragma omp parallel for
    for (int i = 0; i < P.size(); i++) {
        auto [F_corotated, Ri] = _sph_calc_particle_corotated_deformation_grad(i, p);

        glm::dmat2 e = _infinitesimal_strain_tensor(F_corotated);

        P[i] = _piola_kirchhoff_stress_tensor(e);
        R[i] = Ri;
    }

    return { P, R };
}

glm::dvec2 ElasticSolidsSolver::_sph_particle_elastic_force(
    int i,
    const std::vector<glm::dmat2> &stress_tensors, 
    const std::vector<glm::dmat2> &rotation
) {
    glm::dvec2 elastic_force(0.0, 0.0);

    for (int j: initial_neighborhoods_[i]) {
        elastic_force += rest_volumes_[i] * rest_volumes_[j] * (
            stress_tensors[i] * rotation[i] * precalc_grads_[i].at(j) - 
            stress_tensors[j] * rotation[j] * precalc_grads_[j].at(i)
        );
    }

    return elastic_force;
}

std::vector<glm::dvec2> ElasticSolidsSolver::calc_forces(const std::vector<Particle> &p) {
    std::vector<glm::dvec2> forces(p.size(), { 0.0, 0.0 });

    if (initial_state_fixed_ == false) {
        // std::vector<Particle> false_init;

        // for (int i = 0; i < p.size(); i++) {
        //     Particle pt(
        //         // glm::dvec2{
        //         //     radius * std::cos(angle), 
        //         //     radius * std::sin(angle)
        //         // },
        //         0.8 * p[i].pos,
        //         { 0, 0 },
        //         { 0, 0 }
        //     );

        //     false_init.push_back(pt);
        // }

        // _fix_initial_state(false_init);

        _fix_initial_state(p);
    }

    if (p_bulk_modulus_ != 0.0 || p_shear_modulus_ != 0.0) {
        auto [P, R] = _sph_calc_stress_tensors(p);

        #pragma omp parallel for
        for (int i = 0; i < p.size(); i++) {
            glm::dvec2 elastic_force = _sph_particle_elastic_force(i, P, R);

            // std::cout << "elastic_force=" << elastic_force.x << ' ' << elastic_force.y << '\n';

            forces[i] = elastic_force;
        }
    }

    return forces;
}