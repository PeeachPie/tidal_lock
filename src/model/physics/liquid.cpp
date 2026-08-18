#include "liquid.hpp"

LiquidsSolver::LiquidsSolver(LiquidsSettings &settings): 
    p_smoothing_length_(settings.p_smoothing_length), 
    p_density_in_rest_(settings.p_density_in_rest),
    p_temperature_k_(settings.p_temperature_k),
    p_viscousity_k_(settings.p_viscousity_k) {};

double LiquidsSolver::_sph_w_density(double dist) {
    if (dist >= p_smoothing_length_) return 0;
    return 4 * pow((1 - (dist * dist) / (p_smoothing_length_ * p_smoothing_length_)), 3) / (PI * p_smoothing_length_ * p_smoothing_length_);
}

double LiquidsSolver::_sph_w_pressure_grad(double dist) {
    if (dist >= p_smoothing_length_) return 0;
    return -30 * pow(1 * (1 - dist / p_smoothing_length_), 2) / (PI * pow(p_smoothing_length_, 3));
}

double LiquidsSolver::_sph_w_viscosity_laplacian(double dist) {
    if (dist >= p_smoothing_length_) return 0;
    return 40 * (1 - dist / p_smoothing_length_) / (PI * pow(p_smoothing_length_, 4));
}

double LiquidsSolver::_sph_calc_particle_density(int i, const std::vector<Particle> &p, const std::vector<int> &neighborhood) {
    double density = 0.0;

    for (int j: neighborhood) {
        density += p[j].m * _sph_w_density(glm::distance(p[i].pos, p[j].pos));
    }

    return density;
}

std::vector<double> LiquidsSolver::_sph_calc_density(const std::vector<Particle> &p, const std::vector<std::vector<int>> &neighborhoods) {
    std::vector<double> density(p.size());

    #pragma omp parallel for
    for (int i = 0; i < p.size(); i++) {
        density[i] = _sph_calc_particle_density(i, p, neighborhoods[i]);
    }

    return density;
}

std::vector<double> LiquidsSolver::_sph_calc_pressure(const std::vector<Particle> &p, const std::vector<double> &density, const std::vector<std::vector<int>> &neighborhoods) {
    std::vector<double> pressure(p.size());

    if (p_temperature_k_ != 0.0) {
        #pragma omp parallel for
        for (int i = 0; i < p.size(); i++) {
            pressure[i] = p_temperature_k_ * glm::max(0.0, (density[i] - p_density_in_rest_));
        }
    }

    return pressure;
}

glm::dvec2 LiquidsSolver::_sph_particle_pressure_force(int i, const std::vector<Particle> &p, const std::vector<int> &neighborhood, const std::vector<double> &density, const std::vector<double> &pressure) {
    glm::dvec2 pressure_force(0);

    for (int j: neighborhood) {
        if (i == j) continue;

        glm::dvec2 r = p[i].pos - p[j].pos;
        glm::dvec2 n = glm::normalize(r);

        pressure_force += p[j].m * (pressure[i] + pressure[j]) / (2 * density[j]) * _sph_w_pressure_grad(glm::length(r)) * n;
    }

    return -pressure_force;
}

// не сохраняет угловой момент
// glm::dvec2 LiquidsSolver::_sph_particle_viscosity_force(int i, const std::vector<Particle> &p, const std::vector<int> &neighborhood, const std::vector<double> &density) {
//     glm::dvec2 viscosity_force(0);

//     for (int j: neighborhood) {
//         if (p[j].v - p[i].v == glm::dvec2{ 0, 0 }) continue;
//         glm::dvec2 r = p[i].pos - p[j].pos;
//         viscosity_force += p[j].m * (p[j].v - p[i].v) / density[j] * _sph_w_viscosity_laplacian(glm::length(r));
//     }

//     return p_viscousity_k_ * viscosity_force;
// }

glm::dvec2 LiquidsSolver::_sph_particle_viscosity_force(
    int i, 
    const std::vector<Particle> &p,
    const std::vector<int> &neighborhood,
    const std::vector<double> &density
) {
    glm::dvec2 viscosity_force(0.0);

    for (int j : neighborhood) {
        if (i == j) continue;
        glm::dvec2 r = p[i].pos - p[j].pos;
        double dist = glm::length(r);
        glm::dvec2 n = r / dist;

        glm::dvec2 u_ij = p[j].v - p[i].v;
        double proj = glm::dot(u_ij, n) / dist;

        viscosity_force += p[j].m * (DIMS + 2.0) * proj * n * _sph_w_pressure_grad(dist) / density[j];
    }

    return p_viscousity_k_ * viscosity_force;
}


std::pair<std::vector<glm::dvec2>, std::vector<glm::dvec2>> LiquidsSolver::calc_forces(const std::vector<Particle> &p) {
    std::vector<glm::dvec2> p_forces(p.size(), { 0.0, 0.0 });
    std::vector<glm::dvec2> v_forces(p.size(), { 0.0, 0.0 });

    if (p_viscousity_k_ != 0.0 || p_temperature_k_ != 0.0) {
        std::vector<std::vector<int>> neighborhoods = get_neighborhoods(p, p_smoothing_length_);

        std::vector<double> density = _sph_calc_density(p, neighborhoods);
        std::vector<double> pressure = _sph_calc_pressure(p, density, neighborhoods);

        #pragma omp parallel for
        for (int i = 0; i < p.size(); i++) {
            glm::dvec2 pressure_force(0, 0);
            if (p_temperature_k_ != 0.0)
                p_forces[i] = _sph_particle_pressure_force(i, p, neighborhoods[i], density, pressure) * p[i].m / density[i];

            glm::dvec2 viscosity_force(0, 0);
            if (p_viscousity_k_ != 0.0)
                v_forces[i] = _sph_particle_viscosity_force(i, p, neighborhoods[i], density) * p[i].m / density[i];
        }
    }

    return {
        p_forces,
        v_forces
    };
}