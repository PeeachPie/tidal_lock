#include "model.hpp"

Model::Model(double h, double time_delta, double scale): 
    p_smoothing_lenght_(h), time_step_(time_delta), scale_(scale) {};

void Model::add_planet(std::string name, Planet planet) {
    planets_[name] = planet;
}

void Model::add_water_to_planet(std::string planet_name, int p_count, double h) {
    const Planet planet = planets_[planet_name];

    for (int i = 0; i < p_count; i++) {
        double angle = gen_angle();
        double height = gen_h(h);

        Particle p(glm::dvec2{
            (planet.r + height) * std::cos(angle) + planet.pos.x, 
            (planet.r + height) * std::sin(angle) + planet.pos.y
        });

        add_particle(p);
    }

    // p_density_in_rest = _calc_density_in_rest();
}

void Model::add_particle(Particle particle) {
    p_.push_back(particle);
}

void Model::_planets_next_step() {
    for (auto &p: planets_) {
        Planet& planet = p.second;

        if (planet.v) {
            double ang = current_time * (planet.v / planet.rotation_r);

            planet.pos = planet.rotation_focus + glm::dvec2{
                planet.rotation_r * glm::cos(ang), 
                planet.rotation_r * glm::sin(ang)
            };
        }
    }
}

std::unordered_map<int, std::unordered_map<int, std::vector<int>>> Model::_get_particle_grid() {
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> particle_grid;

    for (int i = 0; i < p_.size(); i++) {
        particle_grid[
            floor(p_[i].pos.x / p_smoothing_lenght_)
        ][
            floor(p_[i].pos.y / p_smoothing_lenght_)
        ].push_back(i);
    }

    return particle_grid;
}

std::vector<int> Model::_get_particle_neighborhood(int pi, const std::unordered_map<int, std::unordered_map<int, std::vector<int>>> &particle_grid) {
    std::vector<int> neighborhood;

    int x = floor(p_[pi].pos.x / p_smoothing_lenght_);
    int y = floor(p_[pi].pos.y / p_smoothing_lenght_);

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (particle_grid.count(x + i) && particle_grid.at(x + i).count(y + j)) {
                for (int k: particle_grid.at(x + i).at(y + j)) {
                    glm::dvec2 r = p_[pi].pos - p_[k].pos;
                    if (glm::dot(r, r) <= p_smoothing_lenght_ * p_smoothing_lenght_) {
                        neighborhood.push_back(k);
                    }
                }
            }
        }
    }

    return neighborhood;
}

std::vector<std::vector<int>> Model::_get_neighborhoods() {
    std::vector<std::vector<int>> neighborhoods(p_.size());
    
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> particle_grid = _get_particle_grid();
    
    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        neighborhoods[i] = _get_particle_neighborhood(i, particle_grid);
    }

    return neighborhoods;
}

glm::dvec2 Model::_particle_gravity_forces(int i) {
    glm::dvec2 gravity_forces(0);

    for (auto &p: planets_) {
        const Planet& planet = p.second;

        double r = glm::distance(p_[i].pos, planet.pos);
        glm::dvec2 n = (planet.pos - p_[i].pos) / r;
        glm::dvec2 f = (G * (planet.m * p_mass_) / (r * r)) * n;

        gravity_forces += f;
    }

    return gravity_forces;
}

glm::dvec2 Model::_get_barycenter() {
    glm::dvec2 mr(0);
    double tm = 0;

    for (auto &p: planets_) {
        const Planet& planet = p.second;
        mr += planet.m * planet.pos;
        tm += planet.m;
    }

    glm::dvec2 barycenter = mr / tm;

    return barycenter;
}

glm::dvec2 Model::_particle_centrifugal_forces(int i, glm::dvec2 barycenter) {
    // glm::dvec2 barycenter = planets_["Земля"].pos + (planets_["Луна"].pos - planets_["Земля"].pos) / 10.0;
    // double angular_velocity = MOON_VELOCITY / (glm::distance(barycenter, planets_["Луна"].pos));

    glm::dvec2 r = p_[i].pos - barycenter;
    double l = glm::length(r);
    glm::dvec2 n = r / l;

    glm::dvec2 centrifugal_force = l * (planets_["Луна"].a_v * planets_["Луна"].a_v) * p_mass_ * n;

    return centrifugal_force;
}

double Model::_sph_w_density(double dist) {
    if (dist >= p_smoothing_lenght_) return 0;
    return 4 * pow((1 - (dist * dist) / (p_smoothing_lenght_ * p_smoothing_lenght_)), 3) / PI;
}

double Model::_sph_w_pressure_grad(double dist) {
    if (dist >= p_smoothing_lenght_) return 0;
    return -30 * pow(1 * (1 - dist / p_smoothing_lenght_), 2) / PI;
}

double Model::_sph_w_viscosity_laplacian(double dist) {
    if (dist >= p_smoothing_lenght_) return 0;
    return 40 * (1 - dist / p_smoothing_lenght_) / PI;
}

// double Model::_sph_w_density(double dist) {
//     if (dist >= p_smoothing_lenght_) return 0;
//     return 4 * pow((1 - (dist * dist) / (p_smoothing_lenght_ * p_smoothing_lenght_)), 3) / (PI * pow(1, 8));
// }

// double Model::_sph_w_pressure_grad(double dist) {
//     if (dist >= p_smoothing_lenght_) return 0;
//     return -30 * pow(1 - dist / p_smoothing_lenght_, 2) / (PI * pow(1, 5));
// }

// double Model::_sph_w_viscosity_laplacian(double dist) {
//     if (dist >= p_smoothing_lenght_) return 0;
//     return 40 * (1 - dist / p_smoothing_lenght_) / (PI * pow(1, 5));
// }

double Model::_calc_density_in_rest() {
    double density = 0.0;

    auto pg = _get_particle_grid();
    for (int i = 0; i < p_.size(); i++) {
        density += _sph_calc_particle_density(i, _get_particle_neighborhood(i, pg));
    }

    return density / p_.size();
}

double Model::_sph_calc_particle_density(int i, const std::vector<int> &neighborhood) {
    double density = 0.0;

    for (int j: neighborhood) {
        density += p_mass_ * _sph_w_density(glm::distance(p_[i].pos, p_[j].pos));
    }

    return density;
}

std::vector<double> Model::_sph_calc_density(const std::vector<std::vector<int>> &neighborhoods) {
    std::vector<double> density(p_.size());

    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        density[i] = _sph_calc_particle_density(i, neighborhoods[i]);
    }

    return density;
}

std::vector<double> Model::_sph_calc_pressure(std::vector<double> density, const std::vector<std::vector<int>> &neighborhoods) {
    std::vector<double> pressure(p_.size());

    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        pressure[i] = p_temperature_k * glm::max(0.0, (density[i] - p_density_in_rest));
    }

    return pressure;
}

glm::dvec2 Model::_particle_pressure_forces(int i, const std::vector<int> &neighborhood, const std::vector<double> &density, const std::vector<double> pressure) {
    glm::dvec2 pressure_force(0);

    for (int j: neighborhood) {
        if (i == j) continue;

        glm::dvec2 r = p_[i].pos - p_[j].pos;
        glm::dvec2 n = glm::normalize(r);

        pressure_force += p_mass_ * (pressure[i] + pressure[j]) / (2 * density[j]) * _sph_w_pressure_grad(glm::length(r)) * n;
    }

    return -pressure_force;
}

glm::dvec2 Model::_particle_viscosity_forces(int i, const std::vector<int> &neighborhood, const std::vector<double> &density) {
    glm::dvec2 viscosity_force(0);

    for (int j: neighborhood) {
        if (p_[j].v - p_[i].v == glm::dvec2{ 0, 0 }) continue;
        glm::dvec2 r = p_[i].pos - p_[j].pos;
        viscosity_force += p_mass_ * (p_[j].v - p_[i].v) / density[j] * _sph_w_viscosity_laplacian(glm::length(r));
    }

    return VISCOSITY_K * viscosity_force;
}

void Model::next_step() {
    current_time += time_step_;

    _planets_next_step();

    std::vector<std::vector<int>> neighborhoods = _get_neighborhoods();

    std::vector<double> density = _sph_calc_density(neighborhoods);
    std::vector<double> pressure = _sph_calc_pressure(density, neighborhoods);
    glm::dvec2 barycenter = _get_barycenter();

    std::vector<Particle> new_p_(p_.size());

    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        glm::dvec2 gravity_forces = _particle_gravity_forces(i);
        glm::dvec2 centrifugal_forces = _particle_centrifugal_forces(i, barycenter);
        glm::dvec2 pressure_forces = _particle_pressure_forces(i, neighborhoods[i], density, pressure);
        glm::dvec2 viscosity_forces = _particle_viscosity_forces(i, neighborhoods[i], density);
        
        // std::cout << '\n';
        // std::cout << "density=" << density[i] << '\n';
        // std::cout << "mass=" << p_mass_ << '\n';
        // std::cout << "pressure=" << pressure[i] << '\n';
        // std::cout << "pressure_f=" << pressure_forces.x << ' ' << pressure_forces.y << '\n';
        // std::cout << "gravity_f=" << gravity_forces.x << ' ' << gravity_forces.y << '\n';
        // std::cout << "centrifugal_f=" << centrifugal_forces.x << ' ' << centrifugal_forces.y << '\n';
        // std::cout << "viscosity_f=" << viscosity_forces.x << ' ' << viscosity_forces.y << '\n';

        glm::dvec2 a = 
        (
            gravity_forces + 
            centrifugal_forces +
            0.0
        ) 
        / p_mass_ +
        (
            pressure_forces + 
            viscosity_forces
        ) 
        / density[i];

        // glm::dvec2 a(0);

        glm::dvec2 new_pos = p_[i].pos + p_[i].v * time_step_ + 0.5 * p_[i].a * (time_step_ * time_step_);
        glm::dvec2 new_v = p_[i].v + 0.5 * (p_[i].a + a) * time_step_;

        for (auto &p: planets_) {
            const Planet& planet = p.second;

            double dist = point_to_segment_dist(
                p_[i].pos, 
                new_pos, 
                planet.pos
            );

            if (planet.r >= dist) {
                double r = glm::distance(p_[i].pos, planet.pos);
                glm::dvec2 n = glm::normalize(planet.pos - p_[i].pos);

                new_v += 2 * (-glm::dot(n, new_v)) * n;
                new_pos = p_[i].pos + new_v * time_step_ + 0.5 * p_[i].a * (time_step_ * time_step_);
                a = glm::dvec2(0);
            }
        }

        new_p_[i].v = new_v;
        new_p_[i].pos = new_pos;
        new_p_[i].a = a;
    }

    p_ = new_p_;
}

double Model::calc_energy() const {
    double kinetic = 0.0;
    double potential = 0.0;

    for (const auto& particle : p_) {
        kinetic += 0.5 * p_mass_ * glm::dot(particle.v, particle.v);

        for (const auto& [name, planet] : planets_) {
            double r = glm::distance(particle.pos, planet.pos);

            potential -= G * p_mass_ * planet.m / r;
        }
    }

    return kinetic + potential;
}

Frame Model::get_frame() {
    Frame frame;

    frame.particle_r = p_r_;

    for (auto &p: planets_) {
        frame.planets.push_back(
            { 
                p.second.pos / scale_, 
                p.second.r / scale_ 
            }
        );
    }

    for (auto &p: p_) {
        frame.particles.push_back(
            { p.pos / scale_ }
        );
    }

    return frame;
}