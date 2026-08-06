#include "model.hpp"

Model::Model(ModelSettings &settings):
    liquids_solver{settings.liquids_settings},
    elastic_solids_solver{settings.elastic_solids_settings},
    gravity_solver{settings.gravity_settings},
    time_step_(settings.time_step), 
    scale_(settings.scale) {};

void Model::add_planet(std::string name, Planet planet) {
    assert(planet.m != 0.0);
    planets_[name] = planet;
}

void Model::_add_water_to_planet(const Planet &planet, int p_count, double p_mass, double h) {
    for (int i = 0; i < p_count; i++) {
        double angle = gen_angle();
        double radius = gen_radius(planet.r, h);

        glm::dvec2 rotation_speed = {
            -planet.av * radius * std::sin(angle),
            planet.av * radius * std::cos(angle)
        };

        Particle p(
            glm::dvec2{
                radius * std::cos(angle) + planet.pos.x, 
                radius * std::sin(angle) + planet.pos.y
            },
            p_mass,
            planet.v + rotation_speed,
            planet.a
        );

        _add_particle(p);
    }
}

void Model::add_water_to_planet(std::string planet_name, int p_count, double p_mass, double h) {
    const Planet planet = planets_.at(planet_name);
    _add_water_to_planet(planet, p_count, p_mass, h);
}

void Model::add_particle_planet(Planet planet, int p_count) {
    double r = planet.r;
    planet.r = 0;

    _add_water_to_planet(planet, p_count, planet.m/p_count, r);
}

void Model::add_particle_satellite(std::string planet_name, Planet satellite, int p_count, double velocity) {
    const Planet& planet = planets_[planet_name];
    double p_mass = satellite.m / p_count;

    double R = glm::distance(planet.pos, satellite.pos);

    for (int i = 0; i < p_count; i++) {
        double angle = gen_angle();
        double radius = gen_radius(0, satellite.r);

        glm::dvec2 rotation_speed = {
            -satellite.av * radius * std::sin(angle),
            satellite.av * radius * std::cos(angle)
        };

        glm::dvec2 pos {
            radius * std::cos(angle) + satellite.pos.x, 
            radius * std::sin(angle) + satellite.pos.y
        };

        glm::dvec2 r = pos - planet.pos;
        glm::dvec2 n = glm::normalize(r);

        Particle p(
            pos,
            p_mass,
            velocity * (glm::length(r) / R) * glm::dvec2{ -n.y, n.x } + rotation_speed,
            satellite.a
        );

        _add_particle(p);
    }
}

void Model::_add_particle(Particle particle) {
    assert(particle.m != 0.0);
    p_.push_back(particle);
}

void Model::_planets_next_step(const std::map<std::string, glm::dvec2> &forces) {
    std::map<std::string, Planet> new_planets_;

    // #pragma omp parallel for
    for (auto &pr: planets_) {
        if (pr.second.m == 0) continue;

        glm::dvec2 a = forces.at(pr.first) / pr.second.m;

        glm::dvec2 new_pos = pr.second.pos +  pr.second.v * time_step_ + 0.5 *  pr.second.a * (time_step_ * time_step_);
        glm::dvec2 new_v =  pr.second.v + 0.5 * (pr.second.a + a) * time_step_;

        // TODO: скорее всего это неправда
        for (int i = 0; i < p_.size(); i++) {
            double dist = glm::distance(p_[i].pos, new_pos);

            if (pr.second.r >= dist) {
                // double r = glm::distance(p_[i].pos, pr.second.pos);
                glm::dvec2 n = glm::normalize(pr.second.pos - p_[i].pos);

                p_[i].v = glm::dot(n, new_v) * -n;
                p_[i].a = a;
                p_[i].pos = pr.second.pos + pr.second.r * -n;
                p_[i].pos = p_[i].pos + p_[i].v * time_step_ + 0.5 * p_[i].a * (time_step_ * time_step_);
            }
        }

        new_planets_[pr.first] = { 
            new_pos, 
            pr.second.r, 
            pr.second.m, 
            new_v, 
            a,
            pr.second.av
        };
    }

    planets_ = new_planets_;
}

void Model::_particles_next_step(const std::vector<glm::dvec2> &forces) {
    std::vector<Particle> new_p_(p_.size());

    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        glm::dvec2 a = forces[i] / p_[i].m;

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

        new_p_[i] = { new_pos, p_[i].m, new_v, a };
    }

    p_ = new_p_;
}

void Model::next_step() {
    std::vector<glm::dvec2> l_forces = liquids_solver.calc_forces(p_);
    std::vector<glm::dvec2> e_forces = elastic_solids_solver.calc_forces(p_);
    auto [g_forces, pl_forces] = gravity_solver.calc_forces(p_, planets_);

    std::vector<glm::dvec2> p_forces(p_.size(), {0.0, 0.0});
    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        p_forces[i] = e_forces[i] + l_forces[i] + g_forces[i];
    }

    _planets_next_step(pl_forces);
    _particles_next_step(p_forces);
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