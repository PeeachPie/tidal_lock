#include "model.hpp"

Model::Model(ModelSettings &settings):
    liquids_solver{settings.liquids_settings},
    elastic_solids_solver{settings.elastic_solids_settings},
    gravity_solver{settings.gravity_settings},
    time_step_(settings.time_step), 
    scale_(settings.scale),
    planet_scale_(settings.planet_scale),
    debug_(settings.debug) 
{};

void Model::add_planet(std::string name, Planet planet) {
    assert(planet.m != 0.0);
    planets_[name] = planet;
}

glm::dvec2 Model::_get_barycenter() {
    glm::dvec2 s{0, 0};
    double tm = 0.0;

    for (int i = 0; i < p_.size(); i++) {
        tm += p_[i].m;
        s += p_[i].m * p_[i].pos;
    }

    return s / tm;
}

std::pair<double, double> Model::_get_angular_velocity(glm::dvec2 satellite_bc, std::string planet_name) {
    glm::dvec2 v = { 0, 0 };

    for (int i = 0; i < p_.size(); i++) {
        v += p_[i].v;
    }

    v /= p_.size();

    double av = 0.0;

    for (int i = 0; i < p_.size(); i++) {
        glm::dvec2 r = p_[i].pos - satellite_bc;
        glm::dvec2 n = glm::normalize(r);

        av += glm::dot({-n.y, n.x }, p_[i].v - v) / glm::length(r);
    }

    glm::dvec2 R = satellite_bc - planets_[planet_name].pos;
    glm::dvec2 N = glm::normalize(R);

    return {
        av / p_.size(),
        glm::dot({-N.y, N.x }, v) / glm::length(R)
    };
}

// равномерное случайное распределение
// void Model::_add_water_to_planet(const Planet &planet, int p_count, double p_mass, double h) {
//     for (int i = 0; i < p_count; i++) {
//         double angle = gen_angle();
//         double radius = gen_radius(planet.r, h);

//         glm::dvec2 rotation_speed = {
//             -planet.av * radius * std::sin(angle),
//             planet.av * radius * std::cos(angle)
//         };

//         Particle p(
//             glm::dvec2{
//                 radius * std::cos(angle) + planet.pos.x, 
//                 radius * std::sin(angle) + planet.pos.y
//             },
//             p_mass,
//             planet.v + rotation_speed,
//             planet.a
//         );

//         _add_particle(p);
//     }
// }

void Model::_add_water_to_planet(const Planet &planet, double m, double h, double d) {
    std::vector<glm::dvec2> points = gen_hex_ring_points(planet.r, planet.r + h, d);

    std::cout << points.size() << std::endl;

    if (points.empty()) {
        return;
    }

    double p_mass = m / points.size();

    for (const glm::dvec2 &local_pos : points) {
        double radius = std::sqrt(local_pos.x * local_pos.x + local_pos.y * local_pos.y);
        double angle = std::atan2(local_pos.y, local_pos.x);

        glm::dvec2 rotation_speed = {
            -planet.av * radius * std::sin(angle),
             planet.av * radius * std::cos(angle)
        };

        Particle p(
            glm::dvec2{local_pos.x + planet.pos.x, local_pos.y + planet.pos.y},
            p_mass,
            planet.v + rotation_speed,
            planet.a
        );
        _add_particle(p);
    }
}

void Model::add_water_to_planet(std::string planet_name, double p_mass, double h, double d) {
    const Planet planet = planets_.at(planet_name);
    _add_water_to_planet(planet, p_mass, h, d);
}

void Model::add_particle_planet(Planet planet, double d) {
    double r = planet.r;
    planet.r = 0;

    _add_water_to_planet(planet, planet.m, r, d);
}

void Model::add_particle_satellite(std::string planet_name, const Planet &satellite, double d, double velocity) {
    const Planet& planet = planets_[planet_name];

    double R = glm::distance(planet.pos, satellite.pos);

    std::vector<glm::dvec2> points = gen_hex_ring_points(0.0, satellite.r, d);

    std::cout << points.size() << std::endl;

    assert(!points.empty());

    double p_mass = satellite.m / points.size();

    for (const glm::dvec2 &local_pos : points) {
        double radius = glm::length(local_pos);
        double angle = std::atan2(local_pos.y, local_pos.x);

        glm::dvec2 rotation_speed = {
            -satellite.av * radius * std::sin(angle),
             satellite.av * radius * std::cos(angle)
        };

        glm::dvec2 pos {
            local_pos.x + satellite.pos.x,
            local_pos.y + satellite.pos.y
        };

        glm::dvec2 r = pos - planet.pos;
        glm::dvec2 n = glm::normalize(r);

        double dist = glm::length(r);

        Particle p(
            pos,
            p_mass,
            velocity * (dist / R) * glm::dvec2{ -n.y, n.x } + rotation_speed,
            -G * planet.m * n / (dist * dist)
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

    for (auto &pr: planets_) {
        if (pr.second.m == 0) continue;

        glm::dvec2 a = forces.at(pr.first) / pr.second.m;

        glm::dvec2 new_pos = pr.second.pos +  pr.second.v * time_step_ + 0.5 *  pr.second.a * (time_step_ * time_step_);
        glm::dvec2 new_v =  pr.second.v + 0.5 * (pr.second.a + a) * time_step_;

        // TODO: (обработка столкновений частиц и планет) не совсем правда, вынести в CollisionsSolver
        for (int i = 0; i < p_.size(); i++) {
            double dist = glm::distance(p_[i].pos, new_pos);

            if (pr.second.r >= dist) {
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

        // TODO: (обработка столкновений частиц и планет) не совсем правда, вынести в CollisionsSolver
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
    if (initial_p_.size() == 0) {
        initial_p_ = p_;
        initial_bc_ = _get_barycenter();
    }

    auto [p_forces, v_forces] = liquids_solver.calc_forces(p_);
    auto [e_forces, ed_forces] = elastic_solids_solver.calc_forces(p_);
    auto [g_forces, pl_forces] = gravity_solver.calc_forces(p_, planets_);

    std::vector<glm::dvec2> pt_forces(p_.size(), {0.0, 0.0});
    #pragma omp parallel for
    for (int i = 0; i < p_.size(); i++) {
        pt_forces[i] = e_forces[i] + ed_forces[i] + p_forces[i] + v_forces[i] + g_forces[i];
    }

    // логирование
    if (debug_) {
        double p_a = 0.0;
        double v_a = 0.0;
        double e_a = 0.0;
        double ed_a = 0.0;
        double g_a = 0.0;

        for (int i = 0; i < p_.size(); i++) {
            p_a += glm::length(p_forces[i]);
            v_a += glm::length(v_forces[i]);
            e_a += glm::length(e_forces[i]);
            ed_a += glm::length(ed_forces[i]);
            g_a += glm::length(g_forces[i]);
        }
        p_a /= p_.size();
        v_a /= p_.size();
        e_a /= p_.size();
        ed_a /= p_.size();
        g_a /= p_.size();

        glm::dvec2 bc = _get_barycenter();
        auto [ s_own_av, s_av ] = _get_angular_velocity(bc, "Нептун");

        std::cout 
            << "p=" << p_a << ' ' 
            << "v=" << v_a << ' ' 
            << "e=" << e_a << ' ' 
            << "ed=" << ed_a << ' ' 
            << "g=" << g_a << std::endl;

        std::cout << "own_angular_velocity=" << s_own_av << ' ' << "angular_velocity=" << s_av << std::endl;

        if (!debug_csv_.is_open()) {
            debug_csv_.open("debug_log.csv");
            debug_csv_ << "step,p,v,e,ed,g,own_angular_velocity,angular_velocity\n";
        }
        debug_csv_ << step_counter_ << ','
                << p_a << ','
                << v_a << ','
                << e_a << ','
                << ed_a << ','
                << g_a << ','
                << s_own_av << ','
                << s_av << '\n';
        debug_csv_.flush();
        step_counter_++;
    }

    _planets_next_step(pl_forces);
    _particles_next_step(pt_forces);
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

    glm::dvec2 bc = _get_barycenter();

    for (int i = 0; i < p_.size(); i++) {
        glm::dvec2 r = p_[i].pos - bc;
        glm::dvec2 r0 = initial_p_[i].pos - initial_bc_;

        // TODO: визуальное усиление деформаций (нужно посчитать угол поворота относительно init)
        frame.particles.push_back(
            { ((r0 + 1.0 * (r - r0)) * planet_scale_ + bc) / scale_ }
        );
    }

    return frame;
}