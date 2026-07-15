#pragma once

#include <vector>
#include <glm/glm.hpp>

struct Particle {
    glm::dvec2 pos;
    glm::dvec2 v;
    glm::dvec2 a;

    Particle(): pos(0), v(0), a(0) {};
    Particle(glm::dvec2 pos, glm::dvec2 v=glm::dvec2(0, 0), glm::dvec2 a=glm::dvec2(0, 0)): 
        pos(pos), v(v), a(a) {};
};

struct Planet: Particle {
    double r;
    double m;
    double av;

    Planet(): Particle(), r(0), m(0), av(0) {};
    Planet(glm::dvec2 pos, double r, double m, glm::dvec2 v=glm::dvec2(0, 0), glm::dvec2 a=glm::dvec2(0, 0), double av=0.0):
        Particle(pos, v, a), r(r), m(m), av(av) {};
};

struct ParticleFrame {
    glm::dvec2 pos;
};

struct PlanetFrame {
    glm::dvec2 pos;
    double r;
};

struct Frame {
    std::vector<PlanetFrame> planets;
    std::vector<ParticleFrame> particles;

    float particle_r;
};