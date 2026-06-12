#include "utils.hpp"
#include <iostream>

std::mt19937 gen(42);

double gen_angle() {
    std::uniform_real_distribution<> dis(0.0, 2.0 * PI);

    return dis(gen);
}

double gen_h(double h) {
    std::uniform_real_distribution<> dis(0.0, h);

    return dis(gen);
}

double cross_prod(const glm::dvec2 &a, const glm::dvec2 &b) {
    return a.x * b.y - b.x * a.y;
}

double point_to_segment_dist(glm::dvec2 a, glm::dvec2 b, glm::dvec2 p) {
    glm::dvec2 ap = p - a;
    glm::dvec2 bp = p - b;
    glm::dvec2 ab = b - a;
    glm::dvec2 ba = a - b;

    if (glm::dot(ap, ab) >= 0 && glm::dot(bp, ba) >= 0) {
        return 0.5 * glm::abs(cross_prod(ap, ab)) / glm::length(ab);
    } else {
        return glm::min(glm::distance(a, p), glm::distance(b, p));
    }
}