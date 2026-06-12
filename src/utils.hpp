#pragma once

#include <random>
#include <glm/glm.hpp>
#include "constants.hpp"

double gen_angle();

double gen_h(double h);

double point_to_segment_dist(glm::dvec2 a, glm::dvec2 b, glm::dvec2 p);