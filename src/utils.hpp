#pragma once

#include <random>
#include <glm/glm.hpp>
#include <cmath>
#include <algorithm>

#include "constants.hpp"

double gen_angle();

double gen_radius(double r, double h);

double point_to_segment_dist(const glm::dvec2 &a, const glm::dvec2 &b, const glm::dvec2 &p);

void svd(const glm::dmat2& A, glm::dmat2& U, glm::dmat2& D, glm::dmat2& V);

glm::dmat2 extract_rotation(const glm::dmat2& A);

glm::dmat2 moore_penrose_inverse(const glm::dmat2& A, double eps=1e-12);

glm::dmat2 tensor_product(const glm::dvec2 &a, const glm::dvec2 &b);

double trace(const glm::dmat2 &matrix);

double calc_orbital_velocity(double M, double R);