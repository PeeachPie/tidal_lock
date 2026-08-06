#pragma once

#include <vector>
#include <iostream>
#include "particle.hpp"

std::vector<std::vector<int>> get_neighborhoods(const std::vector<Particle> &p, double search_radius);