#pragma once

#include <vector>
#include <iostream>
#include <cmath> 
#include <unordered_map>
#include "particle.hpp"

std::vector<std::vector<int>> get_neighborhoods(const std::vector<Particle> &p, double search_radius);