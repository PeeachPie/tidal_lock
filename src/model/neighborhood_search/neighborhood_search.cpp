#include "neighborhood_search.hpp"

std::unordered_map<int, std::unordered_map<int, std::vector<int>>> get_particle_grid(const std::vector<Particle> &p, double search_radius) {
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> particle_grid;

    for (int i = 0; i < p.size(); i++) {
        particle_grid[
            std::floor(p[i].pos.x / search_radius)
        ][
            std::floor(p[i].pos.y / search_radius)
        ].push_back(i);
    }

    return particle_grid;
}

std::vector<int> get_particle_neighborhood(int pi, const std::vector<Particle> &p, double search_radius, const std::unordered_map<int, std::unordered_map<int, std::vector<int>>> &particle_grid) {
    std::vector<int> neighborhood;

    int x = floor(p[pi].pos.x / search_radius);
    int y = floor(p[pi].pos.y / search_radius);

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (particle_grid.count(x + i) && particle_grid.at(x + i).count(y + j)) {
                for (int k: particle_grid.at(x + i).at(y + j)) {
                    glm::dvec2 r = p[pi].pos - p[k].pos;
                    if (glm::dot(r, r) <= search_radius * search_radius) {
                        neighborhood.push_back(k);
                    }
                }
            }
        }
    }

    return neighborhood;
}

std::vector<std::vector<int>> get_neighborhoods(const std::vector<Particle> &p, double search_radius) {
    std::vector<std::vector<int>> neighborhoods(p.size());
    
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> particle_grid = get_particle_grid(p, search_radius);
    
    #pragma omp parallel for
    for (int i = 0; i < p.size(); i++) {
        neighborhoods[i] = get_particle_neighborhood(i, p, search_radius, particle_grid);
        // для проверки кол-ва соседей
        // std::cout << neighborhoods[i].size() << '\n';
    }

    return neighborhoods;
}