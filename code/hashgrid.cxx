#include <iostream>
#include <span>
#include "hashgrid.h"

HashGrid::HashGrid(unsigned int num_particles, double radius):
    num_particles{num_particles * 2},
    radius{radius},
    start_indices{num_particles + 1},
    particle_indices{num_particles}
{}


void HashGrid::construct(std::vector<Particle*> particles) {
    fill(start_indices.begin(), start_indices.end(), 0);
    fill(particle_indices.begin(), particle_indices.end(), 0);

    for (Particle* p : particles) {
        auto i = pos_to_hash(p->pos);
        start_indices[i]++;
    }

    for (unsigned int i = 0, sum = 0; i < start_indices.size(); ++i) {
        sum += start_indices[i];
        start_indices[i] = sum;
    }

    for (unsigned int p = 0; p < particles.size(); ++p) {
        auto i = pos_to_hash(particles[p]->pos);
        auto j = --start_indices[i];
        particle_indices[j] = p;
    }

}

std::span<unsigned int> HashGrid::neighbours(Vec3 pos) {
    auto i = pos_to_hash(pos);
    auto si = start_indices[i];
    auto ei = start_indices[i + 1];
    auto sp = std::span{particle_indices}.subspan(si, ei - si);
    return sp;
}

void HashGrid::all_cells(std::vector<std::span<unsigned int>>& output) {
    std::span<unsigned int> s{particle_indices};
    for (int i = 0; i < start_indices.size() - 1; ++i) {
        auto si = start_indices[i];
        auto ei = start_indices[i + 1];
        if (ei > si) output.push_back(s.subspan(si, ei - si));
    }
}
