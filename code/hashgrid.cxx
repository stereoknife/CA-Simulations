#include <span>
#include "hashgrid.h"

HashGrid::HashGrid(unsigned int num_particles, double radius):
    num_particles{num_particles},
    radius{radius},
    start_indices{num_particles + 1},
    particle_indices{num_particles},
    particle_span{particle_indices}
{}


void HashGrid::construct(std::vector<Particle*> particles) {
    start_indices.clear();
    particle_indices.clear();

    for (Particle* p : particles) {
        auto pos = p->pos;
        auto i = coords_to_hash(pos.x(), pos.y(), pos.z());
        start_indices[i]++;
    }

    for (unsigned int i = 0, sum = 0; i < start_indices.size(); ++i) {
        sum += start_indices[i];
        start_indices[i] = sum;
    }

    for (unsigned int p = 0; p < particles.size(); ++p) {
        auto i = pos_to_hash(particles[p]->pos);
        i = --start_indices[i];
        particle_indices[i] = p;
    }
}

std::span<unsigned int> HashGrid::neighbours(Vec3 pos) {
    auto i = pos_to_hash(pos);
    auto si = start_indices[i];
    auto ei = start_indices[i + 1];
    return particle_span.subspan(si, ei);
}
