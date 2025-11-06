#ifndef HASHGRID_H
#define HASHGRID_H

#include "particle.h"
#include <vector>
#include <span>

class HashGrid
{
private:
    unsigned int num_particles;
    double radius;
    std::vector<unsigned int> start_indices;
    std::vector<unsigned int> particle_indices;

public:
    HashGrid(unsigned int num_particles, double radius);
    void construct(std::vector<Particle*> particles);
    std::span<unsigned int> neighbours(Vec3 pos);
    void all_cells(std::vector<std::span<unsigned int>>& output);

    inline void resize(int size) {
        num_particles = size;
        start_indices.resize(num_particles + 1);
        particle_indices.resize(num_particles);
    };

    inline void set_radius(double r) {
        radius = r;
    };


private:
    // Hash function from 10 Minute Physics's video
    inline unsigned int hash(int x, int y, int z) {
        auto h = (x * 92837111) ^ (y * 689287499) ^ (z * 283923481);
        return std::abs(h) % num_particles;
    }

    inline unsigned int pos_to_hash(Vec3 pos) {
        return coords_to_hash(pos.x(), pos.y(), pos.z());
    };

    inline unsigned int coords_to_hash(double x, double y, double z) {
        return hash(std::floor(x / radius), std::floor(y / radius), std::floor(z / radius));
    };

};

#endif // HASHGRID_H
