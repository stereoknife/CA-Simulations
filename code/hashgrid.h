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
    std::span<unsigned int> particle_span;

public:
    HashGrid(unsigned int num_particles, double radius);
    void construct(std::vector<Particle*> particles);
    std::span<unsigned int> neighbours(Vec3 pos);


private:
    inline unsigned int hash(int x, int y, int z) {
        auto h = (x * 92837111) ^ (y * 689287499) ^ (z * 283923481);
        return std::abs(h) % num_particles;
    }

    inline unsigned int coords_to_hash(double x, double y, double z) {
        return hash(std::floor(x / radius), std::floor(y / radius), std::floor(z / radius));
    };

    inline unsigned int pos_to_hash(Vec3 pos) {
        return hash(pos.x(), pos.y(), pos.z());
    };

};

#endif // HASHGRID_H
