#include "forces.h"

void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        // F = m * a
        p->force += p->mass * acceleration;
    }
}

void ForceDrag::apply() {
    for (Particle* p : particles) {
        p->force -= getLinearCoefficient() * p->vel;
        p->force -= getQuadraticCoefficient() * p->vel.norm() * p->vel;
    }
}

void ForceSpring::apply() {
    if (particles.size() < 2) return;
    Particle* p1 = getParticle1();
    Particle* p2 = getParticle2();

    auto distv = p2->pos - p1->pos;
    auto dist = distv.norm();

    if (dist < 1e-8) return;

    auto distn = distv.normalized();

    auto f1 = (getSpringConstant() * (dist - getRestLength()) + getDampingCoeff() * (p2->vel - p1->vel).dot(distn)) * distn;

    p1->force += f1;
    p2->force -= f1;
}


void ForceGravitation::apply() {
    for (Particle* p : particles) {
        auto at = getAttractor();
        auto dist = at->pos - p->pos;
        double dd = dist.squaredNorm();

        auto g = ((getConstant() * at->mass * p->mass)/dd)*dist.normalized();

        double exponent = -a * (dd/(b*b));
        double denom = 1.0 + std::exp(exponent);
        double smoothing = ((2.0/denom)-1.0);
        //smoothing = 1;

        p->force += g * smoothing;
    }
}
