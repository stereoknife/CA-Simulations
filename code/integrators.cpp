#include "integrators.h"


void IntegratorEuler::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd dx = system.getDerivative();
    Vecd x1 = x0 + dt*dx;
    system.setState(x1);
    system.setTime(t0+dt);
    system.updateForces();
}


void IntegratorSymplecticEuler::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd v0 = system.getVelocities();
    Vecd a0 = system.getAccelerations();
    Vecd v1 = v0 + a0 * dt;
    Vecd p0 = system.getPositions();
    Vecd p1 = p0 + v1 * dt;
    system.setPositions(p1);
    system.setVelocities(v1);
    system.setTime(t0 + dt);
    system.updateForces();
}

void IntegratorMidpoint::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd d0 = system.getDerivative();

    // Get derivative at t1/2
    Vecd x05 = x0 + dt * 0.5 * d0;
    system.setState(x05);
    system.setTime(t0 + dt * 0.5);
    system.updateForces();
    Vecd d05 = system.getDerivative();

    // Calculate x1 using the derivative at t1/2
    Vecd x1 = x0 + dt * d05;
    system.setState(x1);
    system.setTime(t0 + dt);
    system.updateForces();
}

void IntegratorRK2::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd k1 = system.getDerivative();

    // Get derivative at t1
    Vecd xk1 = x0 + dt * k1;
    system.setState(xk1);
    system.setTime(t0 + dt);
    system.updateForces();
    Vecd k2 = system.getDerivative();

    // Calculate x1 using the derivative at t1/2
    // Time was already advanced so no need to do it again
    Vecd x1 = x0 + dt * 0.5 * (k1 + k2);
    system.setState(x1);
    system.updateForces();
}


void IntegratorRK4::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd k1 = system.getDerivative();
    double hdt = dt / 2;

    Vecd xk1 = x0 + hdt * k1;
    system.setState(xk1);
    system.setTime(t0 + hdt);
    system.updateForces();
    Vecd k2 = system.getDerivative();

    Vecd xk2 = x0 + hdt * k2;
    system.setState(xk2);
    system.updateForces();
    Vecd k3 = system.getDerivative();

    Vecd xk3 = x0 + dt * k3;
    system.setState(xk3);
    system.setTime(t0 + dt);
    system.updateForces();
    Vecd k4 = system.getDerivative();

    // Calculate x1 using the derivative at t1/2
    // Time was already advanced so no need to do it again
    Vecd x1 = x0 + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    system.setState(x1);
    system.updateForces();
}


void IntegratorVerlet::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd p0 = system.getPositions();
    Vecd v0 = system.getVelocities();
    Vecd a0 = system.getAccelerations();
    Vecd pn1 = system.getPreviousPositions();

    double k = 0.99;
    Vecd p1 = p0 + k * (p0 - pn1) + dt * dt * a0;
    Vecd v1 = (p1 -  p0) /  dt;

    system.setPreviousPositions(p0);
    system.setPositions(p1);
    system.setVelocities(v1);
    system.setTime(t0 + dt);

    system.updateForces();
}
