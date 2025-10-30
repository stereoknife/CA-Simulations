#ifndef SCENEFLUID_H
#define SCENEFLUID_H

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <list>
#include "scene.h"
#include "widgetfluid.h"
#include "particlesystem.h"
#include "integrators.h"
#include "colliders.h"

class SceneFluid : public Scene
{
    Q_OBJECT

    enum struct ColourMode {
        None, Neighbourhood, Density, Pressure
    };

public:
    SceneFluid();
    virtual ~SceneFluid();

    virtual void initialize();
    virtual void reset();
    virtual void update(double dt);
    virtual void paint(const Camera& cam);

    virtual void mousePressed(const QMouseEvent* e, const Camera& cam);
    virtual void mouseMoved(const QMouseEvent* e, const Camera& cam);

    virtual void getSceneBounds(Vec3& bmin, Vec3& bmax) {
        bmin = Vec3(-110, -10, -110);
        bmax = Vec3( 110, 100,  110);
    }
    virtual unsigned int getNumParticles() { return system.getNumParticles(); }

    virtual QWidget* sceneUI() { return widget; }

public slots:
    void updateSimParams();

protected:
    WidgetFluid* widget = nullptr;

    QOpenGLShaderProgram* shader = nullptr;
    QOpenGLVertexArrayObject* vaoSphereL = nullptr;
    QOpenGLVertexArrayObject* vaoSphereH = nullptr;
    QOpenGLVertexArrayObject* vaoCube    = nullptr;
    QOpenGLVertexArrayObject* vaoFloor   = nullptr;
    unsigned int numFacesSphereL = 0, numFacesSphereH = 0;

    IntegratorEuler integrator;
    ParticleSystem system;
    std::list<Particle*> deadParticles;
    ForceConstAcceleration* fGravity;

    ColliderPlane colliderFloor, colliderNorth, colliderSouth, colliderEast, colliderWest;
    ColliderSphere colliderSphere;

    double kBounce, kFriction;
    double emitRate;
    double maxParticleLife;

    Vec3 fountainPos;
    int mouseX, mouseY;

    bool draw_walls;
    double box_size;
    int num_x, num_y, num_z, num_total;

    double ref_density = 10;
    double particle_mass = 1;
    double dyn_viscosity = 0.001;
    double c = 1;
    double kernel_size = 10;
    int neighbourhood_size = 50;
    double particle_size;
    ColourMode colour;

    std::vector<std::vector<Particle*>*> neighbours;
    std::vector<double> densities, pressures;
};

#endif // SCENEFLUID_H
