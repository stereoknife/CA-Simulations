#include "scenefluid.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>

SceneFluid::SceneFluid() {
    widget = new WidgetFluid();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
}


SceneFluid::~SceneFluid() {
    if (widget)     delete widget;
    if (shader)     delete shader;
    if (vaoFloor)   delete vaoFloor;
    if (vaoSphereH) delete vaoSphereH;
    if (vaoSphereL) delete vaoSphereL;
    if (vaoCube)    delete vaoCube;
    if (fGravity)   delete fGravity;
    for (auto n : neighbours) {
        if (n) delete n;
    }
}


void SceneFluid::initialize() {
    // load shader
    shader = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");

    // create floor VAO
    Model quad = Model::createQuad();
    vaoFloor = glutils::createVAO(shader, &quad, buffers);
    glutils::checkGLError();

    // create particle VAOs
    Model sphereLowres = Model::createIcosphere(1);
    vaoSphereL = glutils::createVAO(shader, &sphereLowres, buffers);
    numFacesSphereL = sphereLowres.numFaces();
    glutils::checkGLError();

    // create sphere VAO
    Model sphere = Model::createIcosphere(3);
    vaoSphereH = glutils::createVAO(shader, &sphere, buffers);
    numFacesSphereH = sphere.numFaces();
    glutils::checkGLError();

    // create box VAO
    Model cube = Model::createCube();
    vaoCube = glutils::createVAO(shader, &cube, buffers);
    glutils::checkGLError();

    // create forces
    //fGravity = new ForceConstAcceleration();
    //system.addForce(fGravity);

    // scene
    box_size = 50.0;
    fountainPos = Vec3(0, 80, 0);
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);
    colliderNorth.setPlane(Vec3(0, 0, 1), box_size);
    colliderSouth.setPlane(Vec3(0, 0, -1), box_size);
    colliderEast.setPlane(Vec3(1, 0, 0), box_size);
    colliderWest.setPlane(Vec3(-1, 0, 0), box_size);
    colliderSphere.setCenter(Vec3(0,0,0));
    colliderSphere.setRadius(20);

    draw_walls = false;

    num_x = num_y = num_z = 10;
    num_total = num_x * num_y * num_z;
    densities.resize(num_total);
    pressures.resize(num_total);
    neighbours.resize(num_total);
    for (int i = 0; i < num_total; ++i) {
        neighbours[i] = new std::vector<Particle*>();
    }
    //updateSimParams();
}


void SceneFluid::reset()
{
    // update values from UI
    updateSimParams();

    // erase all particles
    //fGravity->clearInfluencedParticles();
    system.deleteParticles();
    deadParticles.clear();

    const double box_size = 50;
    double size_x = box_size;
    double size_y = box_size;
    double size_z = box_size;
    double step = particle_size * 4;

    for(int i = 0; i < num_z; ++i) {
        for (int j = 0; j < num_y; ++j) {
            for (int k = 0; k < num_x; ++k) {
                int idx = i * num_z * num_z + j * num_y + k;

                // TODO: you can play here with different start positions and/or fixed particles
                double tx = 1 + k*step - box_size;// * (k > numParticlesX/2 ? 1.0 : -1.0);
                double ty = 1 + j*step;
                double tz = 1 + i*step - box_size;
                Vec3 pos = Vec3(tx, ty, tz);

                //std::cout << pos << std::endl;

                Particle* p = new Particle();
                p->id = idx;
                p->pos = pos;
                p->prevPos = pos;
                p->vel = Vec3(0,0,0);
                p->mass = particle_mass;
                p->color = Vec3(153/255.0, 217/255.0, 234/255.0);
                p->radius = particle_size;
                //p->life = maxParticleLife;

                system.addParticle(p);
                //fGravity->addInfluencedParticle(p);
            }
        }
    }
}


void SceneFluid::updateSimParams()
{
    // get gravity from UI and update force
    //double g = widget->getGravity();
    //fGravity->setAcceleration(Vec3(0, -g, 0));

    // get other relevant UI values and update simulation params
    kBounce = 0.5;
    kFriction = 0.1;
    maxParticleLife = 10.0;
    emitRate = 100;

    colour = (ColourMode)widget->getColour();
    draw_walls = widget->getDrawColliders();
    ref_density = widget->getReferenceDensity();
    particle_mass = widget->getParticleMass();
    dyn_viscosity = widget->getViscosity();
    c = widget->getSpeedOfSound();
    kernel_size = widget->getKernelSize();
    neighbourhood_size = widget->getNeighbourhoodSize();
    particle_size = widget->getParticleSize();
    for (auto n : neighbours) {
        n->resize(neighbourhood_size);
    }
}


void SceneFluid::paint(const Camera& camera) {

    QOpenGLFunctions* glFuncs = nullptr;
    glFuncs = QOpenGLContext::currentContext()->functions();

    shader->bind();

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shader->setUniformValue("ProjMatrix", camProj);
    shader->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 1;
    const QVector3D lightPosWorld[numLights] = {QVector3D(100,500,100)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView.map(lightPosWorld[i]);  // map = matrix * vector
    }
    shader->setUniformValue("numLights", numLights);
    shader->setUniformValueArray("lightPos", lightPosCam, numLights);
    shader->setUniformValueArray("lightColor", lightColor, numLights);

    // draw floor
    vaoFloor->bind();
    QMatrix4x4 modelMat;
    modelMat.scale(100, 1, 100);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw walls
    if (draw_walls) {
        modelMat = QMatrix4x4();
        modelMat.rotate(90.0, QVector3D(1, 0, 0));
        modelMat.translate(0, -box_size, 0);
        modelMat.scale(100, 1, 100);
        shader->setUniformValue("ModelMatrix", modelMat);
        glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        modelMat = QMatrix4x4();
        modelMat.rotate(90.0, QVector3D(1, 0, 0));
        modelMat.translate(0, box_size, 0);
        modelMat.scale(100, 1, 100);
        shader->setUniformValue("ModelMatrix", modelMat);
        glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        modelMat = QMatrix4x4();
        modelMat.rotate(90.0, QVector3D(0, 0, 1));
        modelMat.translate(0, -box_size, 0);
        modelMat.scale(100, 1, 100);
        shader->setUniformValue("ModelMatrix", modelMat);
        glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        modelMat = QMatrix4x4();
        modelMat.rotate(90.0, QVector3D(0, 0, 1));
        modelMat.translate(0, box_size, 0);
        modelMat.scale(100, 1, 100);
        shader->setUniformValue("ModelMatrix", modelMat);
        glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    }

    // draw the particles
    vaoSphereL->bind();
    for (const Particle* particle : system.getParticles()) {
        Vec3   p = particle->pos;
        Vec3   c = particle->color;
        double r = particle->radius;

        modelMat = QMatrix4x4();
        modelMat.translate(p[0], p[1], p[2]);
        modelMat.scale(r);
        shader->setUniformValue("ModelMatrix", modelMat);

        shader->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shader->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereL, GL_UNSIGNED_INT, 0);
    }

    /*
    // draw sphere
    vaoSphereH->bind();
    Vec3 cc = colliderSphere.getCenter();
    modelMat = QMatrix4x4();
    modelMat.translate(cc[0], cc[1], cc[2]);
    modelMat.scale(colliderSphere.getRadius());
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.4f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    // draw box
    vaoCube->bind();
    cc = colliderBox.getCenter();
    Vec3 hs = 0.5*colliderBox.getSize();
    modelMat = QMatrix4x4();
    modelMat.translate(cc[0], cc[1], cc[2]);
    modelMat.scale(hs[0], hs[1], hs[2]);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.4f, 0.8f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);
    vaoCube->release();
    shader->release();
*/
}

double poly6(double r, double h) {
    if (r < 0 || h < r) return 0;
    return (315/(64*M_PI*std::pow(h,9)))*std::pow((h*h)-(r*r),3);
}

double spiky(double r, double h) {
    if (r < 0 || h < r) return 0;
    return (15/(M_PI*std::pow(h, 6)))*std::pow(h-r, 3);
}

Vec3 spiky_grad(Vec3 d, double r, double h) {
    return -d * (45/(M_PI*std::pow(h,6)*r)) * std::pow((h-r),2);
}

double viscosity_lapl(double r, double h){
    return (45/(M_PI * std::pow(h, 5)))*(1-(r/h));
}


void SceneFluid::update(double dt) {
    // SPH
    densities.clear();
    pressures.clear();

    // 1. Find neighbours for each particle and store in a list
    std::vector<std::tuple<double, Particle*>> particles((num_x * num_y * num_z));
    for (Particle* p : system.getParticles()) {
        particles.clear();
        for (Particle* o : system.getParticles()) {
            if(p->id == o->id) continue;
            particles.emplace_back((p->pos - o->pos).squaredNorm(), o);
        }
        std::sort(particles.begin(), particles.end(), [](auto a, auto b){return std::get<0>(a) < std::get<0>(b);});
        for (int i = 0; i < neighbourhood_size; ++i) {
            (*neighbours[p->id])[i] = std::get<1>(particles[i]);
        }
    }

    // 2. Calculate density for each particle
    for (Particle* p : system.getParticles()) {
        double d = p->mass * poly6(0, kernel_size);
        for(auto o : *neighbours[p->id]){
            //std::cout << "Dist: " << (p->pos-o->pos).norm() << ", Range: " << kernel_size << std::endl;
            auto dd = o->mass * poly6((p->pos-o->pos).norm(), kernel_size);
            //if (p->id == 555) std::cout << dd << std::endl;
            d += dd;
        }
        densities[p->id] = d;
    }

    // 3. Calculate pressure for each particle
    for (Particle* p : system.getParticles()) {
        auto pressure = c * c * (densities[p->id] - ref_density);
        pressures[p->id] = pressure;
        //if (p->id == 555) std::cout << "Density: " << densities[p->id] << ", Pressure: " << pressures[p->id] << std::endl;
    }

    // 4. Calculate all type of accelerations for each particle
    for (Particle* i : system.getParticles()) {
        i->force = Vec3(0,0,0);

        Vec3 pressure_acc(0,0,0), viscosity_acc(0,0,0), interactive_acc(0,0,0), gravity_acc(0,-9.8,0);

        //std::cout << "Density: " << densities[i->id] << ", Pressure: " << pressures[i->id] << std::endl;

        for (auto j : *neighbours[i->id]) {
            double pij = j->mass * ((pressures[i->id]/std::pow(densities[i->id],2)) + (pressures[j->id]/std::pow(densities[j->id],2)));
            Vec3 d = j->pos - i->pos;
            double r = d.norm();

            pressure_acc += pij * spiky_grad(d, d.norm(), kernel_size);

            Vec3 vij = dyn_viscosity * j->mass * ((j->vel - i->vel)/(densities[i->id]*densities[j->id]));
            viscosity_acc += vij * viscosity_lapl(r, kernel_size);
        }

        //std::cout << "Pressure: " << pressure_acc.transpose() << ",  Viscosity: " << viscosity_acc.transpose() << std::endl;
        Vec3 total_acc = pressure_acc + viscosity_acc + interactive_acc + gravity_acc;
        i->force = total_acc * i->mass;
    }

    // 5. Find new velocities and positions by using the same integration method as before
    Vecd ppos = system.getPositions();
    integrator.step(system, dt);
    system.setPreviousPositions(ppos);

    // collisions
    Collision colInfo;
    for (Particle* p : system.getParticles()) {
        if (colliderFloor.testCollision(p, colInfo)) {
            colliderFloor.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderNorth.testCollision(p, colInfo)) {
            colliderNorth.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderSouth.testCollision(p, colInfo)) {
            colliderSouth.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderEast.testCollision(p, colInfo)) {
            colliderEast.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderWest.testCollision(p, colInfo)) {
            colliderWest.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderSphere.testCollision(p, colInfo)) {
            colliderSphere.resolveCollision(p, colInfo, kBounce, kFriction);
        }
    }

    //6. Colour particles
    switch(colour) {
    case ColourMode::Neighbourhood: {
        auto p = system.getParticles()[555];
        p->color = Vec3(0,1,0);
        for (auto o : *neighbours[555]) {
            double norm = (p->pos - o->pos).norm();
            if (norm <= neighbourhood_size) {
                o->color = Vec3(1-norm/neighbourhood_size, 0, 0);
            }
        }
    } break;
    case ColourMode::Density:
        for (auto p : system.getParticles()) {
            auto d = densities[p->id];
            p->color = Vec3((d-1)*10,0,1-(d-1)*10);
        }
        break;
    case ColourMode::Pressure:
        for (auto p : system.getParticles()) {
            auto d = densities[p->id];
            if (d < 0) {
                p->color = Vec3(1+d/10,1+d/10,1);
            } else {
                p->color = Vec3(1,1-d/10,1-d/10);
            }
        }
        break;
    case ColourMode::None:
    default:
        break;
    }
}

void SceneFluid::mousePressed(const QMouseEvent* e, const Camera&)
{
    mouseX = e->pos().x();
    mouseY = e->pos().y();
}

void SceneFluid::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - mouseX;
    int dy = e->pos().y() - mouseY;
    mouseX = e->pos().x();
    mouseY = e->pos().y();

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    // example
    if (e->buttons() & Qt::RightButton) {
        if (!e->modifiers()) {
            // move fountain
            fountainPos += disp;
        }
        /*
        else if (e->modifiers() & Qt::ShiftModifier){
            // move box
            colliderBox.setFromCenterSize(
                        colliderBox.getCenter() + disp,
                        colliderBox.getSize());
        }
        */
    }
}
