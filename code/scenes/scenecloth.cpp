#include "scenecloth.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLBuffer>


SceneCloth::SceneCloth() {
    widget = new WidgetCloth();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
    connect(widget, SIGNAL(freeAnchors()), this, SLOT(freeAnchors()));
}

SceneCloth::~SceneCloth() {
    if (widget)      delete widget;
    if (shaderPhong) delete shaderPhong;
    if (vaoSphereS)  delete vaoSphereS;
    if (vaoSphereL)  delete vaoSphereL;
    if (vaoCube)     delete vaoCube;
    if (vaoMesh)     delete vaoMesh;
    if (vboMesh)     delete vboMesh;
    if (iboMesh)     delete iboMesh;

    system.deleteParticles();
    if (fGravity)  delete fGravity;
    for (ForceSpring* f : springsStretch) delete f;
    for (ForceSpring* f : springsShear) delete f;
    for (ForceSpring* f : springsBend) delete f;
}

void SceneCloth::initialize() {

    // load shaders
    shaderPhong = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");
    shaderCloth = glutils::loadShaderProgram(":/shaders/cloth.vert", ":/shaders/cloth.geom", ":/shaders/cloth.frag");

    // create sphere VAOs
    Model sphere = Model::createIcosphere(3);
    vaoSphereL = glutils::createVAO(shaderPhong, &sphere, buffers);
    numFacesSphereL = sphere.numFaces();
    glutils::checkGLError();

    sphere = Model::createIcosphere(1);
    vaoSphereS = glutils::createVAO(shaderPhong, &sphere, buffers);
    numFacesSphereS = sphere.numFaces();
    glutils::checkGLError();

    // create cube VAO
    Model cube = Model::createCube();
    vaoCube = glutils::createVAO(shaderPhong, &cube, buffers);
    glutils::checkGLError();


    // create cloth mesh VAO
    vaoMesh = new QOpenGLVertexArrayObject();
    vaoMesh->create();
    vaoMesh->bind();
    vboMesh = new QOpenGLBuffer(QOpenGLBuffer::Type::VertexBuffer);
    vboMesh->create();
    vboMesh->bind();
    vboMesh->setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
    vboMesh->allocate(1000*1000*3*3*sizeof(float)); // sync with widget max particles
    shaderCloth->setAttributeBuffer("vertex", GL_FLOAT, 0, 3, 0);
    shaderCloth->enableAttributeArray("vertex");
    iboMesh = new QOpenGLBuffer(QOpenGLBuffer::Type::IndexBuffer);
    iboMesh->create();
    iboMesh->bind();
    iboMesh->setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
    iboMesh->allocate(1000*1000*2*3*sizeof(unsigned int));
    vaoMesh->release();

    // create gravity force
    fGravity = new ForceConstAcceleration();
    system.addForce(fGravity);

    // TODO: in my solution setup, these were the colliders
    colliderBall.setCenter(Vec3(40,-20,0));
    colliderBall.setRadius(30);
    //colliderCube.setFromCenterSize(Vec3(-60,30,0), Vec3(60, 40, 60));
    //colliderWalls.setFromCenterSize(Vec3(0, 0, 0), Vec3(200, 200, 200));
}

void SceneCloth::reset()
{
    // we only update numParticles on resets
    updateSimParams();

    // reset particles
    system.deleteParticles();

    // reset forces
    system.clearForces();
    fGravity->clearInfluencedParticles();
    for (ForceSpring* f : springsStretch) delete f;
    springsStretch.clear();
    for (ForceSpring* f : springsShear) delete f;
    springsShear.clear();
    for (ForceSpring* f : springsBend) delete f;
    springsBend.clear();

    // cloth props
    Vec2 dims = widget->getDimensions();
    Vec2i dimParticles = widget->getNumParticles();
    numParticlesX = dimParticles.x();
    numParticlesY = dimParticles.y();
    clothWidth  = dims[0];
    clothHeight = dims[1];
    double edgeX = dims[0]/numParticlesX;
    double edgeY = dims[1]/numParticlesY;
    particleRadius = widget->getParticleRadius();

    // create particles
    numParticles = numParticlesX * numParticlesY;
    fixedParticle = std::vector<bool>(numParticles, false);

    for (int i = 0; i < numParticlesX; i++) {
        for (int j = 0; j < numParticlesY; j++) {

            int idx = i*numParticlesY + j;

            // TODO: you can play here with different start positions and/or fixed particles
            fixedParticle[idx] = idx < numParticlesY && idx % 8 == 0;
            double tx = i*edgeX - 0.5*clothWidth;
            double ty = j*edgeY - 0.5*clothHeight;
            Vec3 pos = Vec3(ty+edgeY, 30, tx + edgeX);

            Particle* p = new Particle();
            p->id = idx;
            p->pos = pos;
            p->prevPos = pos;
            p->vel = Vec3(0,0,0);
            p->mass = 1;
            p->radius = particleRadius;
            p->color = Vec3(235/255.0, 51/255.0, 36/255.0);

            system.addParticle(p);
            fGravity->addInfluencedParticle(p);
        }
    }

    // forces: gravity
    system.addForce(fGravity);

    // TODO: create spring forces
    double ks = widget->getStiffness();
    double kd = widget->getDamping();
    auto particles = system.getParticles();
    double widthLength = clothWidth / numParticlesX;
    double heightLength = clothHeight / numParticlesY;
    double diagonalLength = sqrt(widthLength * widthLength + heightLength * heightLength);
    for (int i = 0; i < numParticlesX; ++i) {
        for (int j = 0; j < numParticlesY; ++j) {
            ForceSpring* f;

            if (i > 0) {
                f = new ForceSpring();
                f->setParticlePair(particles[(i-1)*numParticlesY + j], particles[i*numParticlesY + j]);
                f->setRestLength(widthLength);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsStretch.push_back(f);
            }

            if (j > 0) {
                f = new ForceSpring();
                f->setParticlePair(particles[i*numParticlesY + (j-1)], particles[i*numParticlesY + j]);
                f->setRestLength(heightLength);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsStretch.push_back(f);
            }

            if (i > 0 && j > 0) {
                f = new ForceSpring();
                f->setParticlePair(particles[(i-1)*numParticlesY + (j-1)], particles[i*numParticlesY + j]);
                f->setRestLength(diagonalLength);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsShear.push_back(f);
            }

            if (i > 0 && j < numParticlesX - 1) {
                f = new ForceSpring();
                f->setParticlePair(particles[(i-1)*numParticlesY + (j+1)], particles[i*numParticlesY + j]);
                f->setRestLength(diagonalLength);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsShear.push_back(f);
            }

            if (i > 1) {
                f = new ForceSpring();
                f->setParticlePair(particles[(i-2)*numParticlesY + j], particles[i*numParticlesY + j]);
                f->setRestLength(heightLength * 2);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsBend.push_back(f);
            }

            if (j > 1) {
                f = new ForceSpring();
                f->setParticlePair(particles[i*numParticlesY + (j-2)], particles[i*numParticlesY + j]);
                f->setRestLength(widthLength * 2);
                f->setSpringConstant(ks);
                f->setDampingCoeff(kd);
                system.addForce(f);
                springsBend.push_back(f);
            }
        }
    }


    // Code for PROVOT layout
    updateSprings();

    // update index buffer
    iboMesh->bind();
    numMeshIndices = (numParticlesX - 1)*(numParticlesY - 1)*2*3;
    int* indices = new int[numMeshIndices];
    int idx = 0;
    for (int i = 0; i < numParticlesX-1; i++) {
        for (int j = 0; j < numParticlesY-1; j++) {
            indices[idx  ] = i*numParticlesY + j;
            indices[idx+1] = (i+1)*numParticlesY + j;
            indices[idx+2] = i*numParticlesY + j + 1;
            indices[idx+3] = i*numParticlesY + j + 1;
            indices[idx+4] = (i+1)*numParticlesY + j;
            indices[idx+5] = (i+1)*numParticlesY + j + 1;
            idx += 6;
        }
    }
    void* bufptr = iboMesh->mapRange(0, numMeshIndices*sizeof(int),
                                     QOpenGLBuffer::RangeInvalidateBuffer | QOpenGLBuffer::RangeWrite);
    memcpy(bufptr, (void*)(indices), numMeshIndices*sizeof(int));
    iboMesh->unmap();
    iboMesh->release();
    delete[] indices;
    glutils::checkGLError();
}


void SceneCloth::updateSprings()
{
    double ks = widget->getStiffness();
    double kd = widget->getDamping();

    // here I update all ks and kd parameters.
    // idea: if you want to enable/disable a spring type, you can set ks to 0 for these
    for (ForceSpring* f : springsStretch) {
        f->setSpringConstant(ks);
        f->setDampingCoeff(kd);
    }
    for (ForceSpring* f : springsShear) {
        f->setSpringConstant(ks);
        f->setDampingCoeff(kd);
    }
    for (ForceSpring* f : springsBend) {
        f->setSpringConstant(ks);
        f->setDampingCoeff(kd);
    }
}

void SceneCloth::updateSimParams()
{
    double g = widget->getGravity();
    fGravity->setAcceleration(Vec3(0, -g, 0));

    updateSprings();

    for (Particle* p : system.getParticles()) {
        p->radius = widget->getParticleRadius();
    }

    showParticles = widget->showParticles();
}

void SceneCloth::freeAnchors()
{
    fixedParticle = std::vector<bool>(numParticles, false);
}

void SceneCloth::paint(const Camera& camera)
{
    QOpenGLFunctions* glFuncs = nullptr;
    glFuncs = QOpenGLContext::currentContext()->functions();

    shaderPhong->bind();
    shaderPhong->setUniformValue("normalSign", 1.0f);

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shaderPhong->setUniformValue("ProjMatrix", camProj);
    shaderPhong->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 1;
    const QVector3D lightPosWorld[numLights] = {QVector3D(80,80,80)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView.map(lightPosWorld[i]);
    }
    shaderPhong->setUniformValue("numLights", numLights);
    shaderPhong->setUniformValueArray("lightPos", lightPosCam, numLights);
    shaderPhong->setUniformValueArray("lightColor", lightColor, numLights);

    // draw the particle spheres
    QMatrix4x4 modelMat;
    if (showParticles) {
        vaoSphereS->bind();
        shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shaderPhong->setUniformValue("matshin", 100.f);
        for (int i = 0; i < numParticles; i++) {
            const Particle* particle = system.getParticle(i);
            Vec3   p = particle->pos;
            Vec3   c = particle->color;
            if (fixedParticle[i])      c = Vec3(63/255.0, 72/255.0, 204/255.0);
            if (i == selectedParticle) c = Vec3(1.0,0.9,0);

            modelMat = QMatrix4x4();
            modelMat.translate(p[0], p[1], p[2]);
            modelMat.scale(particle->radius);
            shaderPhong->setUniformValue("ModelMatrix", modelMat);
            shaderPhong->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
            glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereS, GL_UNSIGNED_INT, 0);
        }
    }

    // TODO: draw colliders and walls

    vaoSphereL->bind();
    Vec3 cc = colliderBall.getCenter();
    modelMat = QMatrix4x4();
    modelMat.translate(cc[0], cc[1], cc[2]);
    modelMat.scale(colliderBall.getRadius());
    shaderPhong->setUniformValue("ModelMatrix", modelMat);
    shaderPhong->setUniformValue("matdiff", 0.8f, 0.4f, 0.4f);
    shaderPhong->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shaderPhong->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereL, GL_UNSIGNED_INT, 0);

    shaderPhong->release();


    // update cloth mesh VBO coords
    vboMesh->bind();
    float* pos = new float[3*numParticles];
    for (int i = 0; i < numParticles; i++) {
        pos[3*i  ] = system.getParticle(i)->pos.x();
        pos[3*i+1] = system.getParticle(i)->pos.y();
        pos[3*i+2] = system.getParticle(i)->pos.z();
    }
    void* bufptr = vboMesh->mapRange(0, 3*numParticles*sizeof(float),
                       QOpenGLBuffer::RangeInvalidateBuffer | QOpenGLBuffer::RangeWrite);
    memcpy(bufptr, (void*)(pos), 3*numParticles*sizeof(float));
    vboMesh->unmap();
    vboMesh->release();
    delete[] pos;

    // draw mesh
    shaderCloth->bind();
    shaderCloth->setUniformValue("ProjMatrix", camProj);
    shaderCloth->setUniformValue("ViewMatrix", camView);
    shaderCloth->setUniformValue("NormalMatrix", camView.normalMatrix());
    shaderCloth->setUniformValue("matdiffFront", 0.7f, 0.0f, 0.0f);
    shaderCloth->setUniformValue("matspecFront", 1.0f, 1.0f, 1.0f);
    shaderCloth->setUniformValue("matshinFront", 100.0f);
    shaderCloth->setUniformValue("matdiffBack", 0.7f, 0.3f, 0.0f);
    shaderCloth->setUniformValue("matspecBack", 0.0f, 0.0f, 0.0f);
    shaderCloth->setUniformValue("matshinBack", 0.0f);
    shaderCloth->setUniformValue("numLights", numLights);
    shaderCloth->setUniformValueArray("lightPos", lightPosCam, numLights);
    shaderCloth->setUniformValueArray("lightColor", lightColor, numLights);
    vaoMesh->bind();
    glFuncs->glDrawElements(GL_TRIANGLES, numMeshIndices, GL_UNSIGNED_INT, 0);
    vaoMesh->release();
    shaderCloth->release();

    glutils::checkGLError();
}


void SceneCloth::update(double dt)
{
    // fixed particles: no velocity, no force acting
    for (int i = 0; i < numParticles; i++) {
        if (fixedParticle[i]) {
            Particle* p = system.getParticle(i);
            p->vel = Vec3(0,0,0);
            p->force = Vec3(0,0,0);
        }
    }

    // integration step
    Vecd ppos = system.getPositions();
    integrator.step(system, dt);
    system.setPreviousPositions(ppos);

    // user interaction
    if (selectedParticle >= 0) {
        Particle* p = system.getParticle(selectedParticle);
        // p->pos = ?; TODO: assign cursor world position (see code, it's already computed)
        p->vel = Vec3(0,0,0);

        // TODO: test and resolve for collisions during user movement
    }

    //* TODO: Fix relaxation steps
    for (int i = 0; i < 5; ++i) {
        for (ForceSpring* f : springsStretch) {
            auto p1 = f->getParticle1();
            auto p2 = f->getParticle2();

            auto d = (p2->pos - p1->pos);
            auto dir = d.normalized();

            auto dist = (f->getRestLength() - d.norm()) / 2;

            if (fixedParticle[p1->id]) {
                p2->pos += dir * dist * 2;
            } else if (fixedParticle[p2->id]) {
                p1->pos -= dir * dist * 2;
            } else {
                p1->pos -= dir * dist;
                p2->pos += dir * dist;
            }
        }
        /*
        for (ForceSpring* f : springsShear) {
            auto p1 = f->getParticle1();
            auto p2 = f->getParticle2();

            auto d2 = (p2->pos - p1->pos).normalized();
            auto d1 = -d2;

            auto midpoint = (p1->pos + p2->pos) * 0.5;
            auto d = f->getRestLength() * 0.5;

            if (fixedParticle[p1->id]) {
                p2->pos = p1->pos + d1 * f->getRestLength();
            } else if (fixedParticle[p2->id]) {
                p1->pos = p2->pos + d2 * f->getRestLength();
            } else {
                p1->pos = midpoint + d1 * d;
                p2->pos = midpoint + d2 * d;
            }
        }
        for (ForceSpring* f : springsBend) {
            auto p1 = f->getParticle1();
            auto p2 = f->getParticle2();

            auto d2 = (p2->pos - p1->pos).normalized();
            auto d1 = -d2;

            auto midpoint = (p1->pos + p2->pos) * 0.5;
            auto d = f->getRestLength() * 0.5;

            if (fixedParticle[p1->id]) {
                p2->pos = p1->pos + d1 * f->getRestLength();
            } else if (fixedParticle[p2->id]) {
                p1->pos = p2->pos + d2 * f->getRestLength();
            } else {
                p1->pos = midpoint + d1 * d;
                p2->pos = midpoint + d2 * d;
            }
        }
        */
    }
    //*/

    //*
    for (Particle* p : system.getParticles()) {
        p->vel = (p->pos - p->prevPos)/dt;
    }
    //*/

    double kBounce, kFriction;
    kBounce = 0.5;
    kFriction = 0.1;

    Collision colInfo;
    // collisions
    for (Particle* p : system.getParticles()) {
        // TODO: test and resolve collisions
        /*if (colliderFloor.testCollision(p, colInfo)) {
            colliderFloor.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        if (colliderRamp.testCollision(p, colInfo)) {
            colliderRamp.resolveCollision(p, colInfo, kBounce, kFriction);
        }*/
        if (colliderBall.testCollision(p, colInfo)) {
            colliderBall.resolveCollision(p, colInfo, kBounce, kFriction);
        }
        /*
        if (colliderBox.testCollision(p, colInfo)) {
            colliderBox.resolveCollision(p, colInfo, kBounce, kFriction);
        }*/
    }

    // needed after we have done collisions and relaxation, since spring forces depend on p and v
    system.updateForces();
}


void SceneCloth::mousePressed(const QMouseEvent* e, const Camera& cam)
{
    grabX = e->pos().x();
    grabY = e->pos().y();

    if (!(e->modifiers() & Qt::ControlModifier)) {

        Vec3 rayDir = cam.getRayDir(grabX, grabY);
        Vec3 origin = cam.getPos();

        selectedParticle = -1;
        for (int i = 0; i < numParticles; i++) {
            // TODO: point-ray dist to check if we select one particle
        }

        if (selectedParticle >= 0) {
            cursorWorldPos = system.getParticle(selectedParticle)->pos;
        }
    }
}

void SceneCloth::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - grabX;
    int dy = e->pos().y() - grabY;
    grabX = e->pos().x();
    grabY = e->pos().y();

    if (e->modifiers() & Qt::ControlModifier) {

    }
    else if (e->modifiers() & Qt::ShiftModifier) {

    }
    else {
        if (selectedParticle >= 0) {
            double d = -(system.getParticle(selectedParticle)->pos - cam.getPos()).dot(cam.zAxis());
            Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, d);
            cursorWorldPos += disp;
        }
    }
}

void SceneCloth::mouseReleased(const QMouseEvent*, const Camera&)
{
    selectedParticle = -1;
}

void SceneCloth::keyPressed(const QKeyEvent* e, const Camera&)
{
    if (selectedParticle >= 0 && e->key() == Qt::Key_F) {
        fixedParticle[selectedParticle] = true;
        Particle* p = system.getParticle(selectedParticle);
        p->prevPos = p->pos;
        p->vel = Vec3(0,0,0);
        p->force = Vec3(0,0,0);
    }
}
