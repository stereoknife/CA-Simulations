#include "widgetfluid.h"
#include "ui_widgetfluid.h"

WidgetFluid::WidgetFluid(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WidgetFluid)
{
    ui->setupUi(this);
    connect(ui->update, &QPushButton::clicked, this,
            [=] (void) { emit updatedParameters(); });
}

WidgetFluid::~WidgetFluid()
{
    delete ui;
}

double WidgetFluid::getParticleMass() const {
    return ui->mass->value();
}

double WidgetFluid::getReferenceDensity() const {
    return ui->density->value();
}

double WidgetFluid::getSpeedOfSound() const {
    return ui->sound_speed->value();
}

double WidgetFluid::getViscosity() const {
    return ui->viscosity->value();
}

double WidgetFluid::getKernelSize() const {
    return ui->kernel_size->value();
}

double WidgetFluid::getSurfaceTension() const {
    return ui->surface_tension->value();
}

bool WidgetFluid::getDrawColliders() const {
    return ui->draw_colliders->isChecked();
}

double WidgetFluid::getParticleSize() const {
    return ui->particle_size->value();
}

int WidgetFluid::getColour() const {
    return ui->colour->currentIndex();
}

double WidgetFluid::getGravity() const {
    return ui->gravity->value();
}

/*
bool WidgetFluid::renderTrajectory() const {
    return ui->trajectory->isChecked();
}

void WidgetFluid::setSolverTypes(const std::vector<std::string>& solvers) {
    ui->solver1->clear();
    ui->solver2->clear();
    for (const std::string& s : solvers) {
        ui->solver1->addItem(QString::fromStdString(s));
        ui->solver2->addItem(QString::fromStdString(s));
    }
}

void WidgetFluid::setSolver1(int idx) {
    ui->solver1->setCurrentIndex(idx);
}

void WidgetFluid::setSolver2(int idx) {
    ui->solver2->setCurrentIndex(idx);
}
*/
