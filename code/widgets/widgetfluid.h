#ifndef WIDGETFLUID_H
#define WIDGETFLUID_H

#include <QWidget>

namespace Ui {
class WidgetFluid;
}

class WidgetFluid : public QWidget
{
    Q_OBJECT

public:
    explicit WidgetFluid(QWidget *parent = nullptr);
    ~WidgetFluid();

    double getViscosity() const;
    double getSpeedOfSound() const;
    double getReferenceDensity() const;
    double getParticleMass() const;
    double getKernelSize() const;
    int getNeighbourhoodSize() const;
    bool getDrawColliders() const;
    double getParticleSize() const;
    int getColour() const;

    /*
    int getSolver1() const;
    int getSolver2() const;
    bool renderSameZ() const;
    bool renderTrajectory() const;

    void setSolverTypes(const std::vector<std::string>& solvers);
    void setSolver1(int idx);
    void setSolver2(int idx);
    */

signals:
    void updatedParameters();

private:
    Ui::WidgetFluid *ui;
};

#endif // WIDGETFLUID_H
