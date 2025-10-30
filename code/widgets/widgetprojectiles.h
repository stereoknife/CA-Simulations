#ifndef WIDGETPROJECTILES_H
#define WIDGETPROJECTILES_H

#include <QWidget>

namespace Ui {
class WidgetFluid;
}

class WidgetProjectiles : public QWidget
{
    Q_OBJECT

public:
    explicit WidgetProjectiles(QWidget *parent = nullptr);
    ~WidgetProjectiles();

    double getGravity() const;
    double getHeight() const;
    double getAngle() const;
    double getSpeed() const;
    int getSolver1() const;
    int getSolver2() const;
    bool renderSameZ() const;
    bool renderTrajectory() const;

    void setSolverTypes(const std::vector<std::string>& solvers);
    void setSolver1(int idx);
    void setSolver2(int idx);

private:
    Ui::WidgetFluid *ui;
};

#endif // WIDGETPROJECTILES_H
