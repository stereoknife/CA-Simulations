#include "colliders.h"
#include <cmath>
#include <iostream>

// Helper function because it will need to be reused a bunch
// Returns the distance along a line at which the collision happens
double rayPlaneCollision(const Vec3 origin, const Vec3 dir, const Vec3 plane_n, const double plane_d){
    auto nr = plane_n.dot(dir);
    return -(plane_n.dot(origin)+plane_d)/nr;
}

double max(double a, double b, double c){
    return std::max(a, std::max(b, c));
}

double min(double a, double b, double c){
    return std::min(a, std::min(b, c));
}

double sign(double x) {
    return x > 0 ? 1 : x < 0 ? -1 : 0;
}

/*
 * Generic function for collision response from contact plane
 */
void Collider::resolveCollision(Particle* p, const Collision& col, double kElastic, double kFriction) const
{
    auto position = p->pos;
    auto collision = col.position;
    auto h = col.normal.dot(position-collision);
    auto corrected_position = position - (1 + kElastic) * h * col.normal;

    p->pos = corrected_position;

    auto vel_normal = col.normal.dot(p->vel) * col.normal;
    auto vel_tangent = p->vel - vel_normal;

    p->vel = vel_normal * -kElastic + (1 - kFriction) * vel_tangent;
}



/*
 * Plane
 */
bool ColliderPlane::isInside(const Particle* p) const
{
    auto o = p->prevPos;
    auto r = p->pos - o;
    auto nr = planeN.dot(r);

    // If the dot product between the normal and the ray is 0 they are parallel and there is no collision
    // Unless the ray is on the plane, in which case the collision is at 0
    if (nr == 0) return false;

    // Get the distance along the ray at which the collision happens and only return true if it's positive, otherwise
    // the collision happens behind the origin of the ray.
    auto l = -(planeN.dot(o)+planeD)/nr;

    return 0 <= l;
}


bool ColliderPlane::testCollision(const Particle* p, Collision& colInfo) const
{
    auto o = p->prevPos;
    auto r = p->pos - o;

    // get lambda
    auto l = rayPlaneCollision(o, r, planeN, planeD);

    // particle trajectory is parallel to plane
    if (std::isinf(l)) return false;

    colInfo.position = o + r * l;
    colInfo.normal = planeN;

    return 0 <= l && l <= 1;
}



/*
 * Sphere
 */
bool ColliderSphere::isInside(const Particle* p) const
{
    // Dist between particle and sphere
    auto d = (p->pos - center).squaredNorm();
    // Check if smaller than radius
    return d < radius * radius;
}


bool ColliderSphere::testCollision(const Particle* p, Collision& colInfo) const
{
    // Collision detection method from 3D Math Primer for Graphics and Game Development 2nd Edition by Fletcher Dunn and Ian Parberry
    auto p0 = p->prevPos;
    auto c = getCenter();
    auto rr = getRadius() * getRadius();
    Vec3 v = p->pos - p->prevPos;
    Vec3 d = v.normalized();
    Vec3 e = c - p0;
    auto ee = e.dot(e);
    if (ee < rr) return false; // Origin inside the sphere
    auto a = e.dot(d);
    auto ff = rr - ee + a * a;
    if (ff < 0) return false;
    auto l = a-sqrt(ff);
    auto pi = p0 + d * l;
    colInfo.position = pi;
    colInfo.normal = (pi - c).normalized();
    return l * l <= v.squaredNorm();
}



/*
 * AABB
 */
bool ColliderAABB::isInside(const Particle* p) const
{
    auto pos = p->pos;
    bool inside = true;
    inside &= bmin.x() < pos.x() && pos.x() < bmax.x();
    inside &= bmin.y() < pos.y() && pos.y() < bmax.y();
    inside &= bmin.z() < pos.z() && pos.z() < bmax.z();
    return inside;
}


bool ColliderAABB::testCollision(const Particle* p, Collision& colInfo) const
{
    const Vec3  nx = Vec3{ 1.0,  0.0,  0.0};
    const Vec3  ny = Vec3{ 0.0,  1.0,  0.0};
    const Vec3  nz = Vec3{ 0.0,  0.0,  1.0};

    auto o = p->prevPos;
    auto d = p->pos - o;

    Vec3 omin, omax;

    if (d.x() > 0) {omin.x() = bmin.x(); omax.x() = bmax.x();} else {omin.x() = bmax.x(); omax.x() = bmin.x();}
    if (d.y() > 0) {omin.y() = bmin.y(); omax.y() = bmax.y();} else {omin.y() = bmax.y(); omax.y() = bmin.y();}
    if (d.z() > 0) {omin.z() = bmin.z(); omax.z() = bmax.z();} else {omin.z() = bmax.z(); omax.z() = bmin.z();}

    // X dimension collisions
    auto mxl = rayPlaneCollision(o, d, -nx, omin.x());
    auto xxl = rayPlaneCollision(o, d, -nx, omax.x());

    auto myl = rayPlaneCollision(o, d, -ny, omin.y());
    auto xyl = rayPlaneCollision(o, d, -ny, omax.y());

    auto mzl = rayPlaneCollision(o, d, -nz, omin.z());
    auto xzl = rayPlaneCollision(o, d, -nz, omax.z());

    Vec3 n;
    double m;
    double x = min(xxl, xyl, xzl);

    if (mxl > myl && mxl > mzl) {
        m = mxl;
        colInfo.normal = nx * sign(d.x());
    } else if (myl > mzl) {
        m = myl;
        colInfo.normal = ny * sign(d.y());
    } else {
        m = mzl;
        colInfo.normal = nz * sign(d.z());
    }

    colInfo.position = o + d * m;

    return m < x && m >= 0 && m <= 1;
}
