/*
 * CollidableSphere.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#include <iomanip>

#include "nuto/base/Exception.h"
#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallVirtual.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallCylinder.h"

#include "nuto/visualize/UnstructuredGrid.h"

NuTo::CollidableParticleSphere::CollidableParticleSphere(Eigen::Vector3d rPosition, Eigen::Vector3d rVelocity,
                                                         double rRadius, double rGrowthRate, int rIndex)
    : CollidableParticleBase(rPosition, rVelocity, rIndex)
    , mRadius(rRadius)
    , mRadiusGrowth(rRadius)
    , mRadius0(rRadius)
    , mGrowthRate(rGrowthRate)
    , mTimeOfLastUpdate(0.)
    , mTimeOfGrowthReset(0.)
{
}

void NuTo::CollidableParticleSphere::MoveAndGrow(double rTime)
{
    double dTUpdate = rTime - mTimeOfLastUpdate;
    mPosition += mVelocity * dTUpdate;

    double dTGrowth = rTime - mTimeOfGrowthReset;
    mRadius = mRadiusGrowth + mGrowthRate * dTGrowth;

    mTimeOfLastUpdate = rTime;
}

void NuTo::CollidableParticleSphere::SetGrowthRate(double rGrowthRateFactor, double rTime)
{
    // calculate new growth radius
    double dTGrowth = rTime - mTimeOfGrowthReset;
    mRadiusGrowth += mGrowthRate * dTGrowth;

    // update fields
    mGrowthRate *= rGrowthRateFactor;
    mTimeOfGrowthReset = rTime;
}

double NuTo::CollidableParticleSphere::SphereCollision1D(double rVelocity1, double rVelocity2, double rMass1,
                                                         double rMass2) const
{
    return (rVelocity1 * (rMass1 - rMass2) + 2 * rVelocity2 * rMass2) / (rMass1 + rMass2);
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableBase& rCollidable)
{
    rCollidable.PerformCollision(*this);
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableParticleSphere& rSphere)
{
    const Eigen::Vector3d& velocity1 = this->mVelocity;
    const Eigen::Vector3d& velocity2 = rSphere.mVelocity;

    double mass1 = this->mRadius * this->mRadius * this->mRadius;
    double mass2 = rSphere.mRadius * rSphere.mRadius * rSphere.mRadius;

    // 1) Get normal direction
    Eigen::Vector3d n = this->mPosition - rSphere.mPosition;
    n.normalize();

    // 2) Velocity split:
    // 2.1) normal velocity
    double velocityNormal1 = n.dot(velocity1);
    double velocityNormal2 = n.dot(velocity2);

    // 2.2) tangential velocity
    Eigen::Vector3d velocityTransversal1 = velocity1 - n * velocityNormal1;
    Eigen::Vector3d velocityTransversal2 = velocity2 - n * velocityNormal2;

    // 3.1) 1D collision
    double velocityNormalNew1 = SphereCollision1D(velocityNormal1, velocityNormal2, mass1, mass2);
    double velocityNormalNew2 = SphereCollision1D(velocityNormal2, velocityNormal1, mass2, mass1);

    // 3.2) make sure, that the new normal velocities are at least
    // fast as the growth rates --> no simultaneous re-collision

    double velocityNormalEffective = velocityNormalNew1 - velocityNormalNew2;

    if (velocityNormalEffective < 0.)
    {
        // add the missing effective velocity to the spheres
        double velocityDifference1 = -velocityNormalEffective * mass2 / (mass1 + mass2);
        double velocityDifference2 = -velocityNormalEffective - velocityDifference1;

        velocityNormalNew1 += velocityDifference1;
        velocityNormalNew2 -= velocityDifference2;
    }

    // 3.3) add normal velocity boost
    double velocityExtra = 2.0;
    velocityNormalNew1 += velocityExtra * this->mGrowthRate;
    velocityNormalNew2 -= velocityExtra * rSphere.mGrowthRate;

    // 4) combine normal and tangential velocities, add growth rate for collision buffer
    this->mVelocity = (velocityNormalNew1)*n + velocityTransversal1;
    rSphere.mVelocity = (velocityNormalNew2)*n + velocityTransversal2;
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableWallBase& rWall)
{
    rWall.PerformCollision(*this);
}

double NuTo::CollidableParticleSphere::PredictCollision(CollidableBase& rCollidable, int& rType)
{
    return rCollidable.PredictCollision(*this, rType);
}

double NuTo::CollidableParticleSphere::PredictCollision(CollidableParticleSphere& rSphere, int& rType)
{
    if (this == &rSphere)
        return Event::EVENTNULL;

    rType = Event::EventType::SphereCollision;

    // sync both spheres to the more recent time
    double baseTime;

    // RAII!
    CollidableParticleSphere& s1 = this->mTimeOfLastUpdate > rSphere.mTimeOfLastUpdate ? *this : rSphere;
    CollidableParticleSphere& s2 = this->mTimeOfLastUpdate > rSphere.mTimeOfLastUpdate ? rSphere : *this;
    baseTime = s1.mTimeOfLastUpdate;

    Eigen::Vector3d dP = s1.mPosition - (s2.mPosition + s2.mVelocity * (s1.mTimeOfLastUpdate - s2.mTimeOfLastUpdate));
    Eigen::Vector3d dV = s1.mVelocity - s2.mVelocity;

    double dR = s1.mRadius + s2.mRadius + s2.mGrowthRate * (s1.mTimeOfLastUpdate - s2.mTimeOfLastUpdate);
    double dG = s1.mGrowthRate + s2.mGrowthRate;

    double a = dV.dot(dV) - dG * dG;
    double b = 2 * (dV.dot(dP) - dR * dG);
    double c = dP.dot(dP) - dR * dR;

    if (c < -2.e-10 * s1.mRadius)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Sphere " + std::to_string(s1.mIndex) + " overlaps with Sphere " +
                                                           std::to_string(s2.mIndex));
    }

    double discriminant = b * b - 4 * a * c;

    // allow small negative values ...
    if (discriminant < -1e-12)
        return Event::EVENTNULL;

    // ... an treat them as zero
    discriminant = std::max(discriminant, 0.);

    if (b < 0.)
        return baseTime + 2 * c / (-b + std::sqrt(discriminant));

    if (a < 0.)
        return baseTime + (-b - std::sqrt(discriminant)) / (2 * a);

    return Event::EVENTNULL;
}

double NuTo::CollidableParticleSphere::PredictCollision(CollidableWallBase& rWall, int& rType)
{
    return rWall.PredictCollision(*this, rType);
}

Eigen::MatrixXd NuTo::CollidableParticleSphere::ExportRow(bool rInitialRadius) const
{
    Eigen::MatrixXd data(1, 4);
    double radius = rInitialRadius ? mRadius0 : mRadius;
    data << mPosition[0], mPosition[1], mPosition[2], radius;
    return data;
}

const Eigen::Vector3d& NuTo::CollidableParticleSphere::GetPosition() const
{
    return mPosition;
}

double NuTo::CollidableParticleSphere::GetRadius() const
{
    return mRadius;
}

double NuTo::CollidableParticleSphere::GetRadius0() const
{
    return mRadius0;
}

void NuTo::CollidableParticleSphere::Print(std::ostream& rReturnStream) const
{
    int indexWidth = 3;
    int precision = 20;
    rReturnStream.setf(std::ios::scientific);
    rReturnStream.precision(precision);
    rReturnStream << "Sphere " << std::setw(indexWidth) << mIndex << " Pos:(" << std::setw(precision + 2)
                  << mPosition[0] << "," << std::setw(precision + 2) << mPosition[1] << "," << std::setw(precision + 2)
                  << mPosition[2] << ") Vel:(" << std::setw(precision + 2) << mVelocity[0] << ","
                  << std::setw(precision + 2) << mVelocity[1] << "," << std::setw(precision + 2) << mVelocity[2]
                  << ") R:" << std::setw(precision + 2) << mRadius;
}

void NuTo::CollidableParticleSphere::GetLocalEventsToDelete(Event::LocalEvents& rEventsToDelete) const
{
    for (auto mLocalEvent : mLocalEvents)
    {
        rEventsToDelete.push_back(mLocalEvent);
    }
}


void NuTo::CollidableParticleSphere::VisualizationDynamic(NuTo::Visualize::UnstructuredGrid& rVisualizer,
                                                          bool rFinal) const
{
    unsigned int index = rVisualizer.AddPoint(mPosition);
    double radius = rFinal ? mRadius0 : mRadius;
    rVisualizer.SetPointData(index, "Radius", radius);
    rVisualizer.SetPointData(index, "Velocity", mVelocity);
}


double NuTo::CollidableParticleSphere::GetKineticEnergy() const
{
    return .5 * mRadius * mRadius * mRadius * mVelocity.dot(mVelocity);
}

double NuTo::CollidableParticleSphere::GetVolume() const
{
    return 4. / 3. * M_PI * mRadius * mRadius * mRadius;
}

void NuTo::CollidableParticleSphere::ResetVelocity()
{
    mVelocity << 0, 0, 0;
}
