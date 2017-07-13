/*
 * CollidableSphere.h
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#pragma once

#include "geometryConcrete/collision/collidables/CollidableParticleBase.h"


namespace NuTo
{

class CollidableBase;
class CollidableWallPhysical;
class CollidableWallVirtual;
class CollidableWallCylinder;
class CollidableWallBase;

namespace Visualize
{
class UnstructuredGrid;
}

//! @brief ... class for spherical collidables
class CollidableParticleSphere : public NuTo::CollidableParticleBase
{
    //! @brief provides faster access for collision checks with walls
    friend class CollidableWallPhysical;
    friend class CollidableWallVirtual;
    friend class CollidableWallCylinder;
    friend class CollidableWallBase;

public:
    //! @brief ... constructor. Create CollidableSphere
    //! @param rPosition ... sphere position
    //! @param rVelocity ... sphere velocity
    //! @param rRadius ... sphere radius > 0
    //! @param rGrowthRate ... sphere growth rate > 0
    //! @param rIndex ... name
    CollidableParticleSphere(Eigen::Vector3d rPosition, Eigen::Vector3d rVelocity, double rRadius, double rGrowthRate,
                             int rIndex);

    //! @brief ... move spheres and apply growth, mTimeOfLastCollision update
    //! @param rTime ... new global time.
    void MoveAndGrow(double rTime) override;

    //! @brief ... calculates and returns the kinetic energy of the sphere
    double GetKineticEnergy() const override;

    //! @brief ... calculates and returns the volume of the sphere
    double GetVolume() const override;

    //! @brief ... collision between CollidableSphere and CollidableBase, resolve double dispatch, forward *this
    //! @param rCollidable ... collision partner
    void PerformCollision(CollidableBase& rCollidable) override;

    //! @brief ... collision between CollidableSphere and CollidableSphere, physics here.
    //! @param rSphere ... collision partner
    void PerformCollision(CollidableParticleSphere& rSphere) override;

    //! @brief ... collision between CollidableSphere and CollidableWall, forward to CollidableWall.
    //! @param rWall ... collision partner
    void PerformCollision(CollidableWallBase& rWall) override;

    //! @brief ... collision check between this and a CollidableBase, resolve double dispatch, forward *this
    //! @param rCollidable ... possible collision partner
    //! @param rType ... return argument, element of enum CollidableBase::EventType
    //! @return ... predicted collision time
    double PredictCollision(CollidableBase& rCollidable, int& rType) override;

    //! @brief ... collision check between this and another CollidableSphere:
    //! Physics lead to quadratic equation.
    //! @param rSphere ... possible collision partner
    //! @param rType ... return argument, element of enum CollidableBase::EventType
    //! @return ... predicted collision time
    double PredictCollision(CollidableParticleSphere& rSphere, int& rType) override;

    //! @brief ... collision check between CollidableSphere and CollidableWall, forward to CollidableWall.
    //! @param rWall ... possible collision partner
    //! @param rType ... return argument, element of enum CollidableBase::EventType
    //! @return ... predicted collision time
    double PredictCollision(CollidableWallBase& rWall, int& rType) override;

    //! @brief ... returns all old events, that need to be deleted
    //! @param rEventsToDelete ... return argument
    void GetLocalEventsToDelete(LocalEvents& rEventsToDelete) const override;

    //! @brief ... exports the sphere position and its radius to as a row in a Nx4-matrix
    //! @param rInitialRadius ... switch to export mRadius or mRadius0
    //! @return ... 1x4-matrix, [posX, posY, posZ, radius/radius0]
    Eigen::MatrixXd ExportRow(bool rInitialRadius = false) const;

#ifdef ENABLE_VISUALIZE
    //! @brief ... visualize all moving collidables
    //! @param rVisualizer ... NuTo object for ascii-export
    void VisualizationDynamic(Visualize::UnstructuredGrid& rVisualizer, bool rFinal) const;
#endif

    //! @brief ... resets the sphere velocity to 0.0
    void ResetVelocity();

    //! @brief ... sets a new growth rate
    //! @param rGrowthRateFactor ... growthrate *= rGRFactor
    //! @param rTime ... global time
    void SetGrowthRate(double rGrowthRateFactor, double rTime);

    //! @brief ... getter for sphere position
    const Eigen::Vector3d& GetPosition() const;

    //! @brief ... getter for sphere radius
    double GetRadius() const;

    //! @brief ... getter for initial sphere radius
    double GetRadius0() const;

private:
    //! @brief ... current sphere radius
    double mRadius;

    //! @brief ... base sphere radius
    double mRadiusGrowth;

    //! @brief ... initial sphere radius
    double mRadius0;

    //! @brief ... sphere growth rate
    double mGrowthRate;

    //! @brief ... initialized with 0, updated in MoveAndGrow()
    double mTimeOfLastUpdate;

    //! @brief ... initialized with 0, updated in SetGrowthRate()
    double mTimeOfGrowthReset;

    //! @brief ... physics of a one-dimensional, fully elastic sphere collision.
    //! @param rVelocity1 ... velocity of sphere 1
    //! @param rVelocity2 ... velocity of sphere 2
    //! @param rMass1 ... mass of sphere 1
    //! @param rMass2 ... mass of sphere 2
    //! @return ... post-collision velocity for sphere 1.
    double SphereCollision1D(double rVelocity1, double rVelocity2, double rMass1, double rMass2) const;

    //! @brief ... prints CollidableSphere
    //! @param rReturnStream ... output stream, that gets modified
    void Print(std::ostream& rReturnStream) const override;
};

} /* namespace NuTo */
