/*
 * CollidableWallCellGrid.h
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#pragma once

#include "geometryConcrete/collision/collidables/CollidableWallBase.h"

namespace NuTo
{
class CollidableParticleSphere;

//! @brief ... virtual, inner, planar wall of the specimen
class CollidableWallVirtual : public NuTo::CollidableWallBase
{
public:
    //! @brief ... constructor, create CollidableWallBase using the point-and-normal-vector plane definition
    //! @param rPosition ... point on the plane
    //! @param rDirection ... normal vector pointing inside the domain, gets normalized.
    //! @param rIndex ... name
    CollidableWallVirtual(Eigen::Vector3d rPosition, Eigen::Vector3d rDirection, int rIndex);

    //! @brief ... collision between CollidableWall and CollidableSphere
    //! Case 1) Sphere is NOT in outside box of this wall
    //!     --> Sphere has just passed this wall and needs to be added in the outside box
    //! Case 2) Sphere IS in the outside box
    //!     --> Sphere has just left the inside box and needs to be removed from there
    //! @param rSphere ... collision partner
    void PerformCollision(CollidableParticleSphere& rSphere) override;

    //! @brief ... collision check between CollidableWall and CollidableSphere
    //! Predict wall transfer events
    //! @param rSphere ... possible collision partner
    //! @return ... predicted collision time
    double PredictCollision(CollidableParticleSphere& rSphere, int& rType) override;

    //! @brief ... returns false
    bool IsPhysical() const override
    {
        return false;
    }

private:
    //! @brief ... returns whether a sphere is in the outside box of this wall
    //! @param rSphere ... sphere to test
    bool IsInOutsideBox(const CollidableParticleSphere& rSphere) const;

    //! @brief ... get the distance between the sphere surface and the wall
    void GetDistanceAligned(double& rDynamicDistance, double& rStaticDistance, bool rIsInOutsideBox,
                            CollidableParticleSphere& rSphere);

    //! @brief ... get the distance between the sphere surface and the wall
    void GetDistanceGeneral(double& rDynamicDistance, double& rStaticDistance, bool rIsInOutsideBox,
                            CollidableParticleSphere& rSphere);
};

} /* namespace NuTo */
