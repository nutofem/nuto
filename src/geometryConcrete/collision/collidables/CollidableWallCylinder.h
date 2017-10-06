/*
 * CollidableWallCylinder.h
 *
 *  Created on: 13 Feb 2014
 *      Author: ttitsche
 */

#pragma once

#include "geometryConcrete/collision/collidables/CollidableWallBase.h"

namespace NuTo
{
class CollidableParticleSphere;

//! @brief class for cylindric walls
class CollidableWallCylinder : public NuTo::CollidableWallBase
{
public:
    //! @brief ... constructor.
    //! @param rPosition ... center of the cylinder
    //! @param rDirection ... vector pointing to the top of the cylinder
    //! @param rRadius ... radius
    //! @param rHeigth ... height
    //! @param rIndex ... name
    CollidableWallCylinder(Eigen::Vector3d rPosition, Eigen::Vector3d rDirection, double rRadius, double rHeigth,
                           int rIndex);

    //! @brief ... collision between this and CollidableSphere
    //! @param rSphere ... collision partner
    void PerformCollision(CollidableParticleSphere& rSphere) override;

    //! @brief ... collision prediction between (this) and a sphere
    //! @param rSphere ... collision partner
    //! @return ... predicted time of collision
    double PredictCollision(CollidableParticleSphere& rSphere, int& rType) override;

    bool IsPhysical() const override
    {
        return true;
    }

    //! @brief ... visualize all non-moving collidables
    //! @param rVisualizer ... NuTo object for ascii-export
    void VisualizationStatic(Visualize::UnstructuredGrid& rVisualizer) const override;

    //! @brief ... returns true, if the sphere is in this cylinder
    //! @param rSphere ... sphere to check
    bool IsInside(const CollidableParticleSphere& rSphere) const override;

protected:
    //! @brief ... prints cylinder informations
    //! @param rReturnStream ... output stream, that gets modified
    void Print(std::ostream& rReturnStream) const override;

    //! @brief ... cylinder radius
    double mRadius;

    //! @brief ... cylinder height
    double mHeigth;
};

} /* namespace NuTo */
