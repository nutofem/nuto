/*
 * CollidableWallPhysical.h
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#pragma once

#include "geometryConcrete/collision/collidables/CollidableWallBase.h"

namespace NuTo
{
class CollidableParticleSphere;

//! @brief class for planar physical walls
class CollidableWallPhysical: public NuTo::CollidableWallBase
{
public:

	//! @brief ... constructor, create CollidableWallBase using the point-and-normal-vector plane definition
	//! @param rPosition ... point on the plane
	//! @param rDirection ... normal vector pointing inside the domain, gets normalized.
	//! @param rIndex ... name
	CollidableWallPhysical(
			Eigen::VectorXd rPosition,
			Eigen::VectorXd rDirection,
			const int rIndex);

	//! @brief ... collision between CollidableWall and CollidableSphere
	//! @param rSphere ... collision partner
	void PerformCollision(CollidableParticleSphere& rSphere) override;

	//! @brief ... collision check between CollidableWall and CollidableSphere, physics here
	//! @param rSphere ... possible collision partner
	//! @return ... predicted collision time
	const double PredictCollision(CollidableParticleSphere& rSphere, int& rType) override;

	//! @brief ... returns true
	const bool IsPhysical() const override;
};

} /* namespace NuTo */
