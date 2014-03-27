/*
 * CollidableWallPhysical.h
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#ifndef COLLIDABLEWALLPHYSICAL_H_
#define COLLIDABLEWALLPHYSICAL_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"

namespace NuTo
{
class CollidableParticleSphere;
class CollidableWallPhysical: public NuTo::CollidableWallBase
{
public:
	CollidableWallPhysical(
			NuTo::FullVector<double, Eigen::Dynamic> rPosition,
			NuTo::FullVector<double, Eigen::Dynamic> rDirection,
			const int rIndex);

	//! @brief ... collision between CollidableWall and CollidableSphere
	//! @param rSphere ... collision partner
	void PerformCollision(CollidableParticleSphere& rSphere);

	//! @brief ... collision check between CollidableWall and CollidableSphere, physics here
	//! @param rSphere ... possible collision partner
	//! @return ... predicted collision time
	const double PredictCollision(CollidableParticleSphere& rSphere, int& rType);

};

} /* namespace NuTo */
#endif /* COLLIDABLEWALLPHYSICAL_H_ */
