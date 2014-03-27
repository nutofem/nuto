/*
 * CollidableParticleBase.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"

NuTo::CollidableParticleBase::CollidableParticleBase(
		FullVector<double, Eigen::Dynamic> rPosition,
		FullVector<double, Eigen::Dynamic> rVelocity,
		const int rIndex)
		: CollidableBase(rIndex)
{

	if (rPosition.GetNumRows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] position vector must have exactly 3 rows");

	if (rVelocity.GetNumRows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] velocity vector must have exactly 3 rows");

	mPosition = rPosition;
	mVelocity = rVelocity;
}

void NuTo::CollidableParticleBase::Print(std::ostream& rReturnStream) const
{
}
