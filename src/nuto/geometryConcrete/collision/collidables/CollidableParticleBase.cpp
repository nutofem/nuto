/*
 * CollidableParticleBase.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/math/FullMatrix.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"

NuTo::CollidableParticleBase::CollidableParticleBase(
		FullVector<double, Eigen::Dynamic> rPosition,
		FullVector<double, Eigen::Dynamic> rVelocity,
		const int rIndex)
		: CollidableBase(rIndex), mPosition(rPosition), mVelocity(rVelocity)
{

	if (mPosition.GetNumRows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] position vector must have exactly 3 rows");

	if (mVelocity.GetNumRows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] velocity vector must have exactly 3 rows");
}
