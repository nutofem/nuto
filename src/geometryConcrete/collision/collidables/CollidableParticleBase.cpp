/*
 * CollidableParticleBase.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "math/FullMatrix.h"
#include "geometryConcrete/collision/collidables/CollidableParticleBase.h"

NuTo::CollidableParticleBase::CollidableParticleBase(
		Eigen::VectorXd rPosition,
		Eigen::VectorXd rVelocity,
		const int rIndex)
		: CollidableBase(rIndex), mPosition(rPosition), mVelocity(rVelocity)
{

	if (mPosition.rows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] position vector must have exactly 3 rows");

	if (mVelocity.rows() != 3)
		throw NuTo::Exception("[NuTo::CollidableParticleBase::CollidableParticleBase] velocity vector must have exactly 3 rows");
}
