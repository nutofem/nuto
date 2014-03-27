/*
 * CollidableWallPhysical.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"

NuTo::CollidableWallPhysical::CollidableWallPhysical(
		FullVector<double, Eigen::Dynamic> rPosition,
		FullVector<double, Eigen::Dynamic> rDirection,
		const int rIndex)
		: NuTo::CollidableWallBase(rPosition, rDirection, rIndex)
{

}

void NuTo::CollidableWallPhysical::PerformCollision(CollidableParticleSphere& rSphere)
{
	double velocityNormal = rSphere.mVelocity.dot(mDirection);

	if (velocityNormal < 0.)
		// sphere moves towards wall:
		rSphere.mVelocity += -2 * velocityNormal * mDirection;
	// else, do nothing

	rSphere.mVelocity += (rSphere.mGrowthRate) * mDirection;
}

const double NuTo::CollidableWallPhysical::PredictCollision(
		CollidableParticleSphere& rSphere, int& rType)
{
	rType = EventType::WallCollision;

	double staticDistanceToWall;
	double dynamicDistanceToWall;

	if (abs(mDirection.Sum()) == 1)
	{
		// mDirection is axis aligned

		// get axis
		int index = GetNonNullAxis();

		// get direction in this axis
		int sign = mDirection[index];

		staticDistanceToWall = sign * (rSphere.mPosition[index] - this->mPosition[index]) - rSphere.mRadius;
		dynamicDistanceToWall = sign * rSphere.mVelocity[index] - rSphere.mGrowthRate;

	}
	else
	{
		staticDistanceToWall = this->mDirection.dot(rSphere.mPosition - this->mPosition) - rSphere.mRadius;
		dynamicDistanceToWall = rSphere.mVelocity.dot(this->mDirection)	- rSphere.mGrowthRate;

	}

	if (staticDistanceToWall < 1e-15 && dynamicDistanceToWall >= 0)
		// sphere touches wall, not penetrating it in the future.
		return Event::EVENTNULL;

	double timeCollision = -staticDistanceToWall / dynamicDistanceToWall;
	if (timeCollision >= 0)
		return rSphere.mTimeOfLastUpdate + timeCollision;

	return Event::EVENTNULL;

}

