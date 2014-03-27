/*
 * CollidableWallCellGrid.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/collidables/CollidableWallVirtual.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"

NuTo::CollidableWallVirtual::CollidableWallVirtual(
		FullVector<double, 3> rPosition,
		FullVector<double, 3> rDirection,
		const int rIndex)
		: CollidableWallBase(rPosition, rDirection, rIndex)
{
}

const bool NuTo::CollidableWallVirtual::IsInOutsideBox(
		const CollidableParticleSphere& rSphere) const
		{
	const auto& boxes = rSphere.GetSubBoxes();
	for (auto box = boxes.begin(); box != boxes.end(); ++box)
		if (*box == this->mOutsideBox)
			return true;
	return false;
}

void NuTo::CollidableWallVirtual::PerformCollision(CollidableParticleSphere& rSphere)
{
	if (IsInOutsideBox(rSphere))
	{
		// sphere leaving current box
		(*mBoxes.begin())->RemoveSphere(rSphere);
	}
	else
	{
		// sphere entering new box
		mOutsideBox->AddSphere(rSphere);
	}
}

void NuTo::CollidableWallVirtual::GetDistanceAligned(double& rDynamicDistance,
		double& rStaticDistance, bool rIsInOutsideBox, CollidableParticleSphere& rSphere)
{
	int rSign = rIsInOutsideBox ? 1 : -1;
	int index = GetNonNullAxis();
	int direction = mDirection[index];
	rDynamicDistance = direction * rSphere.mVelocity[index] + rSign * rSphere.mGrowthRate;
	rStaticDistance = direction * (rSphere.mPosition[index] - this->mPosition[index]) + rSign * rSphere.mRadius;
}

void NuTo::CollidableWallVirtual::GetDistanceGeneral(double& rDynamicDistance,
		double& rStaticDistance, bool rIsInOutsideBox,
		CollidableParticleSphere& rSphere)
{
	int rSign = rIsInOutsideBox ? 1 : -1;
	rDynamicDistance = rSphere.mVelocity.dot(this->mDirection) + rSign * rSphere.mGrowthRate;
	rStaticDistance = this->mDirection.dot(rSphere.mPosition - this->mPosition)	+ rSign * rSphere.mRadius;
}

const double NuTo::CollidableWallVirtual::PredictCollision(
		CollidableParticleSphere& rSphere, int& rType)
{
	rType = EventType::WallTransfer;

	bool isInOutsideBox = IsInOutsideBox(rSphere);
	if (isInOutsideBox)
	{

		double dynamicDistanceToExit;
		double staticDistanceToExit;

		if (abs(mDirection.sum()) == 1)
			GetDistanceAligned(dynamicDistanceToExit, staticDistanceToExit,	isInOutsideBox, rSphere);
		else
			GetDistanceGeneral(dynamicDistanceToExit, staticDistanceToExit,	isInOutsideBox, rSphere);

		if (dynamicDistanceToExit > 0.)
			return Event::EVENTNULL;

		double timeCollisionExit = -staticDistanceToExit / dynamicDistanceToExit;
		if (timeCollisionExit >= 0.)
			return rSphere.mTimeOfLastUpdate + timeCollisionExit;
	}
	else
	{
		double dynamicDistanceToWall;
		double staticDistanceToWall;

		if (abs(mDirection.sum()) == 1)
			GetDistanceAligned(dynamicDistanceToWall, staticDistanceToWall,	isInOutsideBox, rSphere);
		else
			GetDistanceGeneral(dynamicDistanceToWall, staticDistanceToWall,	isInOutsideBox, rSphere);

		if (dynamicDistanceToWall > 0.)
			return Event::EVENTNULL;

		if (staticDistanceToWall < 0. && staticDistanceToWall > -1.e-12 * rSphere.mRadius)
			staticDistanceToWall = 0.;

		double timeCollisionWall = -staticDistanceToWall / dynamicDistanceToWall;
		if (timeCollisionWall >= 0.)
			return rSphere.mTimeOfLastUpdate + timeCollisionWall;

	}
	return Event::EVENTNULL;

}

