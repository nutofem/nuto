/*
 * CollidableWallCellGrid.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallVirtual.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/SubBox.h"

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
	for(unsigned int iBox = 0; iBox < boxes.size(); ++iBox)
		if(boxes[iBox] == this->mOutsideBox)
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
	int direction = mDirection[mNonNullAxis];
	rDynamicDistance = direction * rSphere.mVelocity[mNonNullAxis] + rSign * rSphere.mGrowthRate;
	rStaticDistance = direction * (rSphere.mPosition[mNonNullAxis] - this->mPosition[mNonNullAxis]) + rSign * rSphere.mRadius;
}

const bool NuTo::CollidableWallVirtual::IsPhysical() const
{
	return false;
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
	rType = Event::EventType::WallTransfer;

	bool isInOutsideBox = IsInOutsideBox(rSphere);
	if (isInOutsideBox)
	{

		double dynamicDistanceToExit;
		double staticDistanceToExit;

		if (mIsAxisAligned)
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

		if (mIsAxisAligned)
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

