/*
 * CollidableWall.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#include "math/FullMatrix.h"
#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "geometryConcrete/collision/SubBox.h"
#include "visualize/VisualizeUnstructuredGrid.h"

NuTo::CollidableWallBase::CollidableWallBase(
		Eigen::VectorXd rPosition,
		Eigen::VectorXd rDirection,
		const int rIndex)
		: NuTo::CollidableBase::CollidableBase(rIndex),
		  mPosition(rPosition),
		  mDirection(rDirection),
		  mInsideBox(nullptr),
		  mOutsideBox(nullptr),
		  mNonNullAxis(GetNonNullAxis()),
		  mIsAxisAligned(std::abs(mDirection.Sum()) == 1)
{
	mDirection.normalize();
}

NuTo::CollidableWallBase::~CollidableWallBase()
{
}

void NuTo::CollidableWallBase::PerformCollision(CollidableBase& rCollidable)
{
	rCollidable.PerformCollision(*this);
}
void NuTo::CollidableWallBase::PerformCollision(CollidableWallBase& rWall)
{
// nothing here.
}

const double NuTo::CollidableWallBase::PredictCollision(
		CollidableBase& rCollidable, int& rType)
{
	return rCollidable.PredictCollision(*this, rType);
}

const double NuTo::CollidableWallBase::PredictCollision(
		CollidableWallBase& rWall, int& rType)
{
	return -1;
}

void NuTo::CollidableWallBase::Print(std::ostream& rReturnStream) const
{
	int indexWidth = 3;
	int precision = 3;
	rReturnStream.precision(precision);
	rReturnStream << "Wall   "
			<< std::setw(indexWidth) << mIndex << " Pos:("
			<< std::setw(precision + 2) << mPosition[0] << ","
			<< std::setw(precision + 2) << mPosition[1] << ","
			<< std::setw(precision + 2) << mPosition[2] << ") Dir:("
			<< std::setw(precision + 2) << mDirection[0] << ","
			<< std::setw(precision + 2) << mDirection[1] << ","
			<< std::setw(precision + 2) << mDirection[2] << ")";
}

#ifdef ENABLE_VISUALIZE

void NuTo::CollidableWallBase::VisualizationStatic(
		NuTo::VisualizeUnstructuredGrid& rVisualizer) const
{

	double size = 1.;
	if (*mBoxes.begin() == mOutsideBox)
		size = 2.;

	NuTo::FullMatrix<double, 4, 3> corners;

	// get some vector != mDirection
	NuTo::FullVector<double, 3> random;
	random << 1, 0, 0;
	if (std::abs(random.Dot(mDirection)) == 1)
	{
		random << 0, 1, 0;
	}

	NuTo::FullVector<double, 3> transversal = random.cross(mDirection);
	NuTo::FullVector<double, 3> transversal2 = transversal.cross(mDirection);

	// normalize to size/2;
	transversal.normalize();
	transversal2.normalize();

	transversal *= size / 2;
	transversal2 *= size / 2;

	corners.SetRow(0, (mPosition + transversal + transversal2).transpose());
	corners.SetRow(1, (mPosition + transversal - transversal2).transpose());
	corners.SetRow(2, (mPosition - transversal - transversal2).transpose());
	corners.SetRow(3, (mPosition - transversal + transversal2).transpose());


	unsigned int cornerIndex[4];
	for (int i = 0; i < 4; ++i)
	{
		cornerIndex[i] = rVisualizer.AddPoint(corners.GetRow(i).data());
	}
	unsigned int insertIndex = rVisualizer.AddQuadCell(cornerIndex);

	double tmpDirection[3];
	tmpDirection[0] = mDirection[0];
	tmpDirection[1] = mDirection[1];
	tmpDirection[2] = mDirection[2];
	rVisualizer.SetCellDataVector(insertIndex, "Direction", tmpDirection);
}
#endif

void NuTo::CollidableWallBase::SetBoxes(SubBox& rInsideBox, SubBox& rOutsideBox)
{
	mBoxes.push_back(&rInsideBox);
	mInsideBox = &rInsideBox;
	mOutsideBox = &rOutsideBox;
}

void NuTo::CollidableWallBase::GetLocalEventsToDelete(Event::LocalEvents& rEventsToDelete) const
{
}

void NuTo::CollidableWallBase::MoveAndGrow(const double rTime)
{
}

bool NuTo::CollidableWallBase::IsInside(const CollidableParticleSphere& rSphere) const
{
	const FullVector<double, 3>& dPosition = rSphere.mPosition - mPosition;

	double distanceToWall = mDirection.dot(dPosition);

	return (distanceToWall + rSphere.mRadius >= 0);

}

const Eigen::VectorXd NuTo::CollidableWallBase::GetDirection() const
{
	return mDirection;
}

const Eigen::VectorXd NuTo::CollidableWallBase::GetPosition() const
{
	return mPosition;
}

const int NuTo::CollidableWallBase::GetNonNullAxis()
{
	for (int i = 0; i < 3; i++)
		if (mDirection[i] != 0.)
			return i;
	return -1;
}
