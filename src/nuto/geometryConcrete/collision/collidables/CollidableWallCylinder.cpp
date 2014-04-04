/*
 * CollidableWallCylinder.cpp
 *
 *  Created on: 13 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/collidables/CollidableWallCylinder.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"

NuTo::CollidableWallCylinder::CollidableWallCylinder(
		FullVector<double, 3> rPosition,
		FullVector<double, 3> rDirection,
		const double rRadius,
		const double rHeigth,
		const int rIndex)
		: CollidableWallBase(rPosition, rDirection, rIndex),
				mRadius(rRadius),
				mHeigth(rHeigth)
{
}

void NuTo::CollidableWallCylinder::PerformCollision(CollidableParticleSphere& rSphere)
{
	// get perpendicular from sphere centre to direction vector
	FullVector<double, 3> dP = mPosition - rSphere.mPosition;
	FullVector<double, 3> n = dP - dP.Dot(mDirection) * mDirection;
	if (n.Dot(n) == 0)
		throw NuTo::Exception("[NuTo::CollidableWallCylinder::PerformCollision] n = 0");

	n.normalize();

	// perform plane wall collision:
	double velocityNormal = rSphere.mVelocity.dot(n);

	if (velocityNormal < 0.)
		rSphere.mVelocity += -2 * velocityNormal * n;

	rSphere.mVelocity += (rSphere.mGrowthRate) * n;
}

const double NuTo::CollidableWallCylinder::PredictCollision(
		CollidableParticleSphere& rSphere, int& rType)
{
	rType = EventType::WallCollision;

	FullVector<double, 3> dPos = rSphere.mPosition - mPosition;

	FullVector<double, 3> dP = dPos - dPos.Dot(mDirection) * mDirection;
	FullVector<double, 3> dV = rSphere.mVelocity - rSphere.mVelocity.Dot(mDirection) * mDirection;
	long double dR = mRadius - rSphere.mRadius;
	long double dG = -rSphere.mGrowthRate;

	long double a = dV.dot(dV) - dG * dG;
	long double c = dP.dot(dP) - dR * dR;
	long double b = 2 * (dV.dot(dP) - dR * dG);

	long double timeCollision = 0;

	if (c > 2.0e-10 * rSphere.mRadius)
		throw NuTo::Exception("[NuTo::CollidableSphere::CreateNewEvent] Sphere overlapping with Cylinder!");

	long double discriminant = b * b - 4 * a * c;
	if (discriminant > -1e-11)
	{
		if (discriminant < 0.)
			discriminant = 0;

		if (b > 0.)
		{

			timeCollision = 2 * c / (-b - std::sqrt(discriminant));
			if (timeCollision >= 0.)
				return rSphere.mTimeOfLastUpdate + timeCollision;

		}
		else if (a > 0. && b <= 0.)
		{
			timeCollision = (-b + std::sqrt(discriminant)) / (2 * a);
			return rSphere.mTimeOfLastUpdate + timeCollision;
		}
	}
	return Event::EVENTNULL;
}

#ifdef ENABLE_VISUALIZE
void NuTo::CollidableWallCylinder::VisualizationStatic(
		VisualizeUnstructuredGrid& rVisualizer) const
		{
	// =========================================
	// create a circle out of several line cells
	// =========================================
	const int nPoints = 40; // resolution of the circle;

	FullVector<double, 3> n = mDirection; // rename
	FullVector<double, 3> nonCoLin( { 1, 0, 0 });
	if (nonCoLin == n)
		nonCoLin << 0, 1, 0;

	FullVector<double, 3> sVec = n.cross(nonCoLin); // random vector perpendicular to mDirection = n
	sVec.normalize();

	FullMatrix<double, nPoints + 1, 3> pointCircle;	// store circle coordinates
	double dAlpha = 2. * M_PI / nPoints;
	for (int i = 0; i < nPoints + 1; i++)
	{
		// calculate circle points and add them to points
		// point at nPoints+1 closes the circle
		double alpha = i * dAlpha;
		FullVector<double, 3> rotVec = n * (n.dot(sVec))
				+ std::cos(alpha) * n.cross(sVec).cross(n)
				+ std::sin(alpha) * n.cross(sVec);
		pointCircle.SetRow(i, (mRadius * rotVec).transpose());
	}

	// =========================================
	//   create and add top and bottom points
	// =========================================
	FullVector<double, 3> pointTop = mPosition + n * mHeigth / 2.;
	FullVector<double, 3> pointBot = mPosition - n * mHeigth / 2.;
	unsigned int indexTop = rVisualizer.AddPoint(pointTop.data());
	unsigned int indexBot = rVisualizer.AddPoint(pointBot.data());

	unsigned int pointCircleIndexTop[nPoints + 1];
	unsigned int pointCircleIndexBot[nPoints + 1];
	for (int i = 0; i < nPoints + 1; i++)
	{
		FullVector<double, 3> pointCircleTop = pointCircle.GetRow(i).transpose() + pointTop;
		FullVector<double, 3> pointCircleBot = pointCircle.GetRow(i).transpose() + pointBot;
		pointCircleIndexTop[i] = rVisualizer.AddPoint(pointCircleTop.data());
		pointCircleIndexBot[i] = rVisualizer.AddPoint(pointCircleBot.data());
	}

	// =========================================
	//        create top and bottom planes
	// =========================================
	for (int i = 0; i < nPoints; i++)
	{
		unsigned int triangleIndexTop[3];
		unsigned int triangleIndexBot[3];

		triangleIndexTop[0] = indexTop;
		triangleIndexTop[1] = pointCircleIndexTop[i];
		triangleIndexTop[2] = pointCircleIndexTop[i + 1];

		triangleIndexBot[0] = indexBot;
		triangleIndexBot[1] = pointCircleIndexBot[i];
		triangleIndexBot[2] = pointCircleIndexBot[i + 1];

		unsigned int insertIndexTop = rVisualizer.AddTriangleCell(triangleIndexTop);
		unsigned int insertIndexBot = rVisualizer.AddTriangleCell(triangleIndexBot);

		FullVector<double, 3> dirTop = -n;
		FullVector<double, 3> dirBot = n;
		rVisualizer.SetCellDataVector(insertIndexTop, "Direction", dirTop.data());
		rVisualizer.SetCellDataVector(insertIndexBot, "Direction", dirBot.data());

	}
	// =========================================
	//            create side surface
	// =========================================
	for (int i = 0; i < nPoints; i++)
	{
		unsigned int quadIndex[4];
		quadIndex[0] = pointCircleIndexBot[i];
		quadIndex[1] = pointCircleIndexBot[i + 1];
		quadIndex[2] = pointCircleIndexTop[i + 1];
		quadIndex[3] = pointCircleIndexTop[i];

		unsigned int insertIndex = rVisualizer.AddQuadCell(quadIndex);
		FullVector<double, 3> lineVector = pointCircle.GetRow(i).transpose() - pointCircle.GetRow(i + 1).transpose();
		FullVector<double, 3> dir = lineVector.cross(n);
		dir.normalize();
		rVisualizer.SetCellDataVector(insertIndex, "Direction", dir.data());
	}
}
#endif

bool NuTo::CollidableWallCylinder::IsInside(
		const CollidableParticleSphere& rSphere) const
		{
	FullVector<double, 3> dPos = rSphere.mPosition - mPosition;
	FullVector<double, 3> dP = dPos - dPos.Dot(mDirection) * mDirection;

	return dP.dot(dP) < (mRadius - rSphere.mRadius) * (mRadius - rSphere.mRadius);
}

const bool NuTo::CollidableWallCylinder::IsPhysical() const
{
	return true;
}

void NuTo::CollidableWallCylinder::Print(std::ostream& rReturnStream) const
{
	//TODO: NuTo::CollidableWallCylinder::Print
}
