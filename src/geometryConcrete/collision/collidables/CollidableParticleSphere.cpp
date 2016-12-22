/*
 * CollidableSphere.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#include "math/FullMatrix.h"

#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "geometryConcrete/collision/collidables/CollidableBase.h"
#include "geometryConcrete/collision/collidables/CollidableWallPhysical.h"
#include "geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "geometryConcrete/collision/collidables/CollidableWallVirtual.h"
#include "geometryConcrete/collision/collidables/CollidableWallCylinder.h"

#include "visualize/VisualizeUnstructuredGrid.h"

NuTo::CollidableParticleSphere::CollidableParticleSphere(
		FullVector<double, Eigen::Dynamic> rPosition,
		FullVector<double, Eigen::Dynamic> rVelocity,
		double rRadius,
		double rGrowthRate,
		const int rIndex)
		: CollidableParticleBase(rPosition, rVelocity, rIndex),
				mRadius(rRadius),
				mRadiusGrowth(rRadius),
				mRadius0(rRadius),
				mGrowthRate(rGrowthRate),
				mTimeOfLastUpdate(0.),
				mTimeOfGrowthReset(0.)
{
}

void NuTo::CollidableParticleSphere::MoveAndGrow(const double rTime)
{
	double dTUpdate = rTime - mTimeOfLastUpdate;
	mPosition += mVelocity * dTUpdate;

	double dTGrowth = rTime - mTimeOfGrowthReset;
	mRadius = mRadiusGrowth + mGrowthRate * dTGrowth;

	mTimeOfLastUpdate = rTime;
}

void NuTo::CollidableParticleSphere::SetGrowthRate(const double rGrowthRateFactor, const double rTime)
{
	// calculate new growth radius
	double dTGrowth = rTime - mTimeOfGrowthReset;
	mRadiusGrowth += mGrowthRate * dTGrowth;

	// update fields
	mGrowthRate *= rGrowthRateFactor;
	mTimeOfGrowthReset = rTime;
}

const double NuTo::CollidableParticleSphere::SphereCollision1D(
		const double rVelocity1,
		const double rVelocity2,
		const double rMass1,
		const double rMass2) const
		{
	return (rVelocity1 * (rMass1 - rMass2) + 2 * rVelocity2 * rMass2) / (rMass1 + rMass2);
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableBase& rCollidable)
{
	rCollidable.PerformCollision(*this);
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableParticleSphere& rSphere)
{
	FullVector<double, 3> velocity1 = this->mVelocity;
	FullVector<double, 3> velocity2 = rSphere.mVelocity;

	double mass1 = pow(this->mRadius, 3);
	double mass2 = pow(rSphere.mRadius, 3);

	// 1) Get normal direction
	FullVector<double, 3> n = this->mPosition - rSphere.mPosition;
	n.normalize();

	// 2) Velocity split:
	// 2.1) normal velocity
	double velocityNormal1 = n.dot(velocity1);
	double velocityNormal2 = n.dot(velocity2);

	// 2.2) tangential velocity
	FullVector<double, 3> velocityTransversal1 = velocity1 - n * velocityNormal1;
	FullVector<double, 3> velocityTransversal2 = velocity2 - n * velocityNormal2;

	// 3.1) 1D collision
	double velocityNormalNew1 = SphereCollision1D(velocityNormal1, velocityNormal2, mass1, mass2);
	double velocityNormalNew2 = SphereCollision1D(velocityNormal2, velocityNormal1, mass2, mass1);

	// 3.2) make sure, that the new normal velocities are at least
	// fast as the growth rates --> no simultaneous re-collision

	double velocityNormalEffective = velocityNormalNew1 - velocityNormalNew2;

	if (velocityNormalEffective < 0.)
	{
		// add the missing effective velocity to the spheres
		double velocityDifference1 = -velocityNormalEffective * mass2 / (mass1 + mass2);
		double velocityDifference2 = -velocityNormalEffective - velocityDifference1;

		velocityNormalNew1 += velocityDifference1;
		velocityNormalNew2 -= velocityDifference2;
	}

	// 3.3) add normal velocity boost
	double velocityExtra = 2.0;
	velocityNormalNew1 += velocityExtra * this->mGrowthRate;
	velocityNormalNew2 -= velocityExtra * rSphere.mGrowthRate;

	// 4) combine normal and tangential velocities, add growth rate for collision buffer
	this->mVelocity = (velocityNormalNew1) * n + velocityTransversal1;
	rSphere.mVelocity = (velocityNormalNew2) * n + velocityTransversal2;
}

void NuTo::CollidableParticleSphere::PerformCollision(CollidableWallBase& rWall)
{
	rWall.PerformCollision(*this);
}

const double NuTo::CollidableParticleSphere::PredictCollision(CollidableBase& rCollidable, int& rType)
{
	return rCollidable.PredictCollision(*this, rType);
}

const double NuTo::CollidableParticleSphere::PredictCollision(CollidableParticleSphere& rSphere, int& rType)
{
	if (this == &rSphere)
		return Event::EVENTNULL;

	rType = Event::EventType::SphereCollision;

	// sync both spheres to the more recent time
	double baseTime;

	CollidableParticleSphere* s1 = this;
	CollidableParticleSphere* s2 = &rSphere;

	if (s1->mTimeOfLastUpdate <= s2->mTimeOfLastUpdate)
	{
	    s1 = &rSphere;
	    s2 = this;
	}
	baseTime = s1->mTimeOfLastUpdate;

	FullVector<double, 3> dP = s1->mPosition - (s2->mPosition + s2->mVelocity * (s1->mTimeOfLastUpdate - s2->mTimeOfLastUpdate));
    FullVector<double, 3> dV = s1->mVelocity - s2->mVelocity;

    double dR = s1->mRadius + s2->mRadius + s2->mGrowthRate * (s1->mTimeOfLastUpdate - s2->mTimeOfLastUpdate);
    double dG = s1->mGrowthRate + s2->mGrowthRate;

	double a = dV.dot(dV) - dG * dG;
	double b = 2 * (dV.dot(dP) - dR * dG);
	double c = dP.dot(dP) - dR * dR;

//	std::cout << "a: " << a << " b: " << b << " c: " << c <<std::endl;
	double timeCollision = 0.;

	if (c < -2.e-10 * s1->mRadius)
	{
		std::stringstream exceptionStream;
		exceptionStream << "[NuTo::CollidableSphere::CreateNewEvent] "
				<< "Sphere " << this->mIndex << " overlapping with Sphere "
				<< rSphere.mIndex << "!";
		throw NuTo::Exception(exceptionStream.str());
	}


	double discriminant = b * b - 4 * a * c;

	// allow small negative values ...
	if (discriminant < -1e-12)
		return Event::EVENTNULL;

	// ... an treat them as zero
	if (discriminant < 0.)
		discriminant = 0.;

	if (b < 0.)
	{
		timeCollision = 2 * c / (-b + std::sqrt(discriminant));
		return baseTime + timeCollision;
	}
	else
		if (a < 0. && b >= 0.)
		{
			timeCollision = (-b - std::sqrt(discriminant)) / (2 * a);
			return baseTime + timeCollision;
		}

	return Event::EVENTNULL;
}

const double NuTo::CollidableParticleSphere::PredictCollision(CollidableWallBase& rWall, int& rType)
{
	return rWall.PredictCollision(*this, rType);
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::CollidableParticleSphere::ExportRow(bool rInitialRadius) const
{
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> data(1, 4);
	double radius = rInitialRadius? mRadius0 : mRadius;
	data << mPosition[0], mPosition[1], mPosition[2], radius;
	return data;
}

const NuTo::FullVector<double, Eigen::Dynamic> NuTo::CollidableParticleSphere::GetPosition() const
{
	return mPosition;
}

const double NuTo::CollidableParticleSphere::GetRadius() const
{
	return mRadius;
}

const double NuTo::CollidableParticleSphere::GetRadius0() const
{
	return mRadius0;
}

void NuTo::CollidableParticleSphere::Print(std::ostream& rReturnStream) const
		{
	int indexWidth = 3;
	int precision = 20;
	rReturnStream.setf(std::ios::scientific);
	rReturnStream.precision(precision);
	rReturnStream << "Sphere "
			<< std::setw(indexWidth) << mIndex << " Pos:("
			<< std::setw(precision + 2) << mPosition[0] << ","
			<< std::setw(precision + 2) << mPosition[1] << ","
			<< std::setw(precision + 2) << mPosition[2] << ") Vel:("
			<< std::setw(precision + 2) << mVelocity[0] << ","
			<< std::setw(precision + 2) << mVelocity[1] << ","
			<< std::setw(precision + 2) << mVelocity[2] << ") R:"
			<< std::setw(precision + 2) << mRadius;
}

void NuTo::CollidableParticleSphere::GetLocalEventsToDelete(Event::LocalEvents& rEventsToDelete) const
{
	for (unsigned int iEvent = 0; iEvent < mLocalEvents.size(); ++iEvent) {
		rEventsToDelete.push_back(mLocalEvents[iEvent]);
	}
}

#ifdef ENABLE_VISUALIZE

void NuTo::CollidableParticleSphere::VisualizationDynamic(
		NuTo::VisualizeUnstructuredGrid& rVisualizer,
		bool rFinal) const
		{
	double const* coords = mPosition.data();
	unsigned int index = rVisualizer.AddPoint(coords);
	double radius = rFinal ? mRadius0 : mRadius;
	rVisualizer.SetPointDataScalar(index, "Radius", radius);
	// wierd cast as SetPointDataVector does not handle const data. TODO ?
	// no known conversion for argument 3 from ‘const Scalar* {aka const double*}’ to ‘double*’
	double tmpVelocity[3];
	tmpVelocity[0] = mVelocity[0];
	tmpVelocity[1] = mVelocity[1];
	tmpVelocity[2] = mVelocity[2];
	rVisualizer.SetPointDataVector(index, "Velocity", tmpVelocity);
}
#endif

const double NuTo::CollidableParticleSphere::GetKineticEnergy() const
{
	return .5 * std::pow(mRadius, 3) * mVelocity.dot(mVelocity);
}

const double NuTo::CollidableParticleSphere::GetVolume() const
{
	return 4. / 3. * M_PI * std::pow(mRadius, 3);
}

void NuTo::CollidableParticleSphere::ResetVelocity()
{
//	mVelocity *= 0.5;
	mVelocity << 0, 0, 0;
}
