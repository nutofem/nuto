/*
 * CollidableSphere.h
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#ifndef COLLIDABLEPARTICLESPHERE_H_
#define COLLIDABLEPARTICLESPHERE_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"


namespace NuTo
{

class CollidableBase;
class CollidableWallPhysical;
class CollidableWallVirtual;
class CollidableWallCylinder;
class CollidableWallBase;
class VisualizeUnstructuredGrid;

//! @brief ... class for spherical collidables
class CollidableParticleSphere: public NuTo::CollidableParticleBase
{
	//! @brief provides faster access for collision checks with walls
	friend class CollidableWallPhysical;
	friend class CollidableWallVirtual;
	friend class CollidableWallCylinder;
	friend class CollidableWallBase;

public:

	//! @brief ... constructor. Create CollidableSphere
	//! @param rPosition ... sphere position
	//! @param rVelocity ... sphere velocity
	//! @param rRadius ... sphere radius > 0
	//! @param rGrowthRate ... sphere growth rate > 0
	//! @param rIndex ... name
	CollidableParticleSphere(
			NuTo::FullVector<double, Eigen::Dynamic> rPosition,
			NuTo::FullVector<double, Eigen::Dynamic> rVelocity,
			double rRadius,
			double rGrowthRate,
			const int rIndex);

	//! @brief ... move spheres and apply growth, mTimeOfLastCollision update
	//! @param rTime ... new global time.
	void MoveAndGrow(const double rTime);

	//! @brief ... calculates and returns the kinetic energy of the sphere
	const double GetKineticEnergy() const;

	//! @brief ... calculates and returns the volume of the sphere
	const double GetVolume() const;

	//! @brief ... collision between CollidableSphere and CollidableBase, resolve double dispatch, forward *this
	//! @param rCollidable ... collision partner
	void PerformCollision(CollidableBase& rCollidable);

	//! @brief ... collision between CollidableSphere and CollidableSphere, physics here.
	//! @param rSphere ... collision partner
	void PerformCollision(CollidableParticleSphere& rSphere);

	//! @brief ... collision between CollidableSphere and CollidableWall, forward to CollidableWall.
	//! @param rWall ... collision partner
	void PerformCollision(CollidableWallBase& rWall);

	//! @brief ... collision check between this and a CollidableBase, resolve double dispatch, forward *this
	//! @param rCollidable ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	const double PredictCollision(CollidableBase& rCollidable, int& rType);

	//! @brief ... collision check between this and another CollidableSphere:
	//! Physics lead to quadratic equation.
	//! @param rSphere ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	const double PredictCollision(CollidableParticleSphere& rSphere, int& rType);

	//! @brief ... collision check between CollidableSphere and CollidableWall, forward to CollidableWall.
	//! @param rWall ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	const double PredictCollision(CollidableWallBase& rWall, int& rType);

	//! @brief ... returns all old events, that need to be deleted
	//! @param rEventsToDelete ... return argument
	void GetLocalEventsToDelete(Event::LocalEvents& rEventsToDelete) const;

	//! @brief ... exports the sphere position and its radius to as a row in a Nx4-matrix
	//! @return ... 1x4-matrix, [posX, posY, posZ, radius]
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ExportRow() const;

#ifdef ENABLE_VISUALIZE
	//! @brief ... visualize all moving collidables
	//! @param rVisualizer ... NuTo object for ascii-export
	void VisualizationDynamic(VisualizeUnstructuredGrid& rVisualizer, bool rFinal) const;
#endif

	//! @brief ... resets the sphere velocity to 0.0
	void ResetVelocity();

	//! @brief ... sets a new growth rate
	//! @param rGrowthRateFactor ... growthrate *= rGRFactor
	//! @param rTime ... global time
	void SetGrowthRate(const double rGrowthRateFactor, const double rTime);

	//! @brief ... getter for sphere position
	const NuTo::FullVector<double, Eigen::Dynamic> GetPosition() const;

	//! @brief ... getter for sphere radius
	const double GetRadius() const;

	//! @brief ... getter for initial sphere radius
	const double GetRadius0() const;

private:

	//! @brief ... current sphere radius
	double mRadius;

	//! @brief ... base sphere radius
	double mRadiusGrowth;

	//! @brief ... initial sphere radius
	const double mRadius0;

	//! @brief ... sphere growth rate
	double mGrowthRate;

	//! @brief ... initialized with 0, updated in MoveAndGrow()
	double mTimeOfLastUpdate;

	//! @brief ... initialized with 0, updated in SetGrowthRate()
	double mTimeOfGrowthReset;

	//! @brief ... physics of a one-dimensional, fully elastic sphere collision.
	//! @param rVelocity1 ... velocity of sphere 1
	//! @param rVelocity2 ... velocity of sphere 2
	//! @param rMass1 ... mass of sphere 1
	//! @param rMass2 ... mass of sphere 2
	//! @return ... post-collision velocity for sphere 1.
	const double SphereCollision1D(
			const double rVelocity1,
			const double rVelocity2,
			const double rMass1,
			const double rMass2) const;

	//! @brief ... prints CollidableSphere
	//! @param rReturnStream ... output stream, that gets modified
	void Print(std::ostream & rReturnStream) const;

};

} /* namespace NuTo */
#endif /* COLLIDABLEPARTICLESPHERE_H_ */
