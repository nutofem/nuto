/*
 * ParticleHandler.h
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#ifndef PARTICLEHANDLER_H_
#define PARTICLEHANDLER_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"

namespace NuTo
{

class CollidableParticleBase;
class CollidableParticleSphere;
class VisualizeUnstructuredGrid;

//! @brief ... handles the particle list
class ParticleHandler
{
public:

	//! @brief ... constructor, builds rNumParticles
	ParticleHandler(
			const int rNumParticles,
			const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox,
			const double rVelocityRange,
			const double rGrowthRate);

	//! @brief ... constructor, uses rSpheres
	ParticleHandler(
			const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres,
			const double rVelocityRange,
			const double rGrowthRate);

	//! @brief ... destructor, deletes all particles
	~ParticleHandler();

	//! @brief ... updates all particles to the same time global time rTime
	void Sync(const double rTime);

	//! @brief ... getter for the kinetic energy of all particles
	const double GetKineticEnergy() const;

	//! @brief ... getter for the volume of all particles
	const double GetVolume() const;

	//! @brief ... writes a sphere visualization file
	//! @param rName ... workdir
	//! @param rTimeStep ... current timestep of the simulation
	//! @param rGlobalTime ... current global time != wall time
	//! @param rFinal ... false: use current radius, true: use initial radius
	void VisualizeSpheres(std::string rName, int rTimeStep, double rGlobalTime, bool rFinal);

	//! @brief ... resets all velocities
	void ResetVelocities();

	//! @brief ... converts the particle list to a Nx4-matrix
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetParticles() const;

	//! @brief ... get a single particle from the particle list
	CollidableParticleSphere* GetParticle(const int rIndex);

	//! @brief ... getter for the particle list size
	const int GetNumParticles() const;

	//! @brief ... calculates the minimal distance between all particles using sub boxes
	double GetAbsoluteMininimalDistance(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox);

	//! @brief ... caluclates approximate sub box length, based on box size
	NuTo::FullVector<int,Eigen::Dynamic> GetSubBoxDivisions(NuTo::FullVector<double, Eigen::Dynamic> rLength);

private:

	CollidableBase::ParticleContainer mParticles;
	int mParticleIndex;

	NuTo::FullVector<double,Eigen::Dynamic> GetRandomVector(const double rStart, const double rEnd);
	NuTo::FullVector<double,Eigen::Dynamic> GetRandomVector(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBounds);

	//! @return ... returns true if all z-positions are equal
	bool Is2DSimulation(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres);


};



} /* namespace NuTo */
#endif /* PARTICLEHANDLER_H_ */
