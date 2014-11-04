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
class Specimen;

//! @brief ... handles the particle list
class ParticleHandler
{
public:


    //! @brief ... constructor, builds rNumParticles equal particles
	ParticleHandler(
			const int rNumParticles,
			const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rParticleBoundingBox,
			const double rVelocityRange,
			const double rGrowthRate);

	//! @brief ... constructor, uses rSpheres
	ParticleHandler(
			const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres,
			const double rVelocityRange,
			const double rRelativeGrowthRate,
			const double rAbsoluteGrowthRate);

    //! @brief ... constructor, builds particles from external file
    ParticleHandler(
            const std::string& rFileName,
            const double rVelocityRange,
            const double rRelativeGrowthRate,
            const double rAbsoluteGrowthRate);



	//! @brief ... destructor, deletes all particles
	~ParticleHandler();

	//! @brief ... updates all particles to the same time global time rTime
	void Sync(const double rTime);

	//! @brief ... getter for the kinetic energy of all particles
	const double GetKineticEnergy() const;

	//! @brief ... getter for the volume of all particles
	const double GetVolume() const;

	//! @brief ... writes a sphere visualization file
    //! @param rOutputDirectory ... workdir
	//! @param rTimeStep ... current timestep of the simulation
	//! @param rGlobalTime ... current global time != wall time
	//! @param rFinal ... false: use current radius, true: use initial radius
    void ExportParticlesToVTU3D(std::string rOutputDirectory, int rTimeStep, double rGlobalTime, bool rFinal);

    //! @brief ... writes a sphere visualization file
    //! @param rOutputDirectory ... workdir
    //! @param rFinal ... false: use current radius, true: use initial radius
    void ExportParticlesToVTU2D(std::string rOutputFile, double rZCoord);


	//! @brief ... resets all velocities
	void ResetVelocities();

	//! @brief ... converts the particle list to a Nx4-matrix
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetParticles(bool rInitialRadius = false) const;

    //! @brief ... exports the particle list to a file
    void ExportParticlesToFile(const std::string& rExportFileName, bool rInitialRadius) const;

	//! @brief ... get a single particle from the particle list
	CollidableParticleSphere* GetParticle(const int rIndex);

	//! @brief ... getter for the particle list size
	const int GetNumParticles() const;

	//! @brief ... calculates the minimal distance between all particles using sub boxes
	double GetAbsoluteMininimalDistance(Specimen& rSpecimen);

	//! @brief ... optional: change the file name, default: "spheres_"
    void SetVisualizationFileName(const std::string& rVisualizationFileName);

	//! @brief ... calculates approximate sub box length, based on box size and the number of particles per sub box
	NuTo::FullVector<int,Eigen::Dynamic> GetSubBoxDivisions(Specimen& rSpecimen, const int rParticlesPerBox);

    //! @brief ... cut spheres at a given z-coordinate to create circles (in 2D)
    //! @param rZCoord z coordinate (where to cut)
    //! @param rMinRadius minimal radius of the circle
    //! @return ... matrix with the circles (x,y,r)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> CutSpheresZ(
            double rZCoord, double rMinRadius) const;

private:

	CollidableBase::ParticleContainer mParticles;

	//! @brief ... returns a random vector with each component in a certain range
	//! @param rStart ... start of value range
	//! @param rEnd ... end of value range
	NuTo::FullVector<double,Eigen::Dynamic> GetRandomVector(const double rStart, const double rEnd);

	//! @brief ... returns a random vector with each component in a certain range
	//! @param rBounds ... rBounds(:,0) start of value range, rBounds(:,1) end of value range
	NuTo::FullVector<double,Eigen::Dynamic> GetRandomVector(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBounds);

    void CreateParticlesFromMatrix(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres, const double rVelocityRange, const double rRelativeGrowthRate, const double rAbsoluteGrowthRate);

	std::string mVisualizationFileName;

};



} /* namespace NuTo */
#endif /* PARTICLEHANDLER_H_ */
