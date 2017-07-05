/*
 * ParticleHandler.h
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#pragma once


#include <vector>
#include <random>
#include <eigen3/Eigen/Core>

namespace NuTo
{

class CollidableParticleBase;
class CollidableParticleSphere;
class VisualizeUnstructuredGrid;
class Specimen;
typedef std::vector<CollidableParticleSphere*> ParticleContainer;

//! @brief ... handles the particle list
class ParticleHandler
{
public:
    //! @brief ... constructor, builds rNumParticles equal particles
    ParticleHandler(int rNumParticles, Eigen::MatrixXd rParticleBoundingBox, double rVelocityRange, double rGrowthRate,
                    int rSeed = 0);

    //! @brief ... constructor, uses rSpheres
    ParticleHandler(Eigen::MatrixXd rSpheres, double rVelocityRange, double rRelativeGrowthRate,
                    double rAbsoluteGrowthRate, int rSeed = 0);

    //! @brief ... constructor, builds particles from external file
    ParticleHandler(const std::string& rFileName, double rVelocityRange, double rRelativeGrowthRate,
                    double rAbsoluteGrowthRate, int rSeed = 0);


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
    void ExportParticlesToVTU3D(std::string rOutputDirectory, int rTimeStep, double rGlobalTime, bool rFinal) const;

    //! @brief ... writes a sphere visualization file
    //! @param rOutputDirectory ... workdir
    //! @param rFinal ... false: use current radius, true: use initial radius
    void ExportParticlesToVTU2D(std::string rOutputFile, double rZCoord) const;

    //! @brief ... writes a gmsh .geo file 3D
    //! @param rOutputFile ... .geo-file name
    //! @param rSpecimen ... specimen
    //! @param rMeshSize ... mesh size
    void ExportParticlesToGmsh3D(std::string rOutputFile, Specimen& rSpecimen, double rMeshSize) const;

    //! @brief ... writes a gmsh .geo file 3D
    //! @param rOutputFile ... .geo-file name
    //! @param rSpecimen ... specimen
    //! @param rMeshSize ... mesh size
    //! @param rZCoord ... z coordinate (where to cut)
    //! @param rMinRadius ... minimal radius of the circle
    void ExportParticlesToGmsh2D(std::string rOutputFile, Specimen& rSpecimen, double rMeshSize, double rZCoord,
                                 double rMinRadius) const;


    //! @brief ... resets all velocities
    void ResetVelocities();

    //! @brief ... converts the particle list to a Nx4-matrix
    Eigen::MatrixXd GetParticles(bool rInitialRadius = false) const;

    //! @brief ... cut spheres at a given z-coordinate to create circles (in 2D)
    //! @param rZCoord z coordinate (where to cut)
    //! @param rMinRadius minimal radius of the circle
    //! @return ... matrix with the circles (x,y,r)
    Eigen::MatrixXd GetParticles2D(double rZCoord, double rMinRadius) const;

    //! @brief ... exports the particle list to a file
    void ExportParticlesToFile(const std::string& rExportFileName, bool rInitialRadius) const;

    //! @brief ... get a single particle from the particle list
    CollidableParticleSphere* GetParticle(const int rIndex) const;

    //! @brief ... getter for the particle list size
    int GetNumParticles() const;

    //! @brief ... calculates the minimal distance between all particles using sub boxes
    double GetAbsoluteMininimalDistance(Specimen& rSpecimen);

    //! @brief ... optional: change the file name, default: "spheres_"
    void SetVisualizationFileName(const std::string& rVisualizationFileName);

    //! @brief ... calculates approximate sub box length, based on box size and the number of particles per sub box
    Eigen::VectorXi GetSubBoxDivisions(Specimen& rSpecimen, const int rParticlesPerBox);

private:
    ParticleContainer mParticles;

    //! @brief ... returns a random vector with each component in a certain range
    //! @param rStart ... start of value range
    //! @param rEnd ... end of value range
    Eigen::VectorXd GetRandomVector(double rStart, double rEnd);

    //! @brief ... returns a random vector with each component in a certain range
    //! @param rBounds ... rBounds(:,0) start of value range, rBounds(:,1) end of value range
    Eigen::VectorXd GetRandomVector(const Eigen::MatrixXd& rBounds);

    void CreateParticlesFromMatrix(const Eigen::MatrixXd& rSpheres, double rVelocityRange, double rRelativeGrowthRate,
                                   double rAbsoluteGrowthRate);

    std::string mVisualizationFileName;

    std::mt19937 mRNG;
};


} /* namespace NuTo */
