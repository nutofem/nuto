/*
 * SubBoxBuilder.h
 *
 *  Created on: 4 Mar 2014
 *      Author: ttitsche
 */

#pragma once


// member
#include "nuto/geometryConcrete/Specimen.h"
#include "nuto/geometryConcrete/collision/SubBox.h"

namespace NuTo
{
class CollidableWallBase;
class CollidableWallCylinder;
class ParticleHandler;


//! @brief builds and handles sub boxes
class SubBoxHandler
{
public:
    //! @brief constructor, initializes the handler with a given number of sub box divisions
    //! @param rSpheres spheres to be added to the sub boxes
    //! @param rDivisions sets divs[0]*divs[1]*divs[2] sub boxes
    SubBoxHandler(ParticleHandler& rSpheres, Specimen& rSpecimen, Eigen::Vector3i rDivisions);

    //! @brief constructor, initializes the handler, calculates the number of sub boxes
    //! @param rSpheres spheres to be added to the sub boxes
    SubBoxHandler(ParticleHandler& rSpheres, Specimen& rSpecimen, int rParticlesPerSubBox);

    //! @brief getter for specimen volume
    double GetVolume() const;

    //! @brief prints all sub boxes and the particles inside of them
    void PrintBoxes();

    //! @brief visualizes the sub boxes
    //! @param rFile output file
    void VisualizeBorders(std::string rFile);

    //! @brief getter for the sub box list
    std::vector<SubBox>& GetSubBoxes();

    //! @brief adds a sub box to the sub box list
    //! @rSubBox sub box to add
    void AddSubBox(SubBox& rSubBox);

private:
    ParticleHandler* mSpheres;
    std::vector<SubBox> mSubBoxes;
    Specimen mSpecimen;
    Eigen::Vector3i mDivisions;

    void Build();

    //! @brief build a cubic specimen
    void BuildBox();

    //! @brief build a cylindric specimen
    void BuildCylinder();

    void AddSpheresToBoxes();

    //! @brief ...
    std::vector<CollidableWallBase*> GetXYWalls(unsigned int rIndex);

    //! @brief ...
    void BuildSubBoxes();

    //! @brief ...
    std::vector<Eigen::VectorXd> GetXYCorners(std::vector<CollidableWallBase*> rWalls);

    //! @brief ...
    unsigned int GetBoxIndex(int rX, int rY, int rZ);

    //! @brief ...
    Eigen::Vector3d GetSubBoxLength();

    //! @brief ...
    std::vector<int> GetNInside(CollidableWallCylinder* rWallCylinder);
};

} /* namespace NuTo */
