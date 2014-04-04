/*
 * SubBoxBuilder.h
 *
 *  Created on: 4 Mar 2014
 *      Author: ttitsche
 */

#ifndef SUBBOXHANDLER_H_
#define SUBBOXHANDLER_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/SubBox.h"

namespace NuTo
{
class ParticleHandler;
class CollidableWallCylinder;

//! @brief ... builds and handles sub boxes
class SubBoxHandler
{
public:

	//! @brief ... constructor, initializes the handler with a given number of sub box divisions
	//! @param rSpheres ... spheres to be added to the sub boxes
	//! @param rDivisions ... sets divs[0]*divs[1]*divs[2] sub boxes
	SubBoxHandler(ParticleHandler& rSpheres,const NuTo::FullVector<int, Eigen::Dynamic> rDivisions);

	//! @brief ... constructor, initializes the handler, calculates the number of sub boxes
	//! @param rSpheres ... spheres to be added to the sub boxes
	//! @param rBox ... sets the number of sub boxes according to the rBox size
	SubBoxHandler(ParticleHandler& rSpheres,const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	//! @brief ... deletes all sub boxes
	~SubBoxHandler();

	//! @brief ... build a cubic specimen
	//! @param rBox ... size of the specimen
	//! @param rNumThreads ... number of threads for parallelization
	void BuildBox(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox, const int rNumThreads = 1);

	//! @brief ... build a cylindric specimen
	//! @param rBox ... size of the specimen
	//! @param rNumThreads ... number of threads for parallelization
	void BuildCylinder(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox, const int rNumThreads = 1);

	//! @brief ... getter for specimen volume
	double GetVolume() const;

	//! @brief ... prints all sub boxes and the particles inside of them
	void PrintBoxes();

	//! @brief ... visualizes the sub boxes
	//! @param rFile ... output file
	void VisualizeBorders(std::string rFile);

	//! @brief ... getter for the sub box list
	const std::vector<SubBox*>& GetSubBoxes() const;

	//! @brief ... adds a sub box to the sub box list
	//! @rSubBox ... sub box to add
	void AddSubBox(SubBox& rSubBox);

private:

	ParticleHandler* mSpheres;
	std::vector<SubBox*> mSubBoxes;
	NuTo::FullVector<int, Eigen::Dynamic> mDivisions;
	double mVolume;


	void AddSpheresToBoxes(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	//! @brief ...
	NuTo::FullVector<double, Eigen::Dynamic> GetLength(const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBox);

	//! @brief ...
	std::vector<CollidableWallBase*> GetXYWalls(unsigned int rIndex);

	//! @brief ...
	void BuildSubBoxes(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox,const int rNumThreads);

	//! @brief ...
	std::vector<NuTo::FullVector<double,Eigen::Dynamic> > GetXYCorners(
			std::vector<CollidableWallBase*> rWalls);

	//! @brief ...
	unsigned int GetBoxIndex(int rX, int rY, int rZ);

	//! @brief ...
	NuTo::FullVector<double, Eigen::Dynamic> GetSubBoxLength(
			const NuTo::FullVector<double, Eigen::Dynamic>& rLength);

	//! @brief ...
	std::vector<int> GetNInside(CollidableWallCylinder* rWallCylinder);
};

} /* namespace NuTo */
#endif /* SUBBOXHANDLER_H_ */
