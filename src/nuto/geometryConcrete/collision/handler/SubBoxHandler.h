/*
 * SubBoxBuilder.h
 *
 *  Created on: 4 Mar 2014
 *      Author: ttitsche
 */

#ifndef SUBBOXHANDLER_H_
#define SUBBOXHANDLER_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"

#include "nuto/geometryConcrete/collision/SubBox.h"

namespace NuTo
{

class CollidableWallCylinder;

class SubBoxHandler
{
public:
	SubBoxHandler(ParticleHandler& rSpheres,const FullVector<int, 3>& rDivisions);
	SubBoxHandler(ParticleHandler& rSpheres,const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);
	~SubBoxHandler();

	std::vector<SubBox*> BuildBox(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	std::vector<SubBox*> BuildCylinder(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	std::vector<SubBox*> BuildCylinder2(const double rRadius, const double rHeight);

	double GetVolume() const;

	void PrintBoxes();

	void VisualizeBorders(std::string rFile);
	static FullVector<double, 3> GetLength(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rBox);

	static int GetIndex(int rX, int rY, int rZ, FullVector<int, 3> rDivisions);

private:
	ParticleHandler mSpheres;
	std::vector<SubBox*> mSubBoxes;
	FullVector<int, 3> mDivisions;
	double mVolume;

	void AddSpheresToBoxes(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	std::vector<CollidableWallBase*> GetXYWalls(unsigned int rIndex);

	std::vector<NuTo::SubBox*> GetSubBoxes(
			FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBox);

	std::vector<NuTo::FullVector<double,3> > GetXYCorners(
			std::vector<CollidableWallBase*> rWalls);

	unsigned int GetBoxIndex(int rX, int rY, int rZ);
	FullVector<double, 3> GetSubBoxLength(
			const FullVector<double, 3>& rLength);

	std::vector<int> GetNInside(CollidableWallCylinder* rWallCylinder);
};





} /* namespace NuTo */
#endif /* SUBBOXHANDLER_H_ */
