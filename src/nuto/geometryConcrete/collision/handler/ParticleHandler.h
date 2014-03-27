/*
 * ParticleHandler.h
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#ifndef PARTICLEHANDLER_H_
#define PARTICLEHANDLER_H_

#include <vector>
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#endif

namespace NuTo
{

class ParticleHandler
{
public:
	ParticleHandler(ParticleHandler& rOther);
	ParticleHandler(
			const int rNumParticles,
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox,
			const double rVelocityRange,
			const double rGrowthRate);
	ParticleHandler(
			const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rSpheres,
			const double rVelocityRange,
			const double rGrowthRate);

	~ParticleHandler();

	void Delete();

	void Sync(const double rTime);
	const double GetKineticEnergy() const;
	const double GetVolume() const;

	void VisualizeSpheres(std::string rName, int rTimeStep, double rGlobalTime, bool rFinal);

	void ResetVelocities();

	void AddParticle(const CollidableParticleSphere& rSphere);

	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetParticles() const;

	CollidableParticleSphere* GetParticle(const int rIndex);

	const int GetNumParticles() const;

	double GetAbsoluteMininimalDistance(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBoundingBox);

	double GetRMax();

private:

	FullVector<double,3> GetRandomVector(const double rStart, const double rEnd);
	FullVector<double,3> GetRandomVector(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rBounds);


	CollidableBase::ParticleContainer mParticles;
	int mParticleIndex;

};

class StatBox
{
public:
	StatBox(FullMatrix<double, 3, 2> rBox):mBox(rBox){};
	std::vector<int> mSphereIndices;
	FullMatrix<double, 3, 2> mBox;
};

} /* namespace NuTo */
#endif /* PARTICLEHANDLER_H_ */
