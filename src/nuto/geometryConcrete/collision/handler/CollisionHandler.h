/*
 * CollisionHandler.h
 *
 *  Created on: 27 Jan 2014
 *      Author: ttitsche
 */

#ifndef COLLISIONHANDLER_H_
#define COLLISIONHANDLER_H_

#include <vector>

#include "nuto/geometryConcrete/collision/Event.h"
#include <string.h>
#include "nuto/math/FullVector.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"

namespace NuTo
{

class SubBox;
class Logger;
class CollidableParticleSphere;
class CollisionHandler
{
public:

	CollisionHandler(ParticleHandler rSpheres, std::vector<SubBox*>& rSubBoxes, const std::string rName);

	virtual ~CollisionHandler();

	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Simulate(
			const long rNumEventsMax,
			const double rTimeMax,
			const double rWTimeMax,
			const double rTimePrintOut,
			const double rInitialTimeBarrier);
	void SetVolume(double volume);

private:

	void SyncSpheres(double rTime); // const;

	void VisualizeSpheres(int rTimeStep, double rGlobalTime, bool rFinal);


	void LogStatus(Logger& rLogger, long rTimeStep, double rGlobalTime) const;

	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ExportSpheres();
	NuTo::Logger InitializeLogger();

	ParticleHandler mSpheres;
	EventListHandler mGlobalEventList;


	std::vector<SubBox*> mSubBoxes;
	const std::string mName;
	double mVolume;
	double mWSTime;

};

} /* namespace NuTo */
#endif /* COLLISIONHANDLER_H_ */
