/*
 * CollisionHandler.h
 *
 *  Created on: 27 Jan 2014
 *      Author: ttitsche
 */

#ifndef COLLISIONHANDLER_H_
#define COLLISIONHANDLER_H_

#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"

namespace NuTo
{
class SubBoxHandler;
class Logger;

//! @brief ... runs the simulation
class CollisionHandler
{
public:

	//! @brief ... constructor
	//! @param rSpheres ... particles for the simulation
	//! @param rSubBoxes ... sub boxes
	//! @param rName ... work dir
	CollisionHandler(ParticleHandler& rSpheres, SubBoxHandler& rSubBoxes, const std::string rName);

	//! @brief ... true: write a sphere visualization file @ every status output
	void SetStatusVisualization(bool rEnableStatusVisualization);

	//! @brief ... performs the simulation loop
	//! @param rNumEventsMax ... abort if reached
	//! @param rTimeMax ... abort if reached
	//! @param rWTimeMax ... abort if reached
	//! @param rTimePrintOut ... print out every rTimePrintOut seconds
	//! @param rInitialTimeBarrier ... doesnt really matter
	void Simulate(
			const long rNumEventsMax,
			const double rTimeMax,
			const double rWTimeMax,
			const double rTimePrintOut,
			const double rInitialTimeBarrier);


private:

	//! @brief ... writes a sphere visualization file
	//! @param rTimeStep ... current timestep of the simulation
	//! @param rGlobalTime ... current global time != wall time
	//! @param rFinal ... false: use current radius, true: use initial radius
	void VisualizeSpheres(int rTimeStep, double rGlobalTime, bool rFinal);

	//! @brief ... writes an entry to the status file
	void LogStatus(Logger& rLogger, const long rTimeStep, const double rGlobalTime, const double rWSTime) const;

	ParticleHandler* mSpheres;
	SubBoxHandler* mSubBoxes;
	EventListHandler mGlobalEventList;

	//! @brief ... workdir
	const std::string mName;

	//! @brief ... true: write a sphere visualization file @ every status output
	bool mEnableStatusVisualization;

};

} /* namespace NuTo */
#endif /* COLLISIONHANDLER_H_ */
