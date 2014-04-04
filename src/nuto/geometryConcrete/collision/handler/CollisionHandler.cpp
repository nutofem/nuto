/*
 * CollisionHandler.cpp
 *
 *  Created on: 27 Jan 2014
 *      Author: ttitsche
 */

#include <iostream>

#include "nuto/geometryConcrete/WallTime.h"
#include "nuto/geometryConcrete/collision/handler/CollisionHandler.h"
#include "nuto/geometryConcrete/collision/SubBox.h"
#include "nuto/geometryConcrete/collision/handler/SubBoxHandler.h"
#include "nuto/geometryConcrete/collision/handler/ParticleHandler.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include "nuto/base/Logger.h"
#include <sstream>
#include <fstream>

NuTo::CollisionHandler::CollisionHandler(
		ParticleHandler& rSpheres,
		SubBoxHandler& rSubBoxes, const std::string rName)
		: mName(rName),
		  mEnableStatusVisualization(false)

{
	mSpheres = &rSpheres;
	mSubBoxes = &rSubBoxes;
}

void NuTo::CollisionHandler::LogStatus(
		Logger& rLogger,
		const long rTimeStep,
		const double rGlobalTime,
		const double rWSTime) const
{
	double eKin = mSpheres->GetKineticEnergy();
	double vol = mSpheres->GetVolume();
	double time = WallTime::Get() - rWSTime;

	std::stringstream toLog;
	toLog.setf(std::ios::scientific);
	toLog << std::setw(8) << rTimeStep << "\t ";
	toLog << std::setprecision(8) << rGlobalTime << "\t";
	toLog << std::setprecision(2) << eKin << "\t";
	toLog << std::setprecision(4) << vol / mSubBoxes->GetVolume() << "\t";
	toLog << std::setprecision(5) << time << "\n";

	rLogger << toLog.str();

}

void NuTo::CollisionHandler::Simulate(
		const long rNumEventsMax,
		const double rTimeMax,
		const double rWTimeMax,
		const double rTimePrintOut,
		const double rInitialTimeBarrier)
{

	NuTo::Logger logger;
	logger.OpenFile(mName + std::string("/log.dat"));
	logger.SetQuiet(false);
	logger << "# Event \t Time \t EKin \t Phi \t WTime \n";

	long numEvents = 0;
	double globalTime = 0.;
	mSpheres->VisualizeSpheres(mName, numEvents, globalTime, false);

	mGlobalEventList.SetTimeBarrier(rInitialTimeBarrier);
	double wTimeEventList = mGlobalEventList.BuildInitialEventList(*mSubBoxes);


	double timePrintOut = rTimePrintOut;

	double globalTimeBarrier = rInitialTimeBarrier;

	double oldGlobalTime = 0.;
	double globalTimeOfBarrierReset = 0.;
	double wallTimeOfBarrierReset = 0.;

	double wStartTime = WallTime::Get();
	while (WallTime::Get() - wStartTime < rTimeMax && numEvents < rNumEventsMax && globalTime < rWTimeMax)
	{
		globalTime = mGlobalEventList.GetNextEventTime();

		if (globalTime == Event::EVENTNULL)
		{
			globalTimeBarrier = oldGlobalTime + 10 * (oldGlobalTime - globalTimeOfBarrierReset);

			mGlobalEventList.SetTimeBarrier(globalTimeBarrier);
			wTimeEventList = mGlobalEventList.BuildEventList(*mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = WallTime::Get() - wStartTime;
		}

		// reset the time barrier
		if (WallTime::Get() - wStartTime > wallTimeOfBarrierReset + 5. * wTimeEventList)
		{

			globalTimeBarrier = oldGlobalTime + 1.1 * (oldGlobalTime - globalTimeOfBarrierReset);
			mGlobalEventList.SetTimeBarrier(globalTimeBarrier);
			wTimeEventList = mGlobalEventList.BuildEventList(*mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = WallTime::Get() - wStartTime;

		}

		// print a status update
		if (rTimePrintOut != 0 && WallTime::Get() - wStartTime > timePrintOut)
		{
			timePrintOut += rTimePrintOut;
			LogStatus(logger, numEvents, globalTime, wStartTime);
			if(mEnableStatusVisualization)
				mSpheres->VisualizeSpheres(mName, numEvents, globalTime, false);
		}


		try
		{
			mGlobalEventList.PerformNextEvent();
		} catch (NuTo::Exception& e)
		{
			std::cout << e.ErrorMessage() << std::endl;
			std::cout << " fail timestep: " << numEvents << std::endl;
			break;
		}

		numEvents++;

		oldGlobalTime = globalTime;
	}

	double wTime = WallTime::Get() - wStartTime;

	mGlobalEventList.PrintStatistics(wTime);

	mSpheres->Sync(globalTime);
	LogStatus(logger, numEvents, globalTime, wStartTime);

	mSpheres->VisualizeSpheres(mName, numEvents, globalTime, false);
	mSpheres->VisualizeSpheres(mName, numEvents + 1, globalTime + 1, true);

	logger.CloseFile();
}

void NuTo::CollisionHandler::SetStatusVisualization(bool rEnableStatusVisualization)
{
	mEnableStatusVisualization = rEnableStatusVisualization;
}
