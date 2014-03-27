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
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include "nuto/base/Logger.h"
#include <sstream>
#include <fstream>

NuTo::CollisionHandler::CollisionHandler(
		ParticleHandler rSpheres,
		std::vector<SubBox*>& rSubBoxes, const std::string rName)
		: mSpheres(rSpheres), mName(rName)

{
	mSubBoxes = rSubBoxes;
	mWSTime = 0.;
	mVolume = 0.;
}


NuTo::CollisionHandler::~CollisionHandler()
{
	// clear event list before spheres!
	// delete event --> remove from local sphere event lists
	mGlobalEventList.Clear();
}

void NuTo::CollisionHandler::SetVolume(double rVolume)
{
	mVolume = rVolume;
}

NuTo::Logger NuTo::CollisionHandler::InitializeLogger()
{
	NuTo::Logger logger;
	logger.OpenFile(mName + std::string("/log.dat"));
	logger.SetQuiet(false);
	logger << "# Event \t Time \t EKin \t Phi \t WTime \n";
	return logger;
}

void NuTo::CollisionHandler::LogStatus(Logger& rLogger, long rTimeStep, double rGlobalTime) const
{
	double eKin = mSpheres.GetKineticEnergy();
	double vol = mSpheres.GetVolume();
	double time = WallTime::Get() - mWSTime;

	std::stringstream toLog;
	toLog.setf(std::ios::scientific);
	toLog << std::setw(8) << rTimeStep << "\t ";
	toLog << std::setprecision(8) << rGlobalTime << "\t";
	toLog << std::setprecision(2) << eKin << "\t";
	toLog << std::setprecision(4) << vol / mVolume << "\t";
	toLog << std::setprecision(5) << time << "\n";

	rLogger << toLog.str();
}

NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::CollisionHandler::Simulate(
		const long rNumEventsMax,
		const double rTimeMax,
		const double rWTimeMax,
		const double rTimePrintOut,
		const double rInitialTimeBarrier)
{

	NuTo::Logger logger = InitializeLogger();


	long numEvents = 0;
	double globalTime = 0.;
	mSpheres.VisualizeSpheres(mName, numEvents, globalTime, false);

	mGlobalEventList.SetTimeBarrier(rInitialTimeBarrier);
	double wTimeEventList = mGlobalEventList.BuildInitialEventList(mSubBoxes);


	double timePrintOut = rTimePrintOut;

	double globalTimeBarrier = rInitialTimeBarrier;

	double oldGlobalTime = 0.;
	double globalTimeOfBarrierReset = 0.;
	double wallTimeOfBarrierReset = 0.;

	mWSTime = WallTime::Get();
	while (WallTime::Get() - mWSTime < rTimeMax && numEvents < rNumEventsMax && globalTime < rWTimeMax)
	{
		globalTime = mGlobalEventList.GetNextEventTime();

		// time barrier reached
		if (globalTime == Event::EVENTNULL)
		{
			globalTimeBarrier = oldGlobalTime + 10 * (oldGlobalTime - globalTimeOfBarrierReset);

			mGlobalEventList.SetTimeBarrier(globalTimeBarrier);
			wTimeEventList = mGlobalEventList.BuildEventList(mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = WallTime::Get() - mWSTime;
		}

		// reset the time barrier
		if (WallTime::Get() - mWSTime > wallTimeOfBarrierReset + 10. * wTimeEventList)
		{

			globalTimeBarrier = oldGlobalTime + 1.1 * (oldGlobalTime - globalTimeOfBarrierReset);
			mGlobalEventList.SetTimeBarrier(globalTimeBarrier);
			wTimeEventList = mGlobalEventList.BuildEventList(mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = WallTime::Get() - mWSTime;

		}

		// print a status update
		if (rTimePrintOut != 0 && WallTime::Get() - mWSTime > timePrintOut)
		{
			timePrintOut += rTimePrintOut;
//			mSpheres.VisualizeSpheres(mName, numEvents, globalTime, false);
			LogStatus(logger, numEvents, globalTime);
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

	double wTime = WallTime::Get() - mWSTime;

	mGlobalEventList.PrintStatistics(wTime);

	mSpheres.Sync(globalTime);
	LogStatus(logger, numEvents, globalTime);

	mSpheres.VisualizeSpheres(mName, numEvents, globalTime, false);
	mSpheres.VisualizeSpheres(mName, numEvents + 1, globalTime + 1, true);

	logger.CloseFile();

	mGlobalEventList.Clear();

	return mSpheres.GetParticles();
}

