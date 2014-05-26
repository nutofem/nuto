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
				mEnableStatusVisualization(false),
				mEnableFileOutput(rName != "")

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

void NuTo::CollisionHandler::InitializeLogger(NuTo::Logger& rLogger)
{
	if (mEnableFileOutput)
		rLogger.OpenFile(mName + std::string("/log.dat"));

	rLogger.SetQuiet(false);
	rLogger << "# Event \t Time \t EKin \t Phi \t WTime \n";
}

void NuTo::CollisionHandler::VisualizeSpheres(long rNumEvents, double rGlobalTime, bool rIsFinal)
{
	if (mEnableFileOutput)
		mSpheres->VisualizeSpheres(mName, rNumEvents, rGlobalTime, rIsFinal);
}

void NuTo::CollisionHandler::Simulate(
		const long rNumEventsMax,
		const double rTimeMax,
		const double rWTimeMax,
		const double rTimePrintOut,
		const double rInitialTimeBarrier)
{
	// Stop criteria for the main loop:
	// 1) rWTimeMax reached --> successful!
	// 2) rWTimeMax reached --> failed
	// 3) rNumEventsMax = number of events without change in particle volume fraction
	// means: no change in current global time
	// if dGlobalTime < epsilon in rNumEventsMax events --> simulation is stuck, failed

	const double epsilon = 1e-6;
	long noChangeCounter = 0;
	double noChangeGlobalTime = 0;

	NuTo::Logger logger;
	InitializeLogger(logger);
	Exception caughtException("");

	long numEvents = 0;
	double globalTime = 0.;
	VisualizeSpheres(numEvents, globalTime, false);


	double wTimeEventList = mGlobalEventList.SetTimeBarrier(rInitialTimeBarrier, *mSubBoxes);

	double timePrintOut = rTimePrintOut;

	double globalTimeBarrier = rInitialTimeBarrier;

	double oldGlobalTime = 0.;
	double globalTimeOfBarrierReset = 0.;
	double wallTimeOfBarrierReset = 0.;

	double wStartTime = WallTime::Get();
	while (WallTime::Get() - wStartTime < rTimeMax && globalTime < rWTimeMax)
	{
		globalTime = mGlobalEventList.GetNextEventTime();

		if (globalTime == Event::EVENTNULL)
		{
			globalTimeBarrier = oldGlobalTime + 10 * (oldGlobalTime - globalTimeOfBarrierReset);


			wTimeEventList = mGlobalEventList.SetTimeBarrier(globalTimeBarrier,*mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = WallTime::Get() - wStartTime;
		}

		// reset the time barrier
		if (WallTime::Get() - wStartTime > wallTimeOfBarrierReset + 5. * wTimeEventList)
		{

			globalTimeBarrier = oldGlobalTime + 1.1 * (oldGlobalTime - globalTimeOfBarrierReset);
			wTimeEventList = mGlobalEventList.SetTimeBarrier(globalTimeBarrier,*mSubBoxes);
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
			if (mEnableStatusVisualization && mEnableFileOutput)
				mSpheres->VisualizeSpheres(mName, numEvents, globalTime, false);
		}

		try
		{
			mGlobalEventList.PerformNextEvent();
		} catch (NuTo::Exception& e)
		{
			std::cout << e.ErrorMessage() << std::endl;
			std::cout << " fail timestep: " << numEvents << std::endl;
			caughtException = e;
			break;
		}

		if (globalTime - noChangeGlobalTime < epsilon)
		{
			noChangeCounter++;
			if (noChangeCounter > rNumEventsMax)
			{
				std::stringstream exceptionStream;
				exceptionStream << "[NuTo::CollisionHandler::Simulate] Simulation stopped after " << rNumEventsMax <<
						" events without a significant change.";
				caughtException = NuTo::Exception(exceptionStream.str());
				break;

			}
		}
		else
		{
			noChangeCounter = 0;
			noChangeGlobalTime = globalTime;
		}



		numEvents++;


		oldGlobalTime = globalTime;
	}

	double wTime = WallTime::Get() - wStartTime;

	mGlobalEventList.PrintStatistics(wTime);

	mSpheres->Sync(globalTime);
	LogStatus(logger, numEvents, globalTime, wStartTime);

	VisualizeSpheres(numEvents, globalTime, false);
	VisualizeSpheres(numEvents + 1, globalTime + 1., true);

	logger.CloseFile();

	// rethrow exceptions for proper test failure
	if (caughtException.ErrorMessage() != "")
		throw Exception("[NuTo::CollisionHandler::Simulate] Simulation ended with the exception: \n" + caughtException.ErrorMessage());
}

void NuTo::CollisionHandler::EnableStatusVisualization(bool rEnableStatusVisualization)
{
	mEnableStatusVisualization = rEnableStatusVisualization;
}

void NuTo::CollisionHandler::EnableFileOutput(bool rEnableFileOutput)
{
	mEnableFileOutput = rEnableFileOutput;
}
