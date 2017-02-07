/*
 * CollisionHandler.cpp
 *
 *  Created on: 27 Jan 2014
 *      Author: ttitsche
 */

#include <iostream>
#include <iomanip>
#include "base/Timer.h"

#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/handler/CollisionHandler.h"
#include "geometryConcrete/collision/SubBox.h"
#include "geometryConcrete/collision/handler/SubBoxHandler.h"
#include "geometryConcrete/collision/handler/ParticleHandler.h"
#include "visualize/VisualizeUnstructuredGrid.h"
#include "base/Logger.h"
#include <sstream>
#include <fstream>

NuTo::CollisionHandler::CollisionHandler(
		ParticleHandler& rSpheres,
		SubBoxHandler& rSubBoxes, const std::string rName)
        : mName(rName),
				mEnableFileOutput(rName != "")

{
	mSpheres = &rSpheres;
	mSubBoxes = &rSubBoxes;
}

void NuTo::CollisionHandler::LogStatus(
		Logger& rLogger,
		const long rTimeStep,
		const double rGlobalTime,
		const double rTime) const
		{
	double eKin = mSpheres->GetKineticEnergy();
	mSpheres->Sync(rGlobalTime);
	double vol = mSpheres->GetVolume();

	std::stringstream toLog;
	toLog.setf(std::ios::scientific);
	toLog << std::setw(8) << rTimeStep << "\t ";
	toLog << std::setprecision(8) << rGlobalTime << "\t";
	toLog << std::setprecision(2) << eKin << "\t";
	toLog << std::setprecision(4) << vol / mSubBoxes->GetVolume() << "\t";
	toLog << std::setprecision(5) << rTime << "\n";

	rLogger << toLog.str();

}

void NuTo::CollisionHandler::InitializeLogger(NuTo::Logger& rLogger)
{
	if (mEnableFileOutput)
		rLogger.OpenFile(mName + std::string("_log.dat"));

	rLogger.SetQuiet(false);
	rLogger << "# Event \t Time \t EKin \t Phi \t WTime \n";
}

void NuTo::CollisionHandler::VisualizeSpheres(long rNumEvents, double rGlobalTime, bool rIsFinal)
{
    // mSpheres->VisualizeSpheres(mName, rNumEvents, rGlobalTime, rIsFinal);
}

double NuTo::CollisionHandler::Simulate(
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
	double globalTimePrint = 0.;

	double wTimeEventList = mGlobalEventList.SetTimeBarrier(rInitialTimeBarrier, *mSubBoxes);

	double timePrintOut = rTimePrintOut;

	double globalTimeBarrier = rInitialTimeBarrier;

	double oldGlobalTime = 0.;
	double globalTimeOfBarrierReset = 0.;
	double wallTimeOfBarrierReset = 0.;

	Timer timer("", false);
    while (timer.GetTimeDifference() < rWTimeMax && globalTime < rTimeMax)
	{
		globalTime = mGlobalEventList.GetNextEventTime();
		if (globalTime == Event::EVENTNULL)
		{
			globalTimeBarrier = oldGlobalTime + 10 * (oldGlobalTime - globalTimeOfBarrierReset);


			wTimeEventList = mGlobalEventList.SetTimeBarrier(globalTimeBarrier,*mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = timer.GetTimeDifference();
		}

		// reset the time barrier
		if (timer.GetTimeDifference() > wallTimeOfBarrierReset + 5. * wTimeEventList)
		{

			globalTimeBarrier = oldGlobalTime + 1.1 * (oldGlobalTime - globalTimeOfBarrierReset);
			wTimeEventList = mGlobalEventList.SetTimeBarrier(globalTimeBarrier,*mSubBoxes);
			wTimeEventList = std::max(wTimeEventList, .05);
			globalTime = mGlobalEventList.GetNextEventTime();

			globalTimeOfBarrierReset = globalTime;
			wallTimeOfBarrierReset = timer.GetTimeDifference();

		}

//		mSpheres->Sync(globalTime);
//		VisualizeSpheres(numEvents, globalTime, false);

		// print a status update
		// a) at time print out
		bool statusPrintOut = rTimePrintOut != 0 && timer.GetTimeDifference() > timePrintOut;
		// b) at every 0.1 global time steps
		bool statusGlobalTime = globalTime > globalTimePrint;

		if (statusPrintOut or statusGlobalTime)
		{
			LogStatus(logger, numEvents, globalTime, timer.GetTimeDifference());

	        if (statusGlobalTime)
	            globalTimePrint = globalTime + 0.1;

	        if (statusPrintOut)
	            timePrintOut += rTimePrintOut;

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

    const double finalTimeDifference = timer.GetTimeDifference();

	mGlobalEventList.PrintStatistics(finalTimeDifference);

	mSpheres->Sync(globalTime);
	LogStatus(logger, numEvents, globalTime, finalTimeDifference);

	logger.CloseFile();


	// rethrow exceptions for proper test failure
	if (caughtException.ErrorMessage() != "")
		throw Exception(__PRETTY_FUNCTION__, "Simulation ended with the exception: \n" + caughtException.ErrorMessage());

    return finalTimeDifference;
}

void NuTo::CollisionHandler::EnableFileOutput(bool rEnableFileOutput)
{
	mEnableFileOutput = rEnableFileOutput;
}
