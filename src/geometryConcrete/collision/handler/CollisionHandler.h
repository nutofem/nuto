/*
 * CollisionHandler.h
 *
 *  Created on: 27 Jan 2014
 *      Author: ttitsche
 */

#pragma once


#include "geometryConcrete/collision/handler/EventListHandler.h"


namespace NuTo
{
class ParticleHandler;
class SubBoxHandler;
class Logger;

//! @brief runs the simulation
class CollisionHandler
{
public:
    //! @brief constructor
    //! @param rSpheres particles for the simulation
    //! @param rSubBoxes sub boxes
    //! @param rName work dir
    CollisionHandler(ParticleHandler& rSpheres, SubBoxHandler& rSubBoxes, const std::string rName);

    void EnableFileOutput(bool rEnableFileOutput);

    //! @brief performs the simulation loop
    //! @param rNumEventsMax abort if reached
    //! @param rTimeMax abort if reached
    //! @param rWTimeMax abort if reached
    //! @param rTimePrintOut print out every rTimePrintOut seconds
    //! @param rInitialTimeBarrier doesnt really matter
    double Simulate(const long rNumEventsMax, const double rTimeMax, const double rWTimeMax, const double rTimePrintOut,
                    const double rInitialTimeBarrier);

private:
    //! @brief writes an entry to the status file
    void LogStatus(Logger& rLogger, const long rTimeStep, const double rGlobalTime, const double rWSTime) const;

    //! @brief initializes the logger and writes a header
    void InitializeLogger(NuTo::Logger& rLogger);

    ParticleHandler* mSpheres;
    SubBoxHandler* mSubBoxes;
    EventListHandler mGlobalEventList;

    //! @brief workdir
    const std::string mName;

    bool mEnableFileOutput;
};

} /* namespace NuTo */
