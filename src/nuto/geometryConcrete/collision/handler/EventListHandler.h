/*
 * EventListHandler.h
 *
 *  Created on: 3 Mar 2014
 *      Author: ttitsche
 */

#ifndef EVENTLISTHANDLER_H_
#define EVENTLISTHANDLER_H_

#include "nuto/geometryConcrete/collision/Event.h"
#include <vector>

namespace NuTo
{
class SubBoxHandler;

//! @brief ... class for event list operations
class EventListHandler
{
public:

	//! @brief ... constructor, provides an empty event list and initializes timers
	EventListHandler();

	//! @brief ... destructor, deletes all events
	~EventListHandler();

	//! @brief ... builds the initial event list with additional statistics
	double BuildInitialEventList(SubBoxHandler& rSubBoxes);

	//! @brief ... rebuilds the event list
	double BuildEventList(SubBoxHandler& rSubBoxes);

	//! @brief ... performs the most recent event
	//! ... 1) get next events and performs it
	//! ... 2) removes all old events
	//! ... 3) finds new events
	void PerformNextEvent();

	//! @brief ... prints all events
	void PrintEvents();


	//! @brief ... adds a new event, of rTime < mTimeBarrier
	//! @param rTime ... collision time
	//! @param rCollidable1 ... collision partner 1
	//! @param rCollidable2 ... collision partner 2
	//! @param rType ... type of the collision
	void AddEvent(const double rTime, CollidableBase& rCollidable1, CollidableBase& rCollidable2, int rType);

	//! @brief ... deletes every event in rOldEvents
	//! @param rOldEvents ... old events that need to be deleted
	void DeleteOldEvents(Event::LocalEvents& rOldEvents);

	//! @brief ... returns the time of the most recent event
	const double GetNextEventTime();

	//! @brief ... setter for the time barrier
	void SetTimeBarrier(double rTimeBarrier);

	//! @brief getter for the event list size
	const int GetEventListSize();

	//! @brief ... deletes all events
	void Clear();

	//! @brief ... prints all time mesurement statistics
	void PrintStatistics(double rTimeTotal);
private:
	Event::GlobalEvents mEvents;

	double mTimeUpdate;
	double mTimeErase;
	double mTimeAdd;
	double mTimeRebuild;
	double mTimeBarrier;

	long mNSphereCollisions;
	long mNWallCollisions;
	long mNWallTransfers;

};

} /* namespace NuTo */
#endif /* EVENTLISTHANDLER_H_ */
