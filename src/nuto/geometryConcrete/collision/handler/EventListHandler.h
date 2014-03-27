/*
 * EventListHandler.h
 *
 *  Created on: 3 Mar 2014
 *      Author: ttitsche
 */

#ifndef EVENTLISTHANDLER_H_
#define EVENTLISTHANDLER_H_

#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/SubBox.h"
#include <vector>

namespace NuTo
{
class EventListHandler
{
public:
	EventListHandler();

	double BuildInitialEventList(std::vector<SubBox*>& rSubBoxes);
	double BuildEventList(std::vector<SubBox*>& rSubBoxes);

	void PrintEvents();

	void Clear();

	void AddEvent(const double rTime, CollidableBase& rCollidable1, CollidableBase& rCollidable2, int rType);

	void DeleteOldEvents(Event::LocalEvents& rOldEvents);

	const double GetNextEventTime();

	void PerformNextEvent();

	void PrintStatistics(double rTimeTotal);
	void SetTimeBarrier(double timeBarrier);

	const int GetEventListSize();

private:
	Event::GlobalEvents mEvents;

	double mTimeUpdate;
	double mTimeErase;
	double mTimeAdd;
	double mTimeRebuild;
	double mTimeAddDetail;
	double mTimeBarrier;

	long mNSphereCollisions;
	long mNWallCollisions;
	long mNWallTransfers;

};

} /* namespace NuTo */
#endif /* EVENTLISTHANDLER_H_ */
