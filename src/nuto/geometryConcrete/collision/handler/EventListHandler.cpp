/*
 * EventListHandler.cpp
 *
 *  Created on: 3 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/handler/SubBoxHandler.h"
#include "nuto/geometryConcrete/WallTime.h"

NuTo::EventListHandler::EventListHandler()
		:
				mTimeUpdate(0.),
				mTimeErase(0.),
				mTimeAdd(0.),
				mTimeRebuild(0.),
				mTimeBarrier(0.),
				mNSphereCollisions(0),
				mNWallCollisions(0),
				mNWallTransfers(0)
{
	mEvents.clear();
}


NuTo::EventListHandler::~EventListHandler()
{
	mEvents.clear();
}

double NuTo::EventListHandler::BuildInitialEventList(
		NuTo::SubBoxHandler& rSubBoxes)
{
	std::cout << "BuildInitialEventList() start." << std::endl;

	BuildEventList(rSubBoxes);

	std::cout << "BuildInitialEventList():" << mTimeRebuild << " s | ";
	std::cout << mEvents.size() << " events." << std::endl;

	double rTime = mTimeRebuild;
	mTimeRebuild = 0.;

	return rTime;
}

double NuTo::EventListHandler::BuildEventList(
		NuTo::SubBoxHandler& rSubBoxes)
{

	double sTime = WallTime::Get();

	mEvents.clear();

//	auto boxes = rSubBoxes.GetSubBoxes();
//	for (unsigned int iBox = 0; iBox < boxes.size(); ++iBox)
//	{
//		auto collidables = boxes[iBox]->GetCollidables();
//		for (unsigned int iCollidable = 0; iCollidable < collidables.size(); ++iCollidable)
//		{
//			boxes[iBox]->CreateEvents(*this, *collidables[iCollidable]);
//		}
//	}

	for (auto box : rSubBoxes.GetSubBoxes())
		for (auto coll : box->GetCollidables())
			box->CreateEvents(*this, *coll);

	mTimeRebuild += WallTime::Get() - sTime;

	std::cout << ", " << WallTime::Get() - sTime << " s , " << mEvents.size() << " events." << std::endl;

	return WallTime::Get() - sTime;
}

void NuTo::EventListHandler::PrintEvents()
{
	std::cout << std::endl << "======== GLOBAL EVENT LIST ========" << std::endl;
	for (auto event : mEvents)
		std::cout << event << std::endl;

	std::cout << "===================================" << std::endl << std::endl;
}

void NuTo::EventListHandler::Clear()
{
	mEvents.clear();
}

void NuTo::EventListHandler::AddEvent(const double rTime,
		CollidableBase& rCollidable1, CollidableBase& rCollidable2, int rType)
{
	if (rTime < 0 || rTime > mTimeBarrier)
		return;

	auto insert = mEvents.insert(new Event(rTime, &rCollidable2, &rCollidable1, rType));
	(*insert.first).AddLocalEvent();

//	Event* tmp = new Event(rTime, &rCollidable2, &rCollidable1, rType);
//	tmp->AddLocalEvent();
//	mEvents.push_back(tmp);
}

void NuTo::EventListHandler::DeleteOldEvents(Event::LocalEvents& rOldEvents)
{
	for (unsigned int iOldEvent = 0; iOldEvent < rOldEvents.size(); ++iOldEvent) {
		mEvents.erase(*rOldEvents[iOldEvent]);
	}
}

const double NuTo::EventListHandler::GetNextEventTime()
{
	if (mEvents.begin() == mEvents.end())
		return Event::EVENTNULL;

	return (*mEvents.begin()).GetTime();
}

void NuTo::EventListHandler::PrintStatistics(double rTimeTotal)
{
	std::ostringstream statistics;
	double timeSum = mTimeUpdate + mTimeErase + mTimeAdd + mTimeRebuild;
	statistics.setf(std::ios_base::fixed, std::ios_base::floatfield);
	statistics.precision(4);
	statistics << "Time Update: " << mTimeUpdate << " | " << mTimeUpdate / rTimeTotal * 100 << "%" << std::endl;
	statistics << "Time Erase:  " << mTimeErase << " | " << mTimeErase / rTimeTotal * 100 << "%" << std::endl;
	statistics << "Time Add:    " << mTimeAdd << " | " << mTimeAdd / rTimeTotal * 100 << "%" << std::endl;
	statistics << "Time Rebuild:" << mTimeRebuild << " | " << mTimeRebuild / rTimeTotal * 100 << "%" << std::endl;
	statistics << "Time Sum:    " << timeSum << " | " << timeSum / rTimeTotal * 100 << "%" << std::endl;
	statistics << "Time Total:  " << std::setw(6) << rTimeTotal << std::endl;

	long nEvents = mNSphereCollisions + mNWallCollisions + mNWallTransfers;

	statistics << "Total Events: " << nEvents << std::endl;
	statistics << 100. * mNSphereCollisions / nEvents << "% SphereCollisions" << std::endl;
	statistics << 100. * mNWallCollisions / nEvents << "% WallCollisions" << std::endl;
	statistics << 100. * mNWallTransfers / nEvents << "% WallTransfers" << std::endl;

	std::cout << statistics.str();
}


const int NuTo::EventListHandler::GetEventListSize()
{
	return mEvents.size();
}

void NuTo::EventListHandler::PerformNextEvent()
{

	const Event* nextEvent = new Event(*(mEvents.begin()));

	switch (nextEvent->GetType())
	{
	case CollidableBase::EventType::SphereCollision:
		mNSphereCollisions++;
		break;
	case CollidableBase::EventType::WallCollision:
		mNWallCollisions++;
		break;
	case CollidableBase::EventType::WallTransfer:
		mNWallTransfers++;
		break;
	default:
		break;
	}

	double timeTmp;

	timeTmp = WallTime::Get();
	nextEvent->PerformCollision();
	mTimeUpdate += WallTime::Get() - timeTmp;

	timeTmp = WallTime::Get();
	nextEvent->EraseOldEvents(*this);
	mTimeErase += WallTime::Get() - timeTmp;

	timeTmp = WallTime::Get();
	nextEvent->AddNewEvents(*this);
	mTimeAdd += WallTime::Get() - timeTmp;

	delete nextEvent;

}

void NuTo::EventListHandler::SetTimeBarrier(double rTimeBarrier)
{
	std::stringstream ostream;
	ostream.setf(std::ios_base::scientific);
	ostream << "Reset time barrier to " << std::setprecision(8) << rTimeBarrier;
	std::cout << ostream.str();
	mTimeBarrier = rTimeBarrier;
}
