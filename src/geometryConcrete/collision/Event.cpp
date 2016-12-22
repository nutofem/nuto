/*
 * EventBase.cpp
 *
 *  Created on: 21 Jan 2014
 *      Author: ttitsche
 */

#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/collidables/CollidableBase.h"
#include "geometryConcrete/collision/handler/EventListHandler.h"
#include "geometryConcrete/collision/SubBox.h"
#include <algorithm>
#include <iomanip>

const double NuTo::Event::EVENTNULL = -1.;

NuTo::Event::Event(const double rTime, CollidableBase* rFirst,
		CollidableBase* rSecond, const int rType)
		: mTime(rTime), mType(rType)
{
	mFirst = rFirst;
	mSecond = rSecond;
}

NuTo::Event::Event(const Event& rEvent)
		: mTime(rEvent.mTime), mType(rEvent.mType)
{
	mFirst = rEvent.mFirst;
	mSecond = rEvent.mSecond;
}

NuTo::Event::~Event()
{
	auto& local1st = mFirst->mLocalEvents;
	auto& local2nd = mSecond->mLocalEvents;

	auto newEnd = std::remove(local1st.begin(), local1st.end(), this);
	local1st.erase(newEnd, local1st.end());

	newEnd = std::remove(local2nd.begin(), local2nd.end(), this);
	local2nd.erase(newEnd, local2nd.end());
}

const double NuTo::Event::GetTime() const
{
	return mTime;
}

void NuTo::Event::PerformCollision() const
{
	mFirst->MoveAndGrow(mTime);
	mSecond->MoveAndGrow(mTime);
	mFirst->PerformCollision(*mSecond);
}

void NuTo::Event::AddNewEvents(EventListHandler& rEvents) const
		{
	// get boxes involved

	auto& boxesFirst = mFirst->GetSubBoxes();
	auto& boxesSecond = mSecond->GetSubBoxes();

	for(unsigned int iBox = 0; iBox < boxesFirst.size(); ++iBox)
		boxesFirst[iBox]->CreateEvents(rEvents, *mFirst);

	for(unsigned int iBox = 0; iBox < boxesSecond.size(); ++iBox)
		boxesSecond[iBox]->CreateEvents(rEvents, *mSecond);
}

void NuTo::Event::EraseOldEvents(
		EventListHandler& rEvents) const
		{

	// rEvents.erase() erases and deletes an event its destructor removes it from the local lists,
	// which will mess up a loop trough the local lists.

	// "toDelete"-list for temporary old event storage
	LocalEvents toDelete;

	mFirst->GetLocalEventsToDelete(toDelete);
	mSecond->GetLocalEventsToDelete(toDelete);

	rEvents.DeleteOldEvents(toDelete);
}

bool NuTo::Event::operator <(const Event& rOther) const
		{
	// most likely case:
	if (this->mTime != rOther.mTime)
		return this->mTime < rOther.mTime;


	// simultaneous events!

	// avoid the same event to be added twice, same event might be: this->mFirst == rOther.mSecond

//	if (*this == rOther)
//		return false;

	CollidableBase* smaller1 = this->mFirst;
	CollidableBase* bigger1 = this->mSecond;

	CollidableBase* smaller2 = rOther.mFirst;
	CollidableBase* bigger2 = rOther.mSecond;

	if (smaller1 > bigger1)
		std::swap(smaller1, bigger1);

	if(smaller2 > bigger2)
		std::swap(smaller2, bigger2);


	// compare smaller pointer:
	if (smaller1 == smaller2)
		return bigger1 > bigger2;

	return smaller1 > smaller2;


	// simultaneous, non-equal events: add somehow distinguish-able

//	if (this->mFirst != rOther.mFirst)
//		return this->mFirst < rOther.mFirst;
//
//	if (this->mSecond != rOther.mSecond)
//		return this->mSecond > rOther.mSecond;
//
//	return this->mFirst < rOther.mSecond;
}

bool NuTo::Event::operator ==(Event const& rRhs) const
		{
	// return false if not simultaneous
	if (mTime != rRhs.mTime)
		return false;

	// allow twisted members
	if (mFirst == rRhs.mFirst && mSecond == rRhs.mSecond)
		return true;

	if (mFirst == rRhs.mSecond && mSecond == rRhs.mFirst)
		return true;

	return false;
}

bool NuTo::Event::operator !=(const Event& rRhs) const
		{
	return !(*this == rRhs);
}

void NuTo::Event::AddLocalEvent()
{
	mFirst->mLocalEvents.push_back(this);
	mSecond->mLocalEvents.push_back(this);
}

const int NuTo::Event::GetType() const
{
	return mType;
}

void NuTo::Event::Print(std::ostream& rOutStream) const
		{
	int indexWidth = 6;
	rOutStream.setf(std::ios::scientific);
	rOutStream << "Time="
			<< std::setprecision(15) << mTime << " | "
			<< std::setw(indexWidth) << mFirst->GetIndex() << " vs. "
			<< std::setw(indexWidth) << mSecond->GetIndex();
}

namespace NuTo
{

std::ostream& operator <<(std::ostream& rOutStream, const Event& rEvent)
{
	rEvent.Print(rOutStream);
	return rOutStream;
}
} // namespace NuTo

