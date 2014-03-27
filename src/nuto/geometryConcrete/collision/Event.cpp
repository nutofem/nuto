/*
 * EventBase.cpp
 *
 *  Created on: 21 Jan 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include "nuto/geometryConcrete/collision/SubBox.h"

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
	RemoveLocalEvent();
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

double NuTo::Event::AddNewEventsLocalAndGlobal(
		EventListHandler& rEvents) const
		{
	// get boxes involved

	auto boxesFirst = mFirst->GetSubBoxes();
	auto boxesSecond = mSecond->GetSubBoxes();

	double time;

	for (auto box : boxesFirst)
		time += box->CreateEvents(rEvents, *mFirst);

	for (auto box : boxesSecond)
		time += box->CreateEvents(rEvents, *mSecond);

	return time;
}

void NuTo::Event::EraseOldEventsLocalAndGlobal(
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

	if (*this == rOther)
		return false;

	// simultaneous, non-equal events: add somehow distinguish-able

	if (this->mFirst != rOther.mFirst)
		return this->mFirst < rOther.mFirst;

	if (this->mSecond != rOther.mSecond)
		return this->mSecond < rOther.mSecond;

	return this->mFirst < rOther.mSecond;
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

std::size_t NuTo::Event::GetHash() const
{
	std::size_t seed = 0;
	boost::hash_combine(seed, mTime);
	int index1 = mFirst->GetIndex();
	int index2 = mSecond->GetIndex();
	if (index1 < index2){
		int& tmp = index1;
		index1 = index2;
		index2 = tmp;
	}


	boost::hash_combine(seed, index1);
	boost::hash_combine(seed, index2);
	return seed;
}


void NuTo::Event::RemoveLocalEvent()
{
	mFirst->mLocalEvents.remove(this);
	mSecond->mLocalEvents.remove(this);
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
std::size_t hash_value(const NuTo::Event& rEvent)
{
	return rEvent.GetHash();
}

std::ostream& operator <<(std::ostream& rOutStream, const Event& rEvent)
{
	rEvent.Print(rOutStream);
	return rOutStream;
}
} // namespace NuTo

