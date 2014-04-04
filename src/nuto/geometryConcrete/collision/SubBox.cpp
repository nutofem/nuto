/*
 * CellGridBox.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/WallTime.h"
#include "nuto/geometryConcrete/collision/SubBox.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include <omp.h>
#include <algorithm>

NuTo::SubBox::SubBox(const int rIndex, const int rNumThreads)
		: mIndex(rIndex), mNumThreads(rNumThreads)
{
}

NuTo::SubBox::~SubBox()
{
	for (auto it = mWalls.begin(); it != mWalls.end(); ++it)
		delete *it;
}
void NuTo::SubBox::AddSphere(CollidableParticleSphere& rSphere)
{
	mCollidables.push_back(&rSphere);
	rSphere.AddBox(*this);
}

void NuTo::SubBox::RemoveSphere(CollidableParticleSphere& rSphere)
{
	auto newEnd = std::remove(mCollidables.begin(), mCollidables.end(), &rSphere);
	mCollidables.erase(newEnd, mCollidables.end());
	rSphere.RemoveBox(*this);
}

void NuTo::SubBox::SetWalls(const std::list<CollidableWallBase*>& rCollidables)
{

	mWalls = rCollidables;

	for (CollidableWallBase* wall : mWalls)
		mCollidables.push_back(wall);
}

void NuTo::SubBox::AddWall(CollidableWallBase& rWall)
{
	mWalls.push_back(&rWall);
	mCollidables.push_back(&rWall);
}

const std::list<NuTo::CollidableWallBase*>& NuTo::SubBox::GetWalls() const
{
	return mWalls;
}

void NuTo::SubBox::Print()
{
	auto it = mCollidables.begin();
	// skip the walls
	for (unsigned int i = 0; i < mWalls.size(); ++i)
		it++;

	std::cout << mCollidables.size();

	for (auto i = it; i != mCollidables.end(); ++i)
		std::cout << (*i) << std::endl;
}

void NuTo::SubBox::CreateEvents(EventListHandler& rEvents,
		CollidableBase& rCollidable)
{

	const unsigned int size(mCollidables.size());

	// temporarily store event time and event type in vectors,
	// otherwise parallel event creation fails due to local list data race
	std::vector<double> newEvents(size);
	std::vector<int> eventType(size);

	// catch exceptions in parallel for
	bool parallelThrow = false;
	Exception parallelException("");

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) num_threads(mNumThreads)
#endif
	for (unsigned int i = 0; i < size; ++i)
	{
		try
		{
			newEvents[i] = rCollidable.PredictCollision(*mCollidables[i], eventType[i]);
		}catch(Exception& e)
		{
			parallelException = e;
			parallelThrow = true;
		}
	}

	if (parallelThrow)
		throw parallelException;

	for (unsigned int i = 0; i < size; ++i)
		if (newEvents[i] != Event::EVENTNULL)
			rEvents.AddEvent(newEvents[i], rCollidable, *mCollidables[i], eventType[i]);
}

const std::vector<NuTo::CollidableBase*>& NuTo::SubBox::GetCollidables() const
{
	return mCollidables;
}

void NuTo::SubBox::RemoveWall(CollidableWallBase& rWall)
{
	auto newEnd = std::remove(mCollidables.begin(), mCollidables.end(), &rWall);
	mCollidables.erase(newEnd, mCollidables.end());
	mWalls.remove(&rWall);
}

const int NuTo::SubBox::GetIndex() const
{
	return mIndex;
}

bool NuTo::SubBox::AddIfInside(CollidableParticleSphere& rSphere)
{
	bool isInside = true;
	// general solution:

	for (auto wall : mWalls)
		if (not wall->IsInside(rSphere))
		{
			isInside = false;
			break;
		}

	if (isInside)
	{
		rSphere.AddBox(*this);
		AddSphere(rSphere);
	}

	return isInside;
}
