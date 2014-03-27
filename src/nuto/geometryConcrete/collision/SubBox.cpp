/*
 * CellGridBox.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/WallTime.h"
#include "nuto/geometryConcrete/collision/SubBox.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"
#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "nuto/geometryConcrete/collision/handler/EventListHandler.h"
#include <omp.h>

NuTo::SubBox::SubBox(const int rIndex)
		: mIndex(rIndex)
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
	mCollidables.remove(&rSphere);
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

double NuTo::SubBox::CreateEvents(EventListHandler& rEvents,
		CollidableBase& rCollidable)
{
	// convert
	// list ( fast add / remove ) to
	// vector ( fast iteration for parallelization )
	// conversion overhead approx. 0
	const int size(mCollidables.size());
	std::vector<CollidableBase*> collidables(size);
	int vecI = 0;

	for (auto coll : mCollidables)
		collidables[vecI++] = coll;

	// temporarily store event time and event type in vectors,
	// because parallel event creation fails due to local list data race
	std::vector<double> newEvents(size);
	std::vector<int> eventType(size);

	// measure parallel part
	double timeTmp = WallTime::Get();

	// catch exceptions in parallel for
	bool parallelThrow = false;
	Exception parallelException("");

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(shared) num_threads(4)
#endif
	for (int i = 0; i < size; ++i)
	{
		try
		{
			newEvents[i] = rCollidable.PredictCollision(*collidables[i], eventType[i]);
		}catch(Exception& e)
		{
			parallelException = e;
			parallelThrow = true;
		}
	}

	if (parallelThrow)
		throw parallelException;

	timeTmp = WallTime::Get() - timeTmp;

	for (int i = 0; i < size; ++i)
		if (newEvents[i] != Event::EVENTNULL)
			rEvents.AddEvent(newEvents[i], rCollidable, *collidables[i], eventType[i]);

	return timeTmp;
}

const std::list<NuTo::CollidableBase*>& NuTo::SubBox::GetCollidables() const
{
	return mCollidables;
}

void NuTo::SubBox::RemoveWall(CollidableWallBase& rWall)
{
	mCollidables.remove(&rWall);
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
