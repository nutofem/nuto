/*
 * CellGridBox.cpp
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#include <iostream>
#include "base/Exception.h"

#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/SubBox.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include "geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "geometryConcrete/collision/handler/EventListHandler.h"

NuTo::SubBox::SubBox(int rIndex)
		: mIndex(rIndex)
{
}

NuTo::SubBox::~SubBox()
{
    for (auto& wall : mWalls)
        delete wall;
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

void NuTo::SubBox::SetWalls(const std::vector<CollidableWallBase*>& rCollidables)
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

const std::vector<NuTo::CollidableWallBase*>& NuTo::SubBox::GetWalls() const
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
	for (auto* collidable : mCollidables)
	{
	    int eventType;
		double collisionTime = rCollidable.PredictCollision(*collidable, eventType);
        if (collisionTime != Event::EVENTNULL)
            rEvents.AddEvent(collisionTime, rCollidable, *collidable, eventType);
	}
}

const std::vector<NuTo::CollidableBase*>& NuTo::SubBox::GetCollidables() const
{
	return mCollidables;
}

void NuTo::SubBox::RemoveWall(CollidableWallBase& rWall)
{
	auto newEnd = std::remove(mCollidables.begin(), mCollidables.end(), &rWall);
	mCollidables.erase(newEnd, mCollidables.end());
    mWalls.erase(std::remove(mWalls.begin(), mWalls.end(), &rWall), mWalls.end());
}

int NuTo::SubBox::GetIndex() const
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
