/*
 * CollidableBase.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */


#include "geometryConcrete/collision/Event.h"
#include "geometryConcrete/collision/collidables/CollidableBase.h"
#include "geometryConcrete/collision/collidables/CollidableWallBase.h"
#include "geometryConcrete/collision/collidables/CollidableParticleSphere.h"
#include <algorithm>


NuTo::CollidableBase::CollidableBase(const int rIndex)
		: mIndex(rIndex)
{
}

NuTo::CollidableBase::~CollidableBase()
{
}

const int NuTo::CollidableBase::GetIndex() const
{
	return mIndex;
}

void NuTo::CollidableBase::PrintLocalEvents() const
{
	for (auto it = mLocalEvents.begin(); it != mLocalEvents.end(); ++it) {
		std::cout << *(*it) << std::endl;
	}
}

const std::vector<NuTo::SubBox*>& NuTo::CollidableBase::GetSubBoxes() const
{
	return mBoxes;
}



void NuTo::CollidableBase::AddBox(SubBox& rBox)
{
	mBoxes.push_back(&rBox);
}

void NuTo::CollidableBase::RemoveBox(SubBox& rBox)
{
	auto newEnd = std::remove(mBoxes.begin(), mBoxes.end(), &rBox);
	mBoxes.erase(newEnd, mBoxes.end());
	if (mBoxes.size() == 0)
		throw NuTo::Exception("[NuTo::CollidableBase::AddBox] Collidable is not handled by any box.");
}


namespace NuTo{


std::ostream& operator <<(std::ostream& rOutStream,
		const CollidableBase* rCollidable)
{
	rCollidable->Print(rOutStream);
	return rOutStream;
}
}
