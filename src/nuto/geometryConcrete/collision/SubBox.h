/*
 * CellGridBox.h
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#ifndef SUBBOX_H_
#define SUBBOX_H_

#include <list>
#include "nuto/geometryConcrete/collision/Event.h"

namespace NuTo
{
class EventListHandler;
class CollidableParticleSphere;
class CollidableBase;
class CollidableWallBase;

class SubBox
{
public:
	SubBox(const int rIndex);
	~SubBox();

	void AddSphere(CollidableParticleSphere& rSphere);
	void RemoveSphere(CollidableParticleSphere& rSphere);
	void SetWalls(const std::list<CollidableWallBase*>& rCollidables);
	void AddWall(CollidableWallBase& rWall);
	void RemoveWall(CollidableWallBase& rWall);
	const std::list<CollidableWallBase*>& GetWalls() const;

	double CreateEvents(EventListHandler& rEvents, CollidableBase& rCollidable);

	void Print();
	const std::list<CollidableBase*>& GetCollidables() const;

	bool AddIfInside(CollidableParticleSphere& rSphere);

	const int GetIndex() const;



private:
	const int mIndex;
	std::list<CollidableWallBase*> mWalls;
	std::list<CollidableBase*> mCollidables;
};

} /* namespace NuTo */
#endif /* SUBBOX_H_ */
