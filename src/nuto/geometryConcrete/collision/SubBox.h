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

//! @brief ... class for sub box handing -> improves the performance without changing the physics
class SubBox
{
public:

	//! @brief ... constructor
	//! @param rIndex ... name for debug
	SubBox(const int rIndex);

	//! @brief ... destructor, Deletes mWalls on destruction
	~SubBox();

	//! @brief ... creates events between mCollidables and rCollidable
	//! @param rEvents ... global event list
	//! @param rCollidable ... collidable involved in this collision
	void CreateEvents(EventListHandler& rEvents, CollidableBase& rCollidable);

	//! @brief ... the sphere is now handled by this sub box
	void AddSphere(CollidableParticleSphere& rSphere);

	//! @brief ... adds a sphere only, if its inside of this box
	bool AddIfInside(CollidableParticleSphere& rSphere);

	//! @brief ... the sphere leaves this sub box
	void RemoveSphere(CollidableParticleSphere& rSphere);

	//! @brief ... sets the virtual and physical walls as sub box boundaries
	void SetWalls(const std::list<CollidableWallBase*>& rCollidables);

	//! @brief ... add a single wall
	void AddWall(CollidableWallBase& rWall);

	//! @brief ... remove a single wall
	void RemoveWall(CollidableWallBase& rWall);

	//! @brief ... print mCollidables without the walls
	void Print();

	//! @brief ... getter for mCollidables
	const std::vector<CollidableBase*>& GetCollidables() const;

	//! @brief ... getter for mWalls
	const std::list<CollidableWallBase*>& GetWalls() const;

	//! @brief ... getter for mIndex
	const int GetIndex() const;

private:
	const int mIndex;
	std::list<CollidableWallBase*> mWalls;
	std::vector<CollidableBase*> mCollidables;
};

} /* namespace NuTo */
#endif /* SUBBOX_H_ */
