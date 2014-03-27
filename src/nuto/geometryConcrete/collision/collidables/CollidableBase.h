/*
 * CollidableBase.h
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#ifndef COLLIDABLEBASE_H_
#define COLLIDABLEBASE_H_

#include <iostream>
#include <list>

#include "nuto/math/FullVector.h"

#include "nuto/geometryConcrete/collision/Event.h"
#include "nuto/geometryConcrete/collision/SubBox.h"

namespace NuTo
{
class CollidableParticleSphere;
class CollidableWallBase;

//! @brief ... Base class for all collidables
class CollidableBase
{
	//! @brief ... friend class for easy local event list access
	friend class Event;
public:

	enum EventType {
		SphereCollision,
		WallCollision,
		WallTransfer
	};

	typedef std::vector<CollidableParticleSphere*> ParticleContainer;

	//! @brief ... constructor, initialized with an index
	CollidableBase(const int rIndex);

	//! @brief ... destructor
	virtual ~CollidableBase();

	//! @brief ... returnes index
	const int GetIndex() const;

	//! @brief ... collision handling, resolve double dispatch
	//! @param rCollidable/rSphere/rWall ... collision partner
	virtual void PerformCollision(CollidableBase& rCollidable) = 0;
	virtual void PerformCollision(CollidableParticleSphere& rSphere) = 0;
	virtual void PerformCollision(CollidableWallBase& rWall) = 0;

	//! @brief ... collision prediction between (this) and collision partner
	//! @param rCollidable/rSphere/rWall ... collision partner
	//! @return ... predicted time of collision
	virtual const double PredictCollision(CollidableBase& rCollidable, int& rType) = 0;
	virtual const double PredictCollision(CollidableParticleSphere& rSphere, int& rType) = 0;
	virtual const double PredictCollision(CollidableWallBase& rWall, int& rType) = 0;

	//! @brief ... adds a SubBox to this collidable
	//! @param rBox ... SubBox to add
	void AddBox(SubBox& rBox);

	//! @brief ... removes a SubBox from this collidable
	//! @param rBox ... SubBox to remove
	void RemoveBox(SubBox& rBox);

	//! @brief ... updates this collidable in time
	//! @param rTime ... new global time
	virtual void MoveAndGrow(const double rTime) = 0;


	//! @brief ... returns all old events, that need to be deleted
	//! @param rEventsToDelete ... return argument
	virtual void GetLocalEventsToDelete(Event::LocalEvents& rEventsToDelete) const = 0;

	//! @brief ... prints the local event list
	void PrintLocalEvents() const;

	//! @brief ... returns all the SubBoxes of this collidable
	const std::list<SubBox*>& GetSubBoxes() const;

protected:

	//! @brief ... index, just a name for each collidable, multiple collidables with the same index possible
	const int mIndex;

	//! @brief ... list of SubBoxes in which this collidable is inside
	//! if the collidable is passing a SubBox wall, multiple mBoxes are possible
	std::list<SubBox*> mBoxes;

	//! @brief ... local event list
	//! adding single event: through EventBase::AddLocalEvents()
	//! remove single event: through EventBase::~EventBase()
	//! clear everything: on collision
	Event::LocalEvents mLocalEvents;

private:

	//! @brief ... standard output for all collidables, calls Print(...) for polymorph behaviour
	//! @param rOutStream ... operator RHS: output stream
	//! @param rCollidable ... operator LHS: object added to the output stream
	//! @return ... modified output stream
#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& rOutStream, const CollidableBase* rCollidable);
#endif
	//! @brief ... prints information based on derived class type.
	//! @param rReturnStream ... output stream, that gets modified
	virtual void Print(std::ostream & rReturnStream) const = 0;

};

} /* namespace NuTo */
#endif /* COLLIDABLEBASE_H_ */
