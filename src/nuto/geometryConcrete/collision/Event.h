/*
 * EventBase.h
 *
 *  Created on: 21 Jan 2014
 *      Author: ttitsche
 */

#ifndef EVENT_H_
#define EVENT_H_

#include <iostream>
#include <list>
#include <boost/ptr_container/ptr_set.hpp>
#include <boost/ptr_container/ptr_unordered_set.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <nuto/base/Exception.h>
#include <iostream>

namespace NuTo
{
class CollidableBase;
class EventListHandler;

class Event
{
public:

	typedef boost::ptr_set<Event> GlobalEvents;
//	typedef boost::ptr_unordered_set<Event> GlobalEvents;

	typedef std::list<Event*> LocalEvents;

	static const double EVENTNULL;

	bool operator<(const Event& rOther) const;
	bool operator==(Event const& rRhs) const;
	bool operator!=(Event const& rRhs) const;


	//! @brief ... constructor
	//! initialized with the two CollidableBase objects  involved in this collision
	//! adds itself (this) to the local event lists of both collidables
	//! @param rTime ... event time
	//! @param rFirst ... first CollidableBase involved
	//! @param rSecond ... second CollidableBase involved
	Event(const double rTime, CollidableBase* rFirst,
			CollidableBase* rSecond, const int rType);

	Event(const Event& rEvent);

	//! @brief ... destructor
	//! removes itself (this) from the local event lists of both collidables
	virtual ~Event();

	//! @brief ... getter for mTime
	//! @return ... time of event
	const double GetTime() const;


	//! @brief ... creates new events for mFirst and mSecond
	//! (--> automatically added to local event lists, see constructor)
	//! stores them to the global event list
	//! @param rEvents ... global event list
	//! @param rCollidables ... global collidable list
	double AddNewEventsLocalAndGlobal(EventListHandler& rEvents) const;


	void EraseOldEventsLocalAndGlobal(EventListHandler& rEvents) const;

	//! @brief ... performs the collision of mFirst vs. mSecond
	void PerformCollision() const;

	//! @brief adds itself (this) to the local event lists of both collidables
	void AddLocalEvent();
	const int GetType() const;


	std::size_t GetHash() const;

#ifndef SWIG
	friend std::size_t hash_value(const Event& rEvent);
#endif


protected:

	//! @brief ... first CollidableBase involved
	CollidableBase* mFirst;

	//! @brief ... second CollidableBase involved
	CollidableBase* mSecond;

private:



	//! @brief ... event time
	const double mTime;
	const int mType;

#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& rOutStream,
			const Event& rEvent);
#endif

	virtual void Print(std::ostream& rOutStream) const;

	//! @brief removes itself (this) from the local event lists of both collidables
	void RemoveLocalEvent();
};

} /* namespace NuTo */


#endif /* EVENT_H_ */
