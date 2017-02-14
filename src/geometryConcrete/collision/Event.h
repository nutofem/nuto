/*
 * EventBase.h
 *
 *  Created on: 21 Jan 2014
 *      Author: ttitsche
 */

#pragma once


#include <boost/ptr_container/ptr_set.hpp>
#include <vector>


namespace NuTo
{
class CollidableBase;
class EventListHandler;


//! @brief ... class for storing events
class Event
{
public:

	typedef std::vector<Event*> LocalEvents;

	//! @brief ... statistics
	enum EventType {
		SphereCollision,
		WallCollision,
		WallTransfer
	};

	//! @brief ... identifier for null events
	static const double EVENTNULL;

	//! @brief ... constructor, initialized with the two CollidableBase objects  involved in this collision
	//! @param rTime ... event time
	//! @param rFirst ... first CollidableBase involved
	//! @param rSecond ... second CollidableBase involved
	Event(const double rTime, CollidableBase* rFirst,
			CollidableBase* rSecond, const int rType);

    //! @brief ... copy constructor
    Event(const Event& rEvent) = default;

    //! @brief ... move constructor
    Event(Event&& rEvent) = default;

    Event& operator=(const NuTo::Event&) = default;

    Event& operator=(Event&&) = default;


	//! @brief ... important operator for the event list sorting
	//! ... sort priority: time >> collidables
	bool operator<(const Event& rOther) const;

	//! @brief ... determines, whether two events are equal
	bool operator==(Event const& rRhs) const;

	//! @brief ... determines, whether two events are unequal
	bool operator!=(Event const& rRhs) const;



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
	void AddNewEvents(EventListHandler& rEvents) const;

	//! @brief ... removes all events in the local event lists of rEvent
	//! @param rEvents ... global event list
	void EraseOldEvents(EventListHandler& rEvents) const;

	//! @brief ... performs the collision of mFirst vs. mSecond
	void PerformCollision() const;

	//! @brief adds itself (this) to the local event lists of both collidables
	void AddLocalEvent();

	//! @brief ... getter for type
	const int GetType() const;

protected:

	//! @brief ... first CollidableBase involved
	CollidableBase* mFirst;

	//! @brief ... second CollidableBase involved
	CollidableBase* mSecond;

private:

	//! @brief ... event time
	double mTime;

	//! @brief ... member of CollidableBase::EventType
	int mType;

#ifndef SWIG
	friend std::ostream& operator<<(std::ostream& rOutStream,
			const Event& rEvent);
#endif

	//! @brief ... output
	//! @param rOutStream ... return argument, gets modified
	virtual void Print(std::ostream& rOutStream) const;

};

} /* namespace NuTo */


