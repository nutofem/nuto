/*
 * CollidableWallBase.h
 *
 *  Created on: 17 Jan 2014
 *      Author: ttitsche
 */

#pragma once

#include "geometryConcrete/collision/collidables/CollidableBase.h"
#include "math/FullVector_Def.h"

namespace NuTo
{
class SubBox;
class CollidableParticleSphere;
class VisualizeUnstructuredGrid;

//! @brief ... base class for walls
class CollidableWallBase: public NuTo::CollidableBase
{
public:
	//! @brief ... constructor, create CollidableWallBase using the point-and-normal-vector plane definition
	//! @param rPosition ... point on the plane
	//! @param rDirection ... normal vector pointing inside the domain, gets normalized.
	//! @param rIndex ... name
	CollidableWallBase(
			NuTo::FullVector<double, Eigen::Dynamic> rPosition,
			NuTo::FullVector<double, Eigen::Dynamic> rDirection,
			const int rIndex);

	//! @brief ... destructor, do nothing
	virtual ~CollidableWallBase();

	//! @brief ... every wall should know its box and the neighbour box
	//! @param rInsideBox ... box on the inner side of the wall
	//! @param rOutsideBox ... box on the outer side of the wall
	void SetBoxes(SubBox& rInsideBox,SubBox& rOutsideBox);

	//! @brief ... collision between CollidableWall and CollidableBase, resolve double dispatch, forward *this
	//! @param rCollidable ... collision partner
	void PerformCollision(CollidableBase& rCollidable);

	//! @brief ... collision between CollidableWall and CollidableSphere, pass to child classes
	//! @param rSphere ... collision partner
	virtual void PerformCollision(CollidableParticleSphere& rSphere) = 0;

	//! @brief ... collision between CollidableWall and CollidableWall, does nothing
	//! @param rWall ... collision partner
	void PerformCollision(CollidableWallBase& rWall);

	//! @brief ... collision check between CollidableWall and CollidableBase, resolve double dispatch, forward *this
	//! @param rCollidable ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	const double PredictCollision(CollidableBase& rCollidable, int& rType);

	//! @brief ... collision check between CollidableWall and CollidableSphere, pass to child classes
	//! @param rSphere ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	virtual const double PredictCollision(CollidableParticleSphere& rSphere, int& rType) = 0;

	//! @brief ... collision check between CollidableWall and CollidableWall, does nothing
	//! @param rWall ... possible collision partner
	//! @param rType ... return argument, element of enum CollidableBase::EventType
	//! @return ... predicted collision time
	const double PredictCollision(CollidableWallBase& rWall, int& rType);

	//! @brief ... walls to neither grow nor move, do nothing
	void MoveAndGrow(const double rTime);

	//! @brief ... ture for physical walls
	virtual const bool IsPhysical() const = 0;

#ifdef ENABLE_VISUALIZE
	//! @brief ... visualize all non-moving collidables
	//! @param rVisualizer ... NuTo object for ascii-export
	virtual void VisualizationStatic(VisualizeUnstructuredGrid& rVisualizer) const;
#endif

	//! @brief ... does nothing as all other wall events are still legal
    void GetLocalEventsToDelete(LocalEvents& rEventsToDelete) const;

	//! @brief ... returns whether a sphere is on the positive side of this wall
	//! @param rSphere ... sphere to test
	virtual bool IsInside(const CollidableParticleSphere& rSphere) const;

	//! @brief getter for the wall direction
	const NuTo::FullVector<double, Eigen::Dynamic> GetDirection() const;

	//! @brief getter fot the wall position
	const NuTo::FullVector<double, Eigen::Dynamic> GetPosition() const;

protected:

	//! @brief ... point of the point-and-normal-vector plane definition
	FullVector<double, 3> mPosition;

	//! @brief ... normal vector of the point-and-normal-vector plane definition
	FullVector<double, 3> mDirection;

	//! @brief ... box on the inner side of the wall
	SubBox* mInsideBox;

	//! @brief ... box on the outer side of the wall
	SubBox* mOutsideBox;

	//! @brief ... prints CollidableWall
	//! @param rReturnStream ... output stream, that gets modified
	virtual void Print(std::ostream & rReturnStream) const;


	//! @brief ... index of the direction component that is != 0
	const int mNonNullAxis;

	//! @brief ... true --> calculations without vector calculations
	const bool mIsAxisAligned;

private:

	//! @brief ... calculates index of the direction component that is != 0
	const int GetNonNullAxis();
};

} /* namespace NuTo */
