/*
 * CellGridBox.h
 *
 *  Created on: 5 Feb 2014
 *      Author: ttitsche
 */

#pragma once

#include <vector>

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
    SubBox(int rIndex);

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
    void SetWalls(const std::vector<CollidableWallBase*>& rCollidables);

    //! @brief ... add a single wall
    void AddWall(CollidableWallBase& rWall);

    //! @brief ... remove a single wall
    void RemoveWall(CollidableWallBase& rWall);

    //! @brief ... print mCollidables without the walls
    void Print();

    //! @brief ... getter for mCollidables
    const std::vector<CollidableBase*>& GetCollidables() const;

    //! @brief ... getter for mWalls
    const std::vector<CollidableWallBase*>& GetWalls() const;

    //! @brief ... getter for mIndex
    int GetIndex() const;

private:
    int mIndex;
    std::vector<CollidableWallBase*> mWalls;
    std::vector<CollidableBase*> mCollidables;
};

} /* namespace NuTo */
