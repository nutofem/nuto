/*
 * CollidableParticleBase.cpp
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#include "geometryConcrete/collision/collidables/CollidableParticleBase.h"
#include "base/Exception.h"

NuTo::CollidableParticleBase::CollidableParticleBase(Eigen::Vector3d rPosition, Eigen::Vector3d rVelocity, int rIndex)
    : CollidableBase(rIndex)
    , mPosition(rPosition)
    , mVelocity(rVelocity)
{
}
