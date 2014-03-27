/*
 * CollidableParticleBase.h
 *
 *  Created on: 10 Mar 2014
 *      Author: ttitsche
 */

#ifndef COLLIDABLEPARTICLEBASE_H_
#define COLLIDABLEPARTICLEBASE_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableBase.h"

namespace NuTo
{

class CollidableParticleBase: public NuTo::CollidableBase
{
public:
	CollidableParticleBase(
			FullVector<double,Eigen::Dynamic> rPosition,
			FullVector<double,Eigen::Dynamic> rVelocity,
			const int rIndex);

	virtual const double GetVolume() const = 0;
	virtual const double GetKineticEnergy() const = 0;
	virtual void MoveAndGrow(const double rTime) = 0;

	virtual void Print(std::ostream & rReturnStream) const;


protected:
	FullVector<double, 3> mPosition;
	FullVector<double, 3> mVelocity;

};

} /* namespace NuTo */
#endif /* COLLIDABLEPARTICLEBASE_H_ */
