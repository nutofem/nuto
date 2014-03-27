/*
 * CollidableWallCylinder.h
 *
 *  Created on: 13 Feb 2014
 *      Author: ttitsche
 */

#ifndef COLLIDABLEWALLCYLINDER_H_
#define COLLIDABLEWALLCYLINDER_H_

#include "nuto/geometryConcrete/collision/collidables/CollidableWallBase.h"

namespace NuTo
{
class CollidableParticleSphere;

//! @brief class for cylindric walls
class CollidableWallCylinder: public NuTo::CollidableWallBase
{
public:

	//! @brief ... constructor.
	//! @param rPosition ... center of the cylinder
	//! @param rDirection ... vector pointing to the top of the cylinder
	//! @param rRadius ... radius
	//! @param rHeigth ... height
	//! @param rIndex ... name
	CollidableWallCylinder(FullVector<double, 3> rPosition,
			FullVector<double, 3> rDirection, const double rRadius, const double rHeigth,
			const int rIndex);

	//! @brief ... collision between this and CollidableSphere
	//! @param rSphere ... collision partner
	void PerformCollision(CollidableParticleSphere& rSphere);

	//! @brief ... collision prediction between (this) and a sphere
	//! @param rSphere ... collision partner
	//! @return ... predicted time of collision
	const double PredictCollision(CollidableParticleSphere& rSphere, int& rType);

#ifdef ENABLE_VISUALIZE

	void VisualizationStatic(VisualizeUnstructuredGrid& rVisualizer) const override;
#endif

	bool IsInside(const CollidableParticleSphere& rSphere) const override;

protected:

	void Print(std::ostream & rReturnStream) const override;

	//! @brief ... cylinder radius
	const double mRadius;

	//! @brief ... cylinder height
	const double mHeigth;
};

} /* namespace NuTo */
#endif /* COLLIDABLEWALLCYLINDER_H_ */
