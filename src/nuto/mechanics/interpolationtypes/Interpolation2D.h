/*
 * Interpolation2D.h
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION2D_H_
#define INTERPOLATION2D_H_

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation2D: public InterpolationBase
{
public:
    Interpolation2D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    const std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    //! @brief return the number of dofs per node depending on dimension
    int GetNumDofsPerNode() const override;

};

} /* namespace NuTo */

#endif /*INTERPOLATION2D_H_ */
