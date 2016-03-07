/*
 * Interpolation1D.h
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION1D_H_
#define INTERPOLATION1D_H_

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation1D: public InterpolationBase
{
public:
    Interpolation1D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    const std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const override;

};

} /* namespace NuTo */

#endif /* INTERPOLATION1D_H_ */
