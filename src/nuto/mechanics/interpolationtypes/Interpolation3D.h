/*
 * Interpolation3D.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION3D_H_
#define INTERPOLATION3D_H_

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation3D: public InterpolationBase
{
public:
    Interpolation3D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const override;

};

} /* namespace NuTo */

#endif /*INTERPOLATION3D_H_ */
