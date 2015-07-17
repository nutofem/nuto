/*
 * Interpolation1DTrussIn2D.h
 *
 *  Created on: 10 July 2015
 *      Author: phuschke
 */

#ifndef INTERPOLATION1DTRUSSIN2D_H_
#define INTERPOLATION1DTRUSSIN2D_H_

#include "nuto/mechanics/interpolationtypes/Interpolation1DTruss.h"

namespace NuTo
{
class Interpolation1DTrussIn2D: public Interpolation1DTruss
{
public:

    Interpolation1DTrussIn2D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief return the number of dofs per node depending on dimension
    int GetNumDofsPerNode() const override;
};

} /* namespace NuTo */

#endif /* INTERPOLATION1DTRUSSIN2D_H_ */
