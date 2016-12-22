/*
 * InterpolationType3D.cpp
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "mechanics/MechanicsException.h"
#include "mechanics/interpolationtypes/Interpolation3D.h"
#include "mechanics/nodes/NodeEnum.h"

NuTo::Interpolation3D::Interpolation3D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        InterpolationBaseFEM::InterpolationBaseFEM(rDofType, rTypeOrder, rDimension)
{

}

int NuTo::Interpolation3D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::eDof::COORDINATES:
        return 3;
    case NuTo::Node::eDof::DISPLACEMENTS:
        return 3;
    case NuTo::Node::eDof::TEMPERATURE:
        return 1;
    case NuTo::Node::eDof::NONLOCALEQSTRAIN:
        return 1;
    case NuTo::Node::eDof::RELATIVEHUMIDITY:
        return 1;
    case NuTo::Node::eDof::WATERVOLUMEFRACTION:
        return 1;
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "dof type not found.");
    }
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation3D)
#endif
