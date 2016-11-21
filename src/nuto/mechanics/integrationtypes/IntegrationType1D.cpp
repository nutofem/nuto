// $Id$ 
// IntegrationType1D.cpp
// created Apr 30, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D.h"

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::IntegrationType1D::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENT:
    case NuTo::Element::eElementType::CONTINUUMCONTACTELEMENT:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
    case NuTo::Element::eElementType::ELEMENT1DINXD:
        return true;
    default:
        return false;
    }
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationType1D)
#endif
