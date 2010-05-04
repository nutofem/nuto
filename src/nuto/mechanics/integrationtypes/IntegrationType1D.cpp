// $ld: $ 
// IntegrationType1D.cpp
// created Apr 30, 2010 by Joerg F. Unger

#include "nuto/mechanics/integrationtypes/IntegrationType1D.h"
#include "nuto/mechanics/elements/ElementBase.h"

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::IntegrationType1D::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    default:
        return false;
    }
}
