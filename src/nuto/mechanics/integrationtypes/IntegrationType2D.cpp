// $ld: $ 
// IntegrationType2D.cpp
// created Apr 30, 2010 by Joerg F. Unger

#include "nuto/mechanics/integrationtypes/IntegrationType2D.h"

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::IntegrationType2D::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    default:
        return false;
    }
}
