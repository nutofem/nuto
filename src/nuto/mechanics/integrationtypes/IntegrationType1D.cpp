// $Id$ 
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
    case NuTo::Element::ELEMENT1D:
        return true;
    case NuTo::Element::BOUNDARYELEMENT2D:
        return true;
    case NuTo::Element::BOUNDARYELEMENT2DADDITIONALNODE:
        return true;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::ELEMENT1DINXD:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return true;
    case NuTo::Element::TRUSS1D4NDISP3NX:
        return true;
    case NuTo::Element::TRUSS1D2NSPECTRALORDER2:
        return true;
    case NuTo::Element::TRUSS1D2NSPECTRALORDER3:
        return true;
    case NuTo::Element::TRUSS1D2NSPECTRALORDER4:
        return true;
    default:
        return false;
    }
}

#ifdef ENABLE_SERIALIZATION
template<class Archive>
void NuTo::IntegrationType1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize IntegrationType0DBoundary" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuTo::IntegrationTypeBase);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize IntegrationType0DBoundary" << std::endl;
#endif
}
#endif
