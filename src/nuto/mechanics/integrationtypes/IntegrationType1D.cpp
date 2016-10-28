#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D.h"

bool NuTo::IntegrationType1D::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENT:
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
