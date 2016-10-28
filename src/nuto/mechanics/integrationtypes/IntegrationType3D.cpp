#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D.h"

bool NuTo::IntegrationType3D::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType3D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationType3D)
#endif
