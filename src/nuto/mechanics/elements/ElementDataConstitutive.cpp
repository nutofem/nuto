// $Id: $

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataConstitutive.h"
#include <assert.h>



NuTo::ElementDataConstitutive::ElementDataConstitutive(const NuTo::IntegrationTypeBase* rIntegrationType)
: NuTo::ElementDataConstitutiveBase::ElementDataConstitutiveBase(rIntegrationType), NuTo::ElementDataBase::ElementDataBase(rIntegrationType)
{
	std::cout << "ElementDataConstitutive constructor " << std::endl;
}

//! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
void NuTo::ElementDataConstitutive::UpdateForModifiedIntegrationType(const ElementBase* rElement)
{
    //nothing to be done
}

//! @brief update the information related to a modification of the constitutive law, e.g. reallocation of the static data
void NuTo::ElementDataConstitutive::UpdateForModifiedConstitutiveLaw(const ElementBase* rElement)
{
    //nothing to be done
}
