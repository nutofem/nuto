// $Id: $

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataConstitutive.h"
#include <assert.h>


NuTo::ElementDataConstitutive::ElementDataConstitutive(const NuTo::IntegrationTypeBase* rIntegrationType) : ElementDataBase(rIntegrationType)
{
	mConstitutiveLaw = 0;
}

void NuTo::ElementDataConstitutive::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	mConstitutiveLaw = rConstitutiveLaw;
}

void NuTo::ElementDataConstitutive::SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataConstitutive::SetConstitutiveLaw] Material can only be set for all integration points.");
}

NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataConstitutive::GetStaticData(int rIp)
{
	// no static data
	return 0;
}

const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataConstitutive::GetStaticData(int rIp)const
{
	return 0;
}

NuTo::ConstitutiveBase* NuTo::ElementDataConstitutive::GetConstitutiveLaw(int rIp)
{
	if (mConstitutiveLaw==0)
		throw MechanicsException("[NuTo::ElementDataConstitutive::GetConstitutiveLaw] No material set for element.");
	return mConstitutiveLaw;
}

const NuTo::ConstitutiveBase* NuTo::ElementDataConstitutive::GetConstitutiveLaw(int rIp)const
{
	if (mConstitutiveLaw==0)
		throw MechanicsException("[NuTo::ElementDataConstitutive::GetConstitutiveLaw] No material set for element.");
	return mConstitutiveLaw;
}

//! @brief update the information related to a modification of the integration type, e.g. reallocation of the static data
void NuTo::ElementDataConstitutive::UpdateForModifiedIntegrationType(const ElementBase* rElement)
{
    //nothing to be done
}
