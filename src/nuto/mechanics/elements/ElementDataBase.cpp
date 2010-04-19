#include <assert.h>
#include "nuto/mechanics/elements/ElementDataBase.h"

    //! @brief constructor
NuTo::ElementDataBase::ElementDataBase(const NuTo::IntegrationTypeBase* rIntegrationType)
{
    //std::cout << "ElementDataBase constructor " << std::endl;
	mIntegrationType = rIntegrationType;
}

//! @brief deconstructor
NuTo::ElementDataBase::~ElementDataBase()
{
}

    //! @brief sets the constitutive law for all integration points of the element
//! @param rConstitutiveLaw constitutive law
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass.");
}

//! @brief sets the constitutive law for a single integration point of the element
//! @param rConstitutiveLaw constitutive law
//! @param rIp integration point
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass.");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the ElementDataClass.");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the ElementDataClass.");
}

//! @brief returns the constitutive law of an integration point
//! @param rIp integration point
//! @return constitutive law
NuTo::ConstitutiveBase* NuTo::ElementDataBase::GetConstitutiveLaw(int rIp)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetConstitutiveLaw] Not implemented for the ElementDataClass.");
}

//! @brief returns the constitutive law of an integration point
//! @param rIp integration point
//! @return constitutive law
const NuTo::ConstitutiveBase* NuTo::ElementDataBase::GetConstitutiveLaw(int rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetConstitutiveLaw] Not implemented for the ElementDataClass.");
}

//! @brief sets the integration type of an element
//! @param rIntegrationType pointer to integration type
void NuTo::ElementDataBase::SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType)
{
	mIntegrationType = rIntegrationType;
	UpdateForModifiedIntegrationType(rElement);
}

//! @brief returns a pointer to the integration type of an element
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataBase::GetIntegrationType()const
{
	assert(mIntegrationType!=0);
	return mIntegrationType;
}

