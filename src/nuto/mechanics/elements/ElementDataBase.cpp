#include <assert.h>
#include "nuto/mechanics/elements/ElementDataBase.h"

//! @brief constructor
NuTo::ElementDataBase::ElementDataBase()
{
    //std::cout << "ElementDataBase constructor " << std::endl;
}

//! @brief deconstructor
NuTo::ElementDataBase::~ElementDataBase()
{
    //std::cout << "NuTo::ElementDataBase::~ElementDataBase()" << std::endl;
}

    //! @brief sets the constitutive law for all integration points of the element
//! @param rConstitutiveLaw constitutive law
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementWithDataBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief sets the constitutive law for a single integration point of the element
//! @param rConstitutiveLaw constitutive law
//! @param rIp integration point
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementWithDataBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the Ele - check the allocated element data type.mentDataClass.");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the constitutive law of an integration point
//! @param rIp integration point
//! @return constitutive law
NuTo::ConstitutiveBase* NuTo::ElementDataBase::GetConstitutiveLaw(int rIp)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the constitutive law of an integration point
//! @param rIp integration point
//! @return constitutive law
const NuTo::ConstitutiveBase* NuTo::ElementDataBase::GetConstitutiveLaw(int rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief sets the integration type of an element
//! @param rIntegrationType pointer to integration type
void NuTo::ElementDataBase::SetIntegrationType(const ElementWithDataBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetIntegrationType] Not implemented for the ElementDataClass - check the allocated element data type.");
}

//! @brief returns a pointer to the integration type of an element
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataBase::GetIntegrationType()const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetIntegrationType] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief adds the nonlocal weight to an integration point
//! @param rLocalIpNumber local Ip
//! @param rConstitutive constitutive model for which nonlocal data is to be calculated
//! @param rNonlocalElement element of the nonlocal ip
//! @param rNonlocalIp local ip number of the nonlocal ip
//! @param rWeight weight
 void NuTo::ElementDataBase::SetNonlocalWeight(int rLocalIpNumber, const ConstitutiveBase* rConstitutive,
		const ElementWithDataBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
    throw MechanicsException("[NuTo::ElementDataBase::AddNonlocalIp] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

//! @brief gets the nonlocal elements for a constitutive model
//! @param rConstitutive constitutive model
//! @return vector to nonlocal elements
const std::vector<const NuTo::ElementWithDataBase*>& NuTo::ElementDataBase::GetNonlocalElements(const ConstitutiveBase* rConstitutive)const
{
    throw MechanicsException("[NuTo::ElementDataBase::GetNonlocalElements] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

//! @brief gets the nonlocal weights
//! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
//! @return vector of weights for all integration points of the nonlocal element
const std::vector<double>& NuTo::ElementDataBase::GetNonlocalWeights(int rIp, int rNonlocalElement, const ConstitutiveBase* rConstitutive)const
{
    throw MechanicsException("[NuTo::ElementDataBase::GetNonlocalWeights] Not implemented for the ElementDataBase class - check the allocated element data type..");
}
