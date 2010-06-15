// $Id: $
#include "nuto/mechanics/elements/ElementDataConstitutiveIp.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveIp::ElementDataConstitutiveIp(const ElementBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType,	NuTo::IpData::eIpDataType rIpDataType) :
		NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase() , ElementDataIpBase(rElement,rIntegrationType,rIpDataType)
{
	//std::cout << "NuTo::ElementDataConstitutiveIp::ElementDataConstitutiveIp" << std::endl;
}

NuTo::ElementDataConstitutiveIp::~ElementDataConstitutiveIp()
{
	//std::cout << "NuTo::ElementDataConstitutiveStaticData::~ElementDataConstitutiveStaticData()" << std::endl;
}

//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
//! @param rElement element
void NuTo::ElementDataConstitutiveIp::InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)
{
	//reinitialize ip data (f.e. if different static data or nonlocal data are required with another constitutive model)
	for (int theIp=0; theIp<(int)mIpData.size();theIp++)
	{
		mIpData[theIp].Initialize(rElement,mConstitutiveLaw);
	}

}
