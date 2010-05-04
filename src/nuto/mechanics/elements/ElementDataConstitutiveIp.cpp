// $Id: $
#include "nuto/mechanics/elements/ElementDataConstitutiveIp.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include <assert.h>


NuTo::ElementDataConstitutiveIp::ElementDataConstitutiveIp(const ElementWithDataBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType,	NuTo::IpData::eIpDataType rIpDataType) :
		NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase() , ElementDataIpBase(rElement,rIntegrationType,rIpDataType)
{
	std::cout << "NuTo::ElementDataConstitutiveIp::ElementDataConstitutiveIp" << std::endl;
}

NuTo::ElementDataConstitutiveIp::~ElementDataConstitutiveIp()
{
	//std::cout << "NuTo::ElementDataConstitutiveStaticData::~ElementDataConstitutiveStaticData()" << std::endl;
}

void NuTo::ElementDataConstitutiveIp::SetConstitutiveLaw(const ElementWithDataBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	std::cout << "NuTo::ElementDataConstitutiveIp::SetConstitutiveLaw" << std::endl;
	NuTo::ElementDataConstitutiveBase::SetConstitutiveLaw(rConstitutiveLaw);

	//reinitialize ip data (f.e. if different static data or nonlocal data are required with another constitutive model)
	for (int theIp=0; theIp<(int)mIpData.size();theIp++)
	{
		mIpData[theIp].Initialize(rElement,theIp);
	}
}
