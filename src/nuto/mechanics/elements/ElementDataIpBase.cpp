// $ld: $ 
// ElementDataIpDataBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/mechanics/elements/IpDataEmpty.h"
#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

//! @brief constructor
NuTo::ElementDataIpBase::ElementDataIpBase(const ElementWithDataBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType) : NuTo::ElementDataBase()
{
    std::cout << "ElementDataIpBase constructor " << std::endl;
	mIntegrationType = rIntegrationType;
	mIpData.clear();
	mIpData.reserve(mIntegrationType->GetNumIntegrationPoints());
	printf("num int ppoints %d",mIntegrationType->GetNumIntegrationPoints());
	for (int theIp=0; theIp<mIntegrationType->GetNumIntegrationPoints();theIp++)
	{
		switch (rIpDataType)
		{
		case NuTo::IpData::NOIPDATA:
			printf("[NuTo::ElementDataIpBase::ElementDataIpBase]empty\n");
			mIpData.push_back(new IpDataEmpty());
			break;
		case NuTo::IpData::STATICDATA:
			printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATA\n");
			mIpData.push_back(new IpDataStaticData());
			break;
		case NuTo::IpData::STATICDATANONLOCAL:
			printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATANONLOCAL\n");
			mIpData.push_back(new IpDataStaticDataNonlocal());
			break;
		default:
			throw MechanicsException("[NuTo::ElementDataIpBase::ElementDataIpBase] Ip data type not known.");
		}
		//initialize data without constitutive law - store zeros as pointers e.g. to static data reallocation, when constitutive model is modified
		mIpData[theIp].Initialize(rElement, (ConstitutiveBase*)0);
	}
}

//! @brief deconstructor
NuTo::ElementDataIpBase::~ElementDataIpBase()
{
    //std::cout << "NuTo::ElementDataBase::~ElementDataBase()" << std::endl;
}

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rElement pointer to element
//! @param rIntegrationType pointer to integration type
void NuTo::ElementDataIpBase::SetIntegrationType(const ElementWithDataBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
	mIntegrationType = rIntegrationType;
	mIpData.clear();
	mIpData.reserve(mIntegrationType->GetNumIntegrationPoints());
	for (int theIp=0; theIp<mIntegrationType->GetNumIntegrationPoints();theIp++)
	{
		switch (rIpDataType)
		{
		case NuTo::IpData::NOIPDATA:
			printf("[NuTo::ElementDataIpBase::SetIntegrationType]empty\n");
			mIpData.push_back(new IpDataEmpty());
			break;
		case NuTo::IpData::STATICDATA:
			printf("[NuTo::ElementDataIpBase::SetIntegrationType]STATICDATA\n");
			mIpData.push_back(new IpDataStaticData());
			break;
		case NuTo::IpData::STATICDATANONLOCAL:
			printf("[NuTo::ElementDataIpBase::SetIntegrationType]STATICDATANONLOCAL\n");
			mIpData.push_back(new IpDataStaticDataNonlocal());
			break;
		default:
			throw MechanicsException("[NuTo::ElementDataIpBase::ElementDataIpBase] Ip data type not known.");
		}
		mIpData[theIp].Initialize(rElement,rElement->GetConstitutiveLaw(theIp));
	}
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataIpBase::GetIntegrationType()const
{
	return mIntegrationType;
}

