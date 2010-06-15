// $Id: $
#include <assert.h>
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIpNonlocal.h"


NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveIpNonlocal(const ElementBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType) :
   NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase(), ElementDataNonlocalBase() , ElementDataIpBase(rElement,rIntegrationType,rIpDataType)
{
	std::cout << "NuTo::ElementDataConstitutiveIp::ElementDataConstitutiveIp" << std::endl;
}

NuTo::ElementDataConstitutiveIpNonlocal::~ElementDataConstitutiveIpNonlocal()
{
	//std::cout << "NuTo::ElementDataConstitutiveIpNonlocal::~ElementDataConstitutiveStaticDataNonlocal()" << std::endl;
}

//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
//! @param rElement element
void NuTo::ElementDataConstitutiveIpNonlocal::InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)
{
	//reinitialize ip data (f.e. if different static data or nonlocal data are required with another constitutive model)
	for (int theIp=0; theIp<(int)mIpData.size();theIp++)
	{
		mIpData[theIp].Initialize(rElement,mConstitutiveLaw);
	}
}

//! @brief adds the nonlocal weight to an integration point
//! @param rLocalIpNumber local Ip
//! @param rConstitutive constitutive model for which nonlocal data is to be calculated
//! @param rNonlocalElement element of the nonlocal ip
//! @param rNonlocalIp local ip number of the nonlocal ip
//! @param rWeight weight
void NuTo::ElementDataConstitutiveIpNonlocal::SetNonlocalWeight(int rLocalIpNumber, const ConstitutiveBase* rConstitutive,
		const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
    int oldNumNonlocalElements(mNonlocalElements.size());
	//search for nonlocal element in existing nonlocal elements (get back local element number)
	int theNonlocalElement = AddNonlocalElement(rNonlocalElement,rConstitutive);

	if (oldNumNonlocalElements!=(int)mNonlocalElements.size())
	{
	    //Got through all integration points and add zero weights
		for (int theIp=0; theIp<(int)mIpData.size(); theIp++)
		{
		    // reallocate the vector of the weights for the additional element *by seting to the zeros integration point zero weight
			mIpData[theIp].SetNonlocalWeight(theNonlocalElement,0,rNonlocalElement->GetNumIntegrationPoints(),0.);
		}
	}
	//add the weight to the IPData
	mIpData[rLocalIpNumber].SetNonlocalWeight(theNonlocalElement,rNonlocalIp,rNonlocalElement->GetNumIntegrationPoints(),rWeight);
}

//! @brief gets the nonlocal weights
//! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
//! @return vector of weights for all integration points of the nonlocal element
const std::vector<double>& NuTo::ElementDataConstitutiveIpNonlocal::GetNonlocalWeights(int rIp, int rNonlocalElement, const ConstitutiveBase* rConstitutive)const
{
    assert(rConstitutive==mConstitutive);
    std::cout<< "rIp " << rIp << "(" << mIpData.size() << ")" << std::endl;
    assert(rIp<(int)mIpData.size());
    return mIpData[rIp].GetNonlocalWeights(rNonlocalElement);
}
