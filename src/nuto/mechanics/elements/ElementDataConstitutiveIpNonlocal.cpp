// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <assert.h>
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIpNonlocal.h"
#include "nuto/mechanics/elements/IpDataBase.h"


NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveIpNonlocal(const ElementBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType) :
   NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase(), ElementDataNonlocalBase() , ElementDataIpBase(rElement,rIntegrationType,rIpDataType)
{
	//std::cout << "NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveStaticDataNonlocal()" << std::endl;
}

NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveIpNonlocal(const ElementBase *rElement,
		int rNumIp, NuTo::IpData::eIpDataType rIpDataType) :
   NuTo::ElementDataBase::ElementDataBase(), ElementDataConstitutiveBase(), ElementDataNonlocalBase() , ElementDataIpBase(rElement,rNumIp,rIpDataType)
{
	//std::cout << "NuTo::ElementDataConstitutiveIpNonlocal::ElementDataConstitutiveStaticDataNonlocal()" << std::endl;
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
void NuTo::ElementDataConstitutiveIpNonlocal::SetNonlocalWeight(int rLocalIpNumber,	const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
    int oldNumNonlocalElements(mNonlocalElements.size());
	//search for nonlocal element in existing nonlocal elements (get back local element number)
	int theNonlocalElement = AddNonlocalElement(rNonlocalElement);

	if (oldNumNonlocalElements!=(int)mNonlocalElements.size())
	{
	    //Go through all integration points and add zero weights
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
const std::vector<double>& NuTo::ElementDataConstitutiveIpNonlocal::GetNonlocalWeights(int rIp, int rNonlocalElement)const
{
    assert(rIp<(int)mIpData.size());
    return mIpData[rIp].GetNonlocalWeights(rNonlocalElement);
}

//! @brief returns the enum of element data type
//! @return enum of ElementDataType
const NuTo::ElementData::eElementDataType NuTo::ElementDataConstitutiveIpNonlocal::GetElementDataType()const
{
    return NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL;
}

//! @brief delete the nonlocal elements
//! @param rConstitutive  constitutive model
void NuTo::ElementDataConstitutiveIpNonlocal::DeleteNonlocalElements()
{
	mNonlocalElements.resize(0);
	for (unsigned int ip=0; ip<mIpData.size(); ip++)
		mIpData[ip].DeleteNonlocalWeights();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataConstitutiveIpNonlocal::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataConstitutiveIpNonlocal::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataConstitutiveIpNonlocal" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataConstitutiveBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataIpBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataNonlocalBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataConstitutiveIpNonlocal" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementDataConstitutiveIpNonlocal)
#endif // ENABLE_SERIALIZATION
