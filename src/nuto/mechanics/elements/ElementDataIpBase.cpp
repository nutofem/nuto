// $Id$ 
// ElementDataIpDataBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataEmpty.h"
#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

//! @brief constructor
NuTo::ElementDataIpBase::ElementDataIpBase(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType) : NuTo::ElementDataBase()
{
    //std::cout << "ElementDataIpBase constructor " << std::endl;
	mIntegrationType = rIntegrationType;
	mIpData.clear();
	mIpData.reserve(mIntegrationType->GetNumIntegrationPoints());
	for (int theIp=0; theIp<mIntegrationType->GetNumIntegrationPoints();theIp++)
	{
		switch (rIpDataType)
		{
		case NuTo::IpData::NOIPDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]empty\n");
			mIpData.push_back(new IpDataEmpty());
			break;
		case NuTo::IpData::STATICDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATA\n");
			mIpData.push_back(new IpDataStaticData());
			break;
		case NuTo::IpData::STATICDATANONLOCAL:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATANONLOCAL\n");
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

//! @brief sets the fine scale model (deserialization from a binary file)
void NuTo::ElementDataIpBase::SetFineScaleModel(int rIp, std::string rFileName)
{
    assert(rIp<mIntegrationType->GetNumIntegrationPoints());
    mIpData[rIp].SetFineScaleModel(rFileName);
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::ElementDataIpBase::SetFineScaleParameter(int rIp, const std::string& rName, double rParameter)
{
    assert(rIp<mIntegrationType->GetNumIntegrationPoints());
    mIpData[rIp].SetFineScaleParameter(rName, rParameter);
}

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rElement pointer to element
//! @param rIntegrationType pointer to integration type
void NuTo::ElementDataIpBase::SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
	mIntegrationType = rIntegrationType;
	mIpData.clear();
	mIpData.reserve(mIntegrationType->GetNumIntegrationPoints());
	for (int theIp=0; theIp<mIntegrationType->GetNumIntegrationPoints();theIp++)
	{
		switch (rIpDataType)
		{
		case NuTo::IpData::NOIPDATA:
			mIpData.push_back(new IpDataEmpty());
			break;
		case NuTo::IpData::STATICDATA:
			mIpData.push_back(new IpDataStaticData());
			break;
		case NuTo::IpData::STATICDATANONLOCAL:
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
//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataIpBase::GetStaticData(int rIp)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    return mIpData[rIp].GetStaticData();
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataIpBase::GetStaticData(int rIp)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	return mIpData[rIp].GetStaticData();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataIpBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataIpBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataIpBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataIpBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataIpBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataIpBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataIpBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataIpBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase)
       & BOOST_SERIALIZATION_NVP(mIntegrationType)
       & BOOST_SERIALIZATION_NVP(mIpData);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataIpBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataIpBase)
#endif // ENABLE_SERIALIZATION
