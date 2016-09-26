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

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/elements/ElementDataIpBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/elements/IpDataEmpty.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"


//! @brief constructor
//! @param rElement			... element for the IP Data
//! @param rIntegrationType	... integration type
//! @param rIpDataType		... the IP Data
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
        case NuTo::IpData::eIpDataType::NOIPDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]empty\n");
			mIpData.push_back(new IpDataEmpty());
			break;
        case NuTo::IpData::eIpDataType::STATICDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATA\n");
			mIpData.push_back(new IpDataStaticData());
			break;
        case NuTo::IpData::eIpDataType::STATICDATANONLOCAL:
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

//! @brief constructor
//! @param rElement			... element for the IP Data
//! @param rIntegrationType	... integration type
//! @param rIpDataType		... the IP Data
NuTo::ElementDataIpBase::ElementDataIpBase(const ElementBase *rElement, int rNumIp, NuTo::IpData::eIpDataType rIpDataType) : NuTo::ElementDataBase()
{
    //std::cout << "ElementDataIpBase constructor " << std::endl;
	mIntegrationType = 0;
	mIpData.clear();
	mIpData.reserve(rNumIp);
	for (int theIp=0; theIp<rNumIp;theIp++)
	{
		switch (rIpDataType)
		{
        case NuTo::IpData::eIpDataType::NOIPDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]empty\n");
			mIpData.push_back(new IpDataEmpty());
			break;
        case NuTo::IpData::eIpDataType::STATICDATA:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATA\n");
			mIpData.push_back(new IpDataStaticData());
			break;
        case NuTo::IpData::eIpDataType::STATICDATANONLOCAL:
			//printf("[NuTo::ElementDataIpBase::ElementDataIpBase]STATICDATANONLOCAL\n");
			mIpData.push_back(new IpDataStaticDataNonlocal());
			break;
		default:
			throw MechanicsException("[NuTo::ElementDataIpBase::ElementDataIpBase] Ip data type not known.");
		}
        //initialize data without constitutive law - store zeros as pointers e.g. to static data reallocation, when
        //constitutive model is modified
		mIpData[theIp].Initialize(rElement, (ConstitutiveBase*)0);
	}
}

//! @brief deconstructor
NuTo::ElementDataIpBase::~ElementDataIpBase()
{
//    std::cout << "NuTo::ElementDataIpBase::~ElementDataIpBase()" << std::endl;
}

//! @brief returns the weight of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIp number of ip
//! @return weight
double NuTo::ElementDataIpBase::GetIntegrationPointWeight(int rIp)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    if (mIntegrationType!=0)
        return mIntegrationType->GetIntegrationPointWeight(rIp);
    else
        return mIpData[rIp].GetIntegrationPointWeight();
}

//! @brief sets the weight of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @param weight
void NuTo::ElementDataIpBase::SetIntegrationPointWeight(int rIp, double rWeight)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    mIpData[rIp].SetIntegrationPointWeight(rWeight);
}


//! @brief Sets the integration type of an element.
//! @note Implemented with an exception for all elements, reimplementation required for those elements which actually
//! need an integration type
//! @param rElement Pointer to the element.
//! @param rIntegrationType Pointer to the integration type.
void NuTo::ElementDataIpBase::SetIntegrationType(const ElementBase* rElement,
        const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
	mIntegrationType = rIntegrationType;
	mIpData.clear();
	mIpData.reserve(mIntegrationType->GetNumIntegrationPoints());
	for (int theIp=0; theIp<mIntegrationType->GetNumIntegrationPoints();theIp++)
	{
		switch (rIpDataType)
		{
        case NuTo::IpData::eIpDataType::NOIPDATA:
			mIpData.push_back(new IpDataEmpty());
			break;
        case NuTo::IpData::eIpDataType::STATICDATA:
			mIpData.push_back(new IpDataStaticData());
			break;
        case NuTo::IpData::eIpDataType::STATICDATANONLOCAL:
			mIpData.push_back(new IpDataStaticDataNonlocal());
			break;
		default:
			throw MechanicsException("[NuTo::ElementDataIpBase::ElementDataIpBase] Ip data type not known.");
		}
		if (rElement->HasConstitutiveLawAssigned(theIp))
		{
			mIpData[theIp].Initialize(rElement,rElement->GetConstitutiveLaw(theIp));
		}
	}
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataIpBase::GetIntegrationType()const
{
	assert(mIntegrationType!=0);
	return mIntegrationType;
}

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementDataIpBase::GetNumIntegrationPoints()const
{
    return mIpData.size();
}

//! @brief returns ip data type of the element
//! implemented with an exception for all element data, reimplementation required for those element data
//! which actually need an integration type
//! @return enum to ip data
NuTo::IpData::eIpDataType NuTo::ElementDataIpBase::GetIpDataType(int  rIp)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	return mIpData[rIp].GetIpDataType();
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::Constitutive::StaticData::Component* NuTo::ElementDataIpBase::GetConstitutiveStaticData(int rIp)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    return mIpData[rIp].GetConstitutiveStaticData();
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::Constitutive::StaticData::Component* NuTo::ElementDataIpBase::GetConstitutiveStaticData(int rIp)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	return mIpData[rIp].GetConstitutiveStaticData();
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::IpDataStaticDataBase& NuTo::ElementDataIpBase::GetIpData(int rIp)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    return mIpData[rIp].GetIpData();
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::IpDataStaticDataBase& NuTo::ElementDataIpBase::GetIpData(int rIp) const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
    return mIpData[rIp].GetIpData();
}

//! @brief sets the static data for an integration point of an element
//! @param rIp integration point
//! @param rStaticData static data
void NuTo::ElementDataIpBase::SetConstitutiveStaticData(int rIp, Constitutive::StaticData::Component* rStaticData)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	mIpData[rIp].SetConstitutiveStaticData(rStaticData);
}

//! @brief returns the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return rLocalCoordinates
void NuTo::ElementDataIpBase::GetLocalIntegrationPointCoordinates2D(int rIp, boost::array<double,2 >& rLocalCoordinates)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	mIpData[rIp].GetLocalIntegrationPointCoordinates2D(rLocalCoordinates);
}

//! @brief sets the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @param rLocalCoordinates
void NuTo::ElementDataIpBase::SetLocalIntegrationPointCoordinates2D(int rIp, const boost::array<double,2 >& rLocalCoordinates)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	mIpData[rIp].SetLocalIntegrationPointCoordinates2D(rLocalCoordinates);
}

//! @brief returns the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return localCoordinatesFacet
void NuTo::ElementDataIpBase::GetLocalIntegrationPointCoordinates3D(int rIp, boost::array<double,3 >& rLocalCoordinates)const
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	mIpData[rIp].GetLocalIntegrationPointCoordinates3D(rLocalCoordinates);
}

//! @brief sets the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @param localCoordinatesFacet
void NuTo::ElementDataIpBase::SetLocalIntegrationPointCoordinates3D(int rIp, const boost::array<double,3 >& rLocalCoordinates)
{
    assert(rIp<(int)mIpData.size() && rIp>=0);
	mIpData[rIp].SetLocalIntegrationPointCoordinates3D(rLocalCoordinates);
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
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementDataIpBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataIpBase)
#endif // ENABLE_SERIALIZATION
