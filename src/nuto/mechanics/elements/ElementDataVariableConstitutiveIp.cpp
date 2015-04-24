// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementDataVariableConstitutiveIp.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include <assert.h>


NuTo::ElementDataVariableConstitutiveIp::ElementDataVariableConstitutiveIp(const ElementBase *rElement,
		const NuTo::IntegrationTypeBase* rIntegrationType,	NuTo::IpData::eIpDataType rIpDataType) :
		NuTo::ElementDataBase::ElementDataBase(), ElementDataVariableConstitutiveBase() , ElementDataIpBase(rElement,rIntegrationType,rIpDataType)
{
	//std::cout << "NuTo::ElementDataVariableConstitutiveIp::ElementDataVariableConstitutiveIp" << std::endl;
}

NuTo::ElementDataVariableConstitutiveIp::ElementDataVariableConstitutiveIp(const ElementBase *rElement,
		int rNumIp,	NuTo::IpData::eIpDataType rIpDataType) :
		NuTo::ElementDataBase::ElementDataBase(), ElementDataVariableConstitutiveBase() , ElementDataIpBase(rElement,rNumIp,rIpDataType)
{
	//std::cout << "NuTo::ElementDataVariableConstitutiveIp::ElementDataVariableConstitutiveIp" << std::endl;
}

NuTo::ElementDataVariableConstitutiveIp::~ElementDataVariableConstitutiveIp()
{
//	std::cout << "NuTo::ElementDataConstitutiveStaticData::~ElementDataConstitutiveStaticData()" << std::endl;
}

//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
//! @param rElement element
void NuTo::ElementDataVariableConstitutiveIp::InitializeUpdatedConstitutiveLaw(const ElementBase* rElement,
		int rNumIp)
{
	//reinitialize ip data (f.e. if different static data or nonlocal data are required with another constitutive model)
	mIpData[rNumIp].Initialize(rElement,mVarConstitutiveLaw[rNumIp]);

}
/*
//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
//! @param rElement element
void NuTo::ElementDataVariableConstitutiveIp::InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)
{
	//reinitialize ip data (f.e. if different static data or nonlocal data are required with another constitutive model)
	for (int theIp=0; theIp<(int)mIpData.size();theIp++)
	{
		mIpData[theIp].Initialize(rElement,&mVarConstitutiveLaw[0]);
		if(mVarConstitutiveLaw[theIp]==0)
			mVarConstitutiveLaw.replace[theIp].Initialize(theIp,&mVarConstitutiveLaw[0]);
	}
}
*/

//! @brief returns the enum of element data type
//! @return enum of ElementDataType
const NuTo::ElementData::eElementDataType NuTo::ElementDataVariableConstitutiveIp::GetElementDataType()const
{
    return NuTo::ElementData::VARIABLECONSTITUTIVELAWIP;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataVariableConstitutiveIp::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataVariableConstitutiveIp::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataVariableConstitutiveIp" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataVariableConstitutiveBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataIpBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataVariableConstitutiveIp" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataVariableConstitutiveIp)
#endif // ENABLE_SERIALIZATION
