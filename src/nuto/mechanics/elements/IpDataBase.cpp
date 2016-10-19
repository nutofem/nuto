#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>
#include <string>
#include <sstream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/IpDataBase.h"

NuTo::IpDataBase::~IpDataBase()
{
	//std::cout << "call Desctructor [NuTo::IpDataBase::~IpDataBase]." << std::endl;
}

//! @brief adds the weight to an integration point, eventually reallocates the data
//! @param rNonlocalElement the Element (local number from the nonlocal elements)
//! @param rNonlocalIp integration point of the nonlocal element
//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
//! @param rWeight nonlocal weight
void NuTo::IpDataBase::SetNonlocalWeight(int rElement,int rNonlocalIp,int rNumIps, double rWeight)
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type cannot store nonlocal weights - check the ip data type.");
}

//! @brief delete the nonlocal elements
//! @param rConstitutive  constitutive model
void NuTo::IpDataBase::DeleteNonlocalWeights()
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type cannot store nonlocal weights - check the ip data type.");
}

//! @brief return the nonlocal weights
//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
//! @return nonlocal weights
const std::vector<double>& NuTo::IpDataBase::GetNonlocalWeights(int rNonlocalElement)const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type cannot store nonlocal weights - check the ip data type.");
}

NuTo::Constitutive::StaticData::Component* NuTo::IpDataBase::GetConstitutiveStaticData()
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This IP data type has no static data.");
}

const NuTo::Constitutive::StaticData::Component* NuTo::IpDataBase::GetConstitutiveStaticData() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This IP data type has no static data.");
}

NuTo::IpDataStaticDataBase& NuTo::IpDataBase::GetIpData()
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This IP data type has no static data.");
}

const NuTo::IpDataStaticDataBase& NuTo::IpDataBase::GetIpData() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This IP data type has no static data.");
}

void NuTo::IpDataBase::SetConstitutiveStaticData(Constitutive::StaticData::Component* rStaticData)
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no static data.");
}

void NuTo::IpDataBase::GetLocalIntegrationPointCoordinates2D(boost::array<double,2 >& rLocalCoordinatesFacet)const
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no static data.");
}

void NuTo::IpDataBase::SetLocalIntegrationPointCoordinates2D(const boost::array<double,2 >& rLocalCoordinatesFacet)
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no static data.");
}

void NuTo::IpDataBase::GetLocalIntegrationPointCoordinates3D(boost::array<double,3 >& rLocalCoordinatesFacet)const
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no static data.");
}

void NuTo::IpDataBase::SetLocalIntegrationPointCoordinates3D(const boost::array<double,3 >& rLocalCoordinatesFacet)
{
	throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no static data.");
}


double NuTo::IpDataBase::GetIntegrationPointWeight() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no weight.");
}


void NuTo::IpDataBase::SetIntegrationPointWeight(double)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This Ip data type has no weight.");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataBase)
#endif // ENABLE_SERIALIZATION
