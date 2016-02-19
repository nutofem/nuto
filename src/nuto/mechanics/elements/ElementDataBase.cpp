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

#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/MechanicsException.h"

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
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief sets the constitutive law for a single integration point of the element
//! @param rConstitutiveLaw constitutive law
//! @param rIp integration point
void NuTo::ElementDataBase::SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetConstitutiveLaw] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the ElementDataClass - check the allocated element data type.");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementDataBase::GetStaticData(int rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetStaticData] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief sets the static data for an integration point of an element
//! @param rIp integration point
//! @param rStaticData static data
void NuTo::ElementDataBase::SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetStaticData] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return rLocalCoordinates
void NuTo::ElementDataBase::GetLocalIntegrationPointCoordinates2D(int rIpNum, boost::array<double,2 >& rLocalCoordinates)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetLocalIntegrationPointCoordinates2D] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @param rLocalCoordinates
void NuTo::ElementDataBase::SetLocalIntegrationPointCoordinates2D(int rIpNum, const boost::array<double,2 >& rLocalCoordinates)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetLocalIntegrationPointCoordinates2D] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return rLocalCoordinates
void NuTo::ElementDataBase::GetLocalIntegrationPointCoordinates3D(int rIpNum, boost::array<double,3 >& rLocalCoordinates)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetLocalIntegrationPointCoordinates3D] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief sets the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return rLocalCoordinates
void NuTo::ElementDataBase::SetLocalIntegrationPointCoordinates3D(int rIpNum, const boost::array<double,3 >& rLocalCoordinates)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetLocalIntegrationPointCoordinates3D] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief gets the weight of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @return weight
double NuTo::ElementDataBase::GetIntegrationPointWeight(int rIpNum)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetIntegrationPointWeight] Not implemented for the ElementDataClass - check the allocated element data type.");

}

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementDataBase::GetNumIntegrationPoints()const
{
    return 0;
}

//! @brief sets the local coordinate of an integration point
//! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
//! there is no standard integration type for the element and the local coordinates are stored directly at the ip
//! @param rIpNum number of the integration point
//! @param weight
void NuTo::ElementDataBase::SetIntegrationPointWeight(int rIpNum, double rWeight)
{
	throw MechanicsException("[NuTo::ElementDataBase::GetIntegrationPointWeight] Not implemented for the ElementDataClass - check the allocated element data type.");
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
void NuTo::ElementDataBase::SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
	throw MechanicsException("[NuTo::ElementDataBase::SetIntegrationType] Not implemented for the ElementDataClass - check the allocated element data type.");
}

//! @brief returns true, if the constitutive law has been assigned
bool NuTo::ElementDataBase::HasConstitutiveLawAssigned(int rIp)const
{
	return false;
}

//! @brief returns a pointer to the integration type of an element
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementDataBase::GetIntegrationType()const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetIntegrationType] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief returns ip data type of the element
//! implemented with an exception for all element data, reimplementation required for those element data
//! which actually need an integration type
//! @return enum to ip data
NuTo::IpData::eIpDataType NuTo::ElementDataBase::GetIpDataType(int  rIp)const
{
	throw MechanicsException("[NuTo::ElementDataBase::GetIpDataTyp] Not implemented for the ElementDataClass - check the allocated element data type..");
}

//! @brief adds the nonlocal weight to an integration point
//! @param rLocalIpNumber local Ip
//! @param rConstitutive constitutive model for which nonlocal data is to be calculated
//! @param rNonlocalElement element of the nonlocal ip
//! @param rNonlocalIp local ip number of the nonlocal ip
//! @param rWeight weight
 void NuTo::ElementDataBase::SetNonlocalWeight(int rLocalIpNumber, const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
    throw MechanicsException("[NuTo::ElementDataBase::AddNonlocalIp] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

//! @brief gets the nonlocal elements for a constitutive model
//! @param rConstitutive constitutive model
//! @return vector to nonlocal elements
const std::vector<const NuTo::ElementBase*>& NuTo::ElementDataBase::GetNonlocalElements()const
{
    throw MechanicsException("[NuTo::ElementDataBase::GetNonlocalElements] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

//! @brief delete the nonlocal elements
void NuTo::ElementDataBase::DeleteNonlocalElements()
{
    throw MechanicsException("[NuTo::ElementDataBase::DeleteNonlocalElements] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

//! @brief gets number of nonlocal elements for a constitutive model
//! @param rConstitutive constitutive model
//! @return number of nonlocal elements
int NuTo::ElementDataBase::GetNumNonlocalElements()const
{
    return 0;
}

//! @brief gets the nonlocal weights
//! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
//! @return vector of weights for all integration points of the nonlocal element
const std::vector<double>& NuTo::ElementDataBase::GetNonlocalWeights(int rIp, int rNonlocalElement)const
{
    throw MechanicsException("[NuTo::ElementDataBase::GetNonlocalWeights] Not implemented for the ElementDataBase class - check the allocated element data type..");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! functions defined in @link ElementDataCrackBase.h @endlink

//! @brief gets the cracks of an element
//! @return vector to cracks
std::vector<NuTo::CrackBase*>& NuTo::ElementDataBase::GetCracks()
{
    throw MechanicsException("[NuTo::ElementDataBase::GetCracks] Not implemented for the ElementDataBase class - check the allocated element data type..");
}
//! @brief gets the number of cracks for an element
//! @return number of cracks
int NuTo::ElementDataBase::GetNumCracks()const
{
    throw MechanicsException("[NuTo::ElementDataBase::GetNumCracks] Not implemented for the ElementDataBase class - check the allocated element data type..");
}
//! @brief adds a crack to the element
//! @param rCrack  crack
//! @return the local crack number, the crack is either append to the list, or the existing local number is returned
unsigned int NuTo::ElementDataBase::AddCrack(NuTo::CrackBase* rCrack)
{
    throw MechanicsException("[NuTo::ElementDataBase::AddCrack] Not implemented for the ElementDataBase class - check the allocated element data type..");
}
//! @brief Set the information that the element is already cracked or not
//! @param bool (Input) cracked or not
void NuTo::ElementDataBase::IsCracked(const bool rIsCracked)
{
    throw MechanicsException("[NuTo::ElementDataBase::IsCracked] Not implemented for the ElementDataBase class - check the allocated element data type..");
}
//! @brief Give the information if the element is already cracked or not
//! @return bool false: elements with this elementDataType are not cracked
const bool NuTo::ElementDataBase::IsCracked() const
{
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! @brief returns the enum of element data type
//! @return enum of ElementDataType
const NuTo::ElementData::eElementDataType NuTo::ElementDataBase::GetElementDataType()const
{
    return NuTo::ElementData::NOELEMENTDATA;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementDataBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementDataBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataBase)
#endif // ENABLE_SERIALIZATION
