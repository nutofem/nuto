// $Id$ 
// IpDataBase.cpp
// created Apr 29, 2010 by Joerg F. Unger

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
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

NuTo::IpDataBase::~IpDataBase()
{
	//std::cout << "call Desctructor [NuTo::IpDataBase::~IpDataBase]." << std::endl;
}

//! @brief sets the fine scale model (deserialization from a binary file)
void NuTo::IpDataBase::SetFineScaleModel(std::string rFileName, double rMacroLength, double rCoordinates[2], std::string rIpName)
{
    throw NuTo::MechanicsException("[NuTo::IpDataBase::IpDataBaseSetFineScaleModel] This Ip data type has no fine scale model.");
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataBase::SetFineScaleParameter(const std::string& rName, double rParameter)
{
    throw NuTo::MechanicsException("[NuTo::IpDataBase::SetFineScaleParameter] This Ip data type has no fine scale model.");
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataBase::SetFineScaleParameter(const std::string& rName, std::string rParameter)
{
    throw NuTo::MechanicsException("[NuTo::IpDataBase::SetFineScaleParameter] This Ip data type has no fine scale model.");
}

#ifdef ENABLE_VISUALIZE
//Visualize for all integration points the fine scale structure
void NuTo::IpDataBase::VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const
{
    //no fine scale structure
}
#endif

//! @brief adds the weight to an integration point, eventually reallocates the data
//! @param rNonlocalElement the Element (local number from the nonlocal elements)
//! @param rNonlocalIp integration point of the nonlocal element
//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
//! @param rWeight nonlocal weight
void NuTo::IpDataBase::SetNonlocalWeight(int rElement,int rNonlocalIp,int rNumIps, double rWeight)
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::SetNonlocalWeight] This Ip data type cannot store nonlocal weights - check the ip data type.");
}

//! @brief delete the nonlocal elements
//! @param rConstitutive  constitutive model
void NuTo::IpDataBase::DeleteNonlocalWeights()
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::DeleteNonlocalWeights] This Ip data type cannot store nonlocal weights - check the ip data type.");
}

//! @brief return the nonlocal weights
//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
//! @return nonlocal weights
const std::vector<double>& NuTo::IpDataBase::GetNonlocalWeights(int rNonlocalElement)const
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::GetNonlocalWeights] This Ip data type cannot store nonlocal weights - check the ip data type.");
}

NuTo::ConstitutiveStaticDataBase* NuTo::IpDataBase::GetStaticData()
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::GetStaticData] This Ip data type has no static data.");
}

const NuTo::ConstitutiveStaticDataBase* NuTo::IpDataBase::GetStaticData()const
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::GetStaticData] This Ip data type has no static data.");
}

void NuTo::IpDataBase::SetStaticData(ConstitutiveStaticDataBase* rStaticData)
{
	throw NuTo::MechanicsException("[NuTo::IpDataBase::SetStaticData] This Ip data type has no static data.");
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
