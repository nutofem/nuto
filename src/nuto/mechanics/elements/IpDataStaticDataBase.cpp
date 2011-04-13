// $Id$ 
// IpDataStaticDataBase.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

NuTo::IpDataStaticDataBase::IpDataStaticDataBase() : IpDataBase()
{
    mStaticData=0;
}

NuTo::IpDataStaticDataBase::~IpDataStaticDataBase()
{
}

//! @brief sets the fine scale model (deserialization from a binary file)
void NuTo::IpDataStaticDataBase::SetFineScaleModel(std::string rFileName, double rMacroLength)
{
    mStaticData->SetFineScaleModel(rFileName, rMacroLength);
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataStaticDataBase::SetFineScaleParameter(const std::string& rName, double rParameter)
{
    mStaticData->SetFineScaleParameter(rName, rParameter);
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataStaticDataBase::SetFineScaleParameter(const std::string& rName, std::string rParameter)
{
    mStaticData->SetFineScaleParameter(rName, rParameter);
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataStaticDataBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataStaticDataBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
       & BOOST_SERIALIZATION_NVP(mStaticData);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataStaticDataBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataStaticDataBase)
#endif // ENABLE_SERIALIZATION
