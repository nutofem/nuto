// $ld: $ 
// IpDataNonlocalBase.cpp
// created Apr 28, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include <assert.h>
#include <iostream>
#include "nuto/mechanics/elements/IpDataNonlocalBase.h"

NuTo::IpDataNonlocalBase::~IpDataNonlocalBase()
{
	//std::cout << std::endl << "call Desctructor [NuTo::IpDataNonlocalBase::~IpDataNonlocalBase]." << std::endl;
}


//! @brief adds the weight to an integration point, eventually reallocates the data
//! @param rNonlocalElement the Element (local number from the nonlocal elements)
//! @param rNonlocalIp integration point of the nonlocal element
//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
//! @param rWeight nonlocal weight
void NuTo::IpDataNonlocalBase::SetNonlocalWeight(int rNonlocalElement,int rNonlocalIp,int rNumIps, double rWeight)
{
    if (rNonlocalElement<(int)mNonlocalWeights.size())
    {
    	assert((int)mNonlocalWeights[rNonlocalElement].size()==rNumIps);
    	assert(rNonlocalIp<rNumIps);
    	mNonlocalWeights[rNonlocalElement][rNonlocalIp]=rWeight;
    }
    else
    {
    	// the nonlocal element has to be added
    	assert(rNonlocalElement==(int)mNonlocalWeights.size());
    	mNonlocalWeights.push_back(std::vector<double>(rNumIps,0));
    	assert(rNonlocalIp<rNumIps);
    	mNonlocalWeights[rNonlocalElement][rNonlocalIp]=rWeight;
    }
}
//! @brief return the nonlocal weights
//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
//! @return nonlocal weights
const std::vector<double>& NuTo::IpDataNonlocalBase::GetNonlocalWeights(int rNonlocalElement)const
{
    assert(rNonlocalElement>=0 && rNonlocalElement<(int)mNonlocalWeights.size());
	return mNonlocalWeights[rNonlocalElement];
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataNonlocalBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataNonlocalBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataNonlocalBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
       & BOOST_SERIALIZATION_NVP(mNonlocalWeights);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataNonlocalBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataNonlocalBase)
#endif // ENABLE_SERIALIZATION
