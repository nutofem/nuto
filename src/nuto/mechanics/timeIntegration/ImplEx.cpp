/*
 * ImplEx.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#include "nuto/mechanics/timeIntegration/ImplEx.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::ImplEx::ImplEx(StructureBase* rStructure) : ImplicitExplicitBase(rStructure)
{}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ImplEx::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ImplEx::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of ImplEx" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase);

    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of ImplEx" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION
