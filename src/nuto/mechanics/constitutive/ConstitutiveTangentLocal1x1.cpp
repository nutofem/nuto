// $Id: ConstitutiveTangentLocal1x1.cpp 102 2009-11-11 10:47:23Z eckardt4 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <cstring>
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"

// constructor
NuTo::ConstitutiveTangentLocal1x1::ConstitutiveTangentLocal1x1() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
    this->mTangent = 0;
}

// destructor
NuTo::ConstitutiveTangentLocal1x1::~ConstitutiveTangentLocal1x1()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal1x1::GetNumberOfRows() const
{
    return 1;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal1x1::GetNumberOfColumns() const
{
    return 1;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal1x1::GetData() const
{
    return &mTangent;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal1x1::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveTangentLocal1x1::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveTangentLocal1x1" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
       & BOOST_SERIALIZATION_NVP(mTangent);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveTangentLocal1x1" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveTangentLocal1x1)
#endif // ENABLE_SERIALIZATION
