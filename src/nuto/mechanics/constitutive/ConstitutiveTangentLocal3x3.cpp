// $Id: ConstitutiveTangentLocal3x3.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>
#include <assert.h>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"

// constructor
NuTo::ConstitutiveTangentLocal3x3::ConstitutiveTangentLocal3x3() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
    mTangent[0] = 0.;
    mTangent[1] = 0.;
    mTangent[2] = 0.;
    mTangent[3] = 0.;
    mTangent[4] = 0.;
    mTangent[5] = 0.;
    mTangent[6] = 0.;
    mTangent[7] = 0.;
    mTangent[8] = 0.;
}

// destructor
NuTo::ConstitutiveTangentLocal3x3::~ConstitutiveTangentLocal3x3()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal3x3::GetNumberOfRows() const
{
    return 3;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal3x3::GetNumberOfColumns() const
{
    return 3;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal3x3::GetData() const
{
    return mTangent;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal3x3::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveTangentLocal3x3::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveTangentLocal3x3" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
       & BOOST_SERIALIZATION_NVP(mTangent);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveTangentLocal3x3" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveTangentLocal3x3)
#endif // ENABLE_SERIALIZATION
