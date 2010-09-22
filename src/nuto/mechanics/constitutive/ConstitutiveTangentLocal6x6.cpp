// $Id: ConstitutiveTangentLocal6x6.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"

// constructor
NuTo::ConstitutiveTangentLocal6x6::ConstitutiveTangentLocal6x6() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
	for (int count=0; count<36; count++)
		mTangent[count]=0.;
}

// destructor
NuTo::ConstitutiveTangentLocal6x6::~ConstitutiveTangentLocal6x6()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal6x6::GetNumberOfRows() const
{
    return 6;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal6x6::GetNumberOfColumns() const
{
    return 6;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal6x6::GetData() const
{
    return mTangent;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal6x6::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveTangentLocal6x6::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveTangentLocal6x6" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
       & BOOST_SERIALIZATION_NVP(mTangent);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveTangentLocal6x6" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveTangentLocal6x6)
#endif // ENABLE_SERIALIZATION
