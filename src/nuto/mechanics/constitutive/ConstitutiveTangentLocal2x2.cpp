// $Id: $
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

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal2x2.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::ConstitutiveTangentLocal2x2::ConstitutiveTangentLocal2x2() : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
    mTangent[0] = 0.;
    mTangent[1] = 0.;
    mTangent[2] = 0.;
    mTangent[3] = 0.;
}

//! @brief ... copy constructor from matrix
NuTo::ConstitutiveTangentLocal2x2& NuTo::ConstitutiveTangentLocal2x2::operator= (const NuTo::FullMatrix<double>& rOtherMatrix)
{
    if (rOtherMatrix.GetNumRows()!=2 || rOtherMatrix.GetNumColumns()!=2)
        throw MechanicsException("[NuTo::ConstitutiveTangentLocal2x2::operator=] matrix has to have dimension of 2x2.");

    memcpy(mTangent,rOtherMatrix.mEigenMatrix.data(),4*sizeof(double));
    return *this;
}

// destructor
NuTo::ConstitutiveTangentLocal2x2::~ConstitutiveTangentLocal2x2()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentLocal2x2::GetNumberOfRows() const
{
    return 2;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentLocal2x2::GetNumberOfColumns() const
{
    return 2;
}

// get tangent
const double* NuTo::ConstitutiveTangentLocal2x2::GetData() const
{
    return mTangent;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentLocal2x2::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveTangentLocal2x2::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveTangentLocal2x2" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
       & BOOST_SERIALIZATION_NVP(mTangent);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveTangentLocal2x2" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveTangentLocal2x2)
#endif // ENABLE_SERIALIZATION
