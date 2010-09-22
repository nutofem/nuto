// $Id: ConstitutiveTangentNonlocal3x3.cpp 102 2009-11-11 10:47:23Z eckardt4 $
#include <cstring>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal3x3.h"

// constructor
NuTo::ConstitutiveTangentNonlocal3x3::ConstitutiveTangentNonlocal3x3(int rNumNonlocalElements) : NuTo::ConstitutiveTangentBase::ConstitutiveTangentBase()
{
	mNonlocalMatrices.resize(rNumNonlocalElements);
}

// destructor
NuTo::ConstitutiveTangentNonlocal3x3::~ConstitutiveTangentNonlocal3x3()
{
}

// get number of rows
unsigned int NuTo::ConstitutiveTangentNonlocal3x3::GetNumberOfRows() const
{
    return mNonlocalMatrices.size() * 3;
}

// get number of columns
unsigned int NuTo::ConstitutiveTangentNonlocal3x3::GetNumberOfColumns() const
{
    return 3;
}

// get tangent
const double* NuTo::ConstitutiveTangentNonlocal3x3::GetData() const
{
    throw MechanicsException("[NuTo::ConstitutiveTangentNonlocal3x3::GetData] Method not implemented for nonlocal matrix.");
}

//! @brief ... get a local submatrix
//! @param ... rSubMatrix number of the submatrix
//! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
NuTo::ConstitutiveTangentLocal3x3* NuTo::ConstitutiveTangentNonlocal3x3::GetSubMatrix(int rSubMatrix)
{
	assert(rSubMatrix<(int)mNonlocalMatrices.size());
	return &mNonlocalMatrices[rSubMatrix];
}

//! @brief ... get a local submatrix
//! @param ... rSubMatrix number of the submatrix
//! @return ... pointer to the tangent submatrix 3x3 matrix (column major storage)
const NuTo::ConstitutiveTangentLocal3x3* NuTo::ConstitutiveTangentNonlocal3x3::GetSubMatrix(int rSubMatrix)const
{
	assert(rSubMatrix<(int)mNonlocalMatrices.size());
	return &mNonlocalMatrices[rSubMatrix];
}

//! @brief ... get a the number of local submatrices
//! @return ... number of submatrices
int NuTo::ConstitutiveTangentNonlocal3x3::GetNumSubMatrices()const
{
	return (int)mNonlocalMatrices.size();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveTangentNonlocal3x3::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveTangentNonlocal3x3::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveTangentNonlocal3x3" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveTangentBase)
       & BOOST_SERIALIZATION_NVP(mNonlocalMatrices)
       & BOOST_SERIALIZATION_NVP(mIsLocal);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveTangentNonlocal3x3" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveTangentNonlocal3x3)
#endif // ENABLE_SERIALIZATION
