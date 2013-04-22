// $Id: FullVector.h 624 2013-04-15 14:35:59Z unger3 $

#ifndef FULL_VECTOR_H
#define FULL_VECTOR_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullVector_Def.h"
#include "nuto/math/FullMatrix.h"

namespace NuTo
{
//! @brief ... default constructor
template <class T, int rows>
FullVector<T,rows>::FullVector(): FullMatrix<T,rows,1>::FullMatrix()
{
}

//! @brief ... constructor
//! @param rNumRows ... number of columns
//! @param rNumColumns ... number of rows
template <class T, int rows>
FullVector<T,rows>::FullVector( int rNumRows): FullMatrix<T,rows,1>::FullMatrix(rNumRows,1)
{
	if ( rNumRows>0 )
		this->setZero ( rNumRows);
	else
		this->resize ( rNumRows);
}

//! @brief ... constructor
//! @param rEigenMatrix ... other matrix
template<class T, int rows>
template<typename OtherDerived>
FullVector<T,rows>::FullVector ( const Eigen::MatrixBase<OtherDerived>& rEigenMatrix): NuTo::FullMatrix<T,rows,1>(rEigenMatrix)
{
}

//! @brief ... constructor
//! @brief ... creates a FullVector(means a FullVector size(m,1) )
//! @param entries_ ... vector containing the matrix in column-major orientation
template<class T, int rows>
FullVector<T,rows>::FullVector ( const std::vector<T>& entries_ ) : FullMatrix<T,rows,1>::FullMatrix()
{
	this->resize ( entries_.size(),1 );
	const T *ptr = &entries_[0];

	for ( int i=0; i<this->rows(); i++, ptr++ )   	// loop over rows
		(*this)( i,0 ) = *ptr;                     	// to access matrix coefficients,
}

//! @brief ... copy constructor
//! @param  rOther ... copied element
template<class T, int rows>
FullVector<T,rows>::FullVector ( const FullVector<T,rows>& rOther ): NuTo::FullMatrix<T,rows,1>(rOther)
{
}

//! @brief ... resize matrix (everything is deleted, attention, the eigenroutine is not setting everything to zero)
//! @param rows ... number of rows
template<class T, int rows>
void FullVector<T,rows>::Resize ( int rRows)
{
	if ( rRows!=0 )
		this->setZero ( rRows);
	else
		this->resize ( rRows);
}


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class T, int rows>
    template<class Archive>
    void FullVector<T,rows>::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of SparseMatrix" << std::endl;
#endif
        ar & boost::serialization::make_nvp("FullMatrix",boost::serialization::base_object< NuTo::FullMatrix<T,rows,1> >(*this));
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of SparseMatrix" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION




} //NAMESPACE NUTO

#endif // FULL_VECTOR_H
