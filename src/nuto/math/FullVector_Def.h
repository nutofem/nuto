// $Id: FullVector.h 623 2013-04-05 08:14:22Z unger3 $

#pragma once

#include "nuto/math/FullMatrix_Def.h"

namespace NuTo
{
class Logger;

//! @author Jörg F. Unger, BAM
//! @date July 2013
//! @brief ... class for full vectors derived from the base class FullVector
template <class T, int rows=Eigen::Dynamic>
class FullVector : public FullMatrix<T,rows,1>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif

public:
    //! @brief ... default constructor
    FullVector<T,rows>();

    //! @brief ... constructor
    //! @param rNumRows ... number of columns
    FullVector<T,rows> ( int rNumRows);

#ifndef SWIG
    //! @brief ... constructor
    //! @param rEigenMatrix ... other matrix
    template<typename OtherDerived>
    FullVector<T,rows> ( const Eigen::MatrixBase<OtherDerived>& rOther);
#endif

    //! @brief ... constructor
    //! @brief ... creates a FullVector(means a FullVector size(m,1) )
    //! @param entries_ ... vector containing the matrix in column-major orientation
    FullVector<T,rows> ( const std::vector<T>& entries_ );

    //! @brief ... copy constructor
    //! @param  rOther ... copied element
    FullVector<T,rows> ( const FullVector<T,rows>& rOther );

    using FullMatrix<T,rows,1>::Resize;

    //! @brief ... resize matrix (everything is deleted)
    //! @param rows ... number of rows
    void Resize ( int rRows);

    //! @brief ... resize matrix (the existing data is kept)
    //! @param rows ... number of rows
    void ConservativeResize ( int rRows);


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

}; //class FullVector



}//namespace nuto

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::FullVector<double,Eigen::Dynamic>)))
BOOST_CLASS_EXPORT_KEY(BOOST_IDENTITY_TYPE((NuTo::FullVector<int,Eigen::Dynamic>)))
#endif // SWIG
#endif // ENABLE_SERIALIZATION

