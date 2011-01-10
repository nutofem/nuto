// $Id: SparseMatrixCSRSymmetric.cpp 195 2009-12-16 09:13:29Z eckardt4 $

#include <string>

#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRVector2Symmetric<double>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2SymmetricDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRVector2Symmetric<int>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2SymmetricInt");
}
}
