// $Id: SparseMatrixCSRVector2.cpp 207 2009-12-18 08:08:29Z eckardt4 $

#include <string>
#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRVector2<double>::GetTypeId() const
{
    return std::string("SparseMatrixCSRVector2Double");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRVector2<int>::GetTypeId() const
{
    return std::string("SparseMatrixCSRVector2Int");
}


}
