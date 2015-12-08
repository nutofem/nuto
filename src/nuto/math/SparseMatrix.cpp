// $Id$

#include "nuto/math/SparseMatrix.h"
namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrix<double>::GetTypeId() const
{
    return std::string("SparseMatrixDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrix<int>::GetTypeId() const
{
    return std::string("SparseMatrixInt");
}

}

