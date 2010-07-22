// $Id: SparseMatrixCSRGeneral.cpp 195 2009-12-16 09:13:29Z eckardt4 $

#include <boost/spirit/include/classic_core.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "nuto/math/Matrix.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/MathException.h"

namespace NuTo
{
//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixDouble
template<>
std::string SparseMatrixCSRVector2General<double>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2GeneralDouble");
}

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name FullMatrixInt
template<>
std::string SparseMatrixCSRVector2General<int>::GetTypeId()const
{
    return std::string("SparseMatrixCSRVector2GeneralInt");
}
}
