// $Id$
#include <string>
#include "nuto/math/MathException.h"

namespace NuTo
{

//! @brief ... return error message of the exception
//! @return ... error message
std::string MathException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Math\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
MathException* MathException::Clone()
{
    return new MathException(*this);
}
} //namespace NuTo
