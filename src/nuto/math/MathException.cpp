// $Id$
#include <string>
#include "nuto/math/MathException.h"

//! @brief ... return error message of the exception
//! @return ... error message
std::string NuTo::MathException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Math\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
NuTo::Exception* NuTo::MathException::Clone()
{
    return new MathException(*this);
}
