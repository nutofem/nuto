// $Id$

#include <string>
#include "nuto/mechanics/MechanicsException.h"

//! @brief ... return error message of the exception
//! @return ... error message
std::string NuTo::MechanicsException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Mechanics\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
NuTo::Exception* NuTo::MechanicsException::Clone()
{
    return new MechanicsException(*this);
}
