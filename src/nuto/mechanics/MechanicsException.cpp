#include <string>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
//! @brief ... return error message of the exception
//! @return ... error message
std::string MechanicsException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Mechanics\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
MechanicsException* MechanicsException::Clone()
{
    return new MechanicsException(*this);
}
} //namespace NuTo
