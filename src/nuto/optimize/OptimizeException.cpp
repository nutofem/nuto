#include "nuto/optimize/OptimizeException.h"

//! @brief ... return error message of the exception
//! @return ... error message
std::string NuTo::OptimizeException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Optimize\n"+message);

    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
NuTo::OptimizeException* NuTo::OptimizeException::Clone()
{
    return new OptimizeException(*this);
}
