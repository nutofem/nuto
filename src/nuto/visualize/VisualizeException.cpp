// $Id$

#include "nuto/visualize/VisualizeException.h"

//! @brief ... return error message of the exception
//! @return ... error message
std::string NuTo::VisualizeException::ErrorMessage() const throw()
{
    std::string tmp_message("Exception in Module Visualize\n"+message);
    return tmp_message;
}

//! @brief ... clone the exception (important, if called from the base class)
//! @return ... a copy of the exception
NuTo::Exception* NuTo::VisualizeException::Clone()
{
    return new VisualizeException(*this);
}
