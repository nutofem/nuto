// $Id$

#include "nuto/base/Exception.h"

namespace NuTo
{
//! @brief ... return error message of the exception
//! @return ... error message
std::string Exception::ErrorMessage() const throw()
{
    return message;
}

//! @brief ... add a message to the exception (to be able to rethrow the exception afterwards)
//! @param message_ ... message to add
void Exception::AddMessage(const std::string &message_)
{
    message = message_ + "\n" + message;
}

//! @brief ... is the exception fatal (inconsistency of the software?)
//! @return ... true or false
bool Exception::IsFatal()const
{
    return FatalFlag;
}
}//NuTo
