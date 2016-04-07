// $Id$

#include "nuto/base/Exception.h"

//! @brief ... return error message of the exception
//! @return ... error message
std::string NuTo::Exception::ErrorMessage() const throw()
{
    return mMessage;
}

//! @brief ... add a message to the exception (to be able to rethrow the exception afterwards)
//! @param message_ ... message to add
void NuTo::Exception::AddMessage(const std::string &rMessage)
{
    mMessage += "\n" + rMessage;
}

//! @brief ...constructor
//! @param rCaller ... name of the method that throws
//! @param rMessage ...error message
void NuTo::Exception::AddMessage(const std::string& rCaller, const std::string& rMessage)
{
    mMessage += std::string("\n[") + rCaller + "]\n" + rMessage;
}

//! @brief ... is the exception fatal (inconsistency of the software?)
//! @return ... true or false
bool NuTo::Exception::IsFatal()const
{
    return mFatalFlag;
}
