// $Id$

#pragma once

#include "base/Exception.h"

namespace NuTo
{
//! @author J�rg F. Unger, ISM
//! @date July 2008
//! @brief ... class for all exceptions thrown in math module of NuTo
class MathException : public NuTo::Exception
{
public:
    //! @brief ...constructor
    //! @param rMessage ...error message
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MathException(const std::string& rMessage, bool rFatalFlag = true) :
        Exception(rMessage, rFatalFlag)
    {}

    //! @brief ...constructor
    //! @param rCaller ... name of the method that throws
    //! @param rMessage ...error message
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MathException(const std::string& rCaller, const std::string& rMessage, bool rFatalFlag = true) :
        Exception(rCaller, rMessage, rFatalFlag)
    {}

    //! @brief ...constructor
    //! @param rCaller ... name of the method that throws
    //! @param rMessage ...error message // overload of const char*, otherwise std::string would be converted to bool and the first ctor is called.
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MathException(const std::string& rCaller, const char* rMessage, bool rFatalFlag = true) :
        Exception(rCaller, rMessage, rFatalFlag)
    {}

    //! @brief ...destructor
    virtual ~MathException() throw() {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw() override
    {
        return "Exception in Module Math\n" + mMessage;
    }

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    virtual Exception* Clone() override
    {
        return new MathException(*this);
    }

};

class out_of_range : public std::out_of_range
{
public:
    out_of_range(std::string what) : std::out_of_range(what) {}
};

class invalid_argument : public std::invalid_argument
{
public:
    invalid_argument(std::string what) : std::invalid_argument(what) {}
};

} //namespace NuTo
