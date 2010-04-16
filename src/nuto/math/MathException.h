// $Id$

#ifndef NUTO_MATH_EXCEPTION_H
#define NUTO_MATH_EXCEPTION_H

#include "nuto/base/Exception.h"

namespace NuTo
{
//! @author J�rg F. Unger, ISM
//! @date July 2008
//! @brief ... class for all exceptions thrown in math module of NuTo
class MathException : public NuTo::Exception
{
public:
    //! @brief ...constructor
    //! @param message_ ...error message
    explicit MathException(const std::string& message_) : NuTo::Exception(message_) {}

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit MathException(const std::string& message_, bool FatalFlag_) : NuTo::Exception(message_, FatalFlag_) {}

    //! @brief ...destructor
    virtual ~MathException() throw() {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    MathException* Clone();
};
} //namespace NuTo
#endif //NUTO_MATH_EXCEPTION_H
