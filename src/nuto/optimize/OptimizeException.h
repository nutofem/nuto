// $Id$

#ifndef NUTO_OPTIMIZE_EXCEPTION_H
#define NUTO_OPTIMIZE_EXCEPTION_H

#include "nuto/base/Exception.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date July 2009
//! @brief ... class for all exceptions thrown in optimize module of NuTo
class OptimizeException : public NuTo::Exception
{
public:
    //! @brief ...constructor
    //! @param message_ ...error message
    explicit OptimizeException(const std::string& message_) : NuTo::Exception(message_) {}

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit OptimizeException(const std::string& message_, bool FatalFlag_) : NuTo::Exception(message_, FatalFlag_) {}

    //! @brief ...destructor
    virtual ~OptimizeException() throw() {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    OptimizeException* Clone();
};
} //namespace NuTo
#endif //NUTO_OPTIMIZE_EXCEPTION_H
