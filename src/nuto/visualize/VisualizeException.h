// $Id$

#ifndef VISUALIZEEXCEPTION_H_
#define VISUALIZEEXCEPTION_H_

#include "nuto/base/Exception.h"

namespace NuTo
{
//! @brief ... class for all exceptions thrown in visualize module of NuTo
//! @author Stefan Eckardt, ISM
//! @date 18.11.2009
class VisualizeException: public NuTo::Exception
{
public:
    //! @brief ...constructor
    //! @param message_ ...error message
    explicit VisualizeException(const std::string& message_) : NuTo::Exception(message_) {}

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit VisualizeException(const std::string& message_, bool FatalFlag_) : NuTo::Exception(message_, FatalFlag_) {}

    //! @brief ...destructor
    virtual ~VisualizeException() throw() {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    Exception* Clone();
};

}

#endif // VISUALIZEEXCEPTION_H_ 
