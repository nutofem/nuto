#ifndef NUTO_EXCEPTION_H
#define NUTO_EXCEPTION_H

#ifdef HAVE_BOOST
#include <boost/serialization/vector.hpp>
#endif

#include <string>
#include <exception>

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date July 2008
//! @brief ... base class for all exceptions thrown in NuTo
class Exception : public std::exception
{
protected:
    std::string message;   //!< error message
    bool FatalFlag;        //!< flag to decide, if throwing this exception leads to inconsistencies of the program -> fatal flag=true

public:
    //! @brief ...constructor
    //! @param message_ ...error message
    explicit Exception(const std::string&  message_) : message(message_)
    {
        FatalFlag = true;
    }

    //! @brief ...constructor
    //! @param message_ ...error message
    //! @param FatalFlag_ ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit Exception(const std::string&  message_, bool FatalFlag_) : message(message_)
    {
        FatalFlag = FatalFlag_;
    }

    //! @brief ...destructor
    virtual ~Exception() throw() {}

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    //! @brief ... add a message to the exception (to be able to rethrow the exception afterwards)
    //! @param message_ ... message to add
    void AddMessage(const std::string &message_);

    //! @brief ... is the exception fatal (inconsistency of the software?)
    //! @return ... true or false
    bool IsFatal()const;

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    virtual Exception* Clone()=0;
};
} //namespace NuTo

#endif  // NUTO_EXCEPTION_H
