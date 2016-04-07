#ifndef NUTO_EXCEPTION_H
#define NUTO_EXCEPTION_H

#ifdef HAVE_BOOST
#include <boost/serialization/vector.hpp>
#endif

#include <string>
#include <exception>
#include <iostream>

namespace NuTo
{
//! @author J�rg F. Unger, ISM
//! @date July 2008
//! @brief ... base class for all exceptions thrown in NuTo
class Exception : public std::exception
{
protected:
    std::string mMessage;   //!< error message
    bool mFatalFlag;        //!< flag to decide, if throwing this exception leads to inconsistencies of the program -> fatal flag=true

public:

    //! @brief ...constructor
    //! @param rMessage ...error message
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rMessage, bool rFatalFlag = true) :
        mMessage(rMessage), mFatalFlag(rFatalFlag)
    {}

    //! @brief ...constructor
    //! @param rCaller ... name of the method that throws
    //! @param rMessage ...error message
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rCaller, const std::string& rMessage, bool rFatalFlag = true) :
        mMessage(std::string("[") + rCaller +"]\n" + rMessage), mFatalFlag(rFatalFlag)
    {}

    //! @brief ...constructor
    //! @param rCaller ... name of the method that throws
    //! @param rMessage ...error message // overload of const char*, otherwise std::string would be converted to bool and the first ctor is called.
    //! @param rFatalFlag ...flag to decide, if the error is fatal or not (exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rCaller, const char* rMessage, bool rFatalFlag = true) :
        mMessage(std::string("[") + rCaller +"]\n" + rMessage), mFatalFlag(rFatalFlag)
    {}


    //! @brief ...destructor
    virtual ~Exception() throw() {}

    Exception(const Exception& ) = default;

    //! @brief ... return error message of the exception
    //! @return ... error message
    virtual std::string ErrorMessage() const throw();

    const char* what() const throw()
    {
        // return ErrorMessage().c_str();
        //
        // For uncaught exceptions mMessage is out of scope. Thus
        // the c_str-pointer to its data points nowhere.
        // A dynamic allocation solves this problem.
        std::string* error = new std::string(ErrorMessage());
        return error->c_str();
    }

    //! @brief ... add a message to the exception (to be able to rethrow the exception afterwards)
    //! @param message_ ... message to add
    void AddMessage(const std::string &message_);

    //! @brief ...constructor
    //! @param rCaller ... name of the method that throws
    //! @param rMessage ...error message
    void AddMessage(const std::string& rCaller, const std::string& rMessage);

    //! @brief ... is the exception fatal (inconsistency of the software?)
    //! @return ... true or false
    bool IsFatal()const;

    //! @brief ... clone the exception (important, if called from the base class)
    //! @return ... a copy of the exception
    virtual Exception* Clone()
    {
        return new Exception(*this);
    }


};
} //namespace NuTo

#endif  // NUTO_EXCEPTION_H
