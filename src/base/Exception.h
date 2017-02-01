#pragma once

#include <string>
#include <stdexcept>
#include <iostream>

namespace NuTo
{
//! @brief Base class for all exceptions thrown in NuTo.
class Exception : public std::exception
{
protected:
    std::string mMessage;   //!< error message
    bool mFatalFlag;        //!< flag to decide, if throwing this exception leads to inconsistencies of the program
                            //!< -> `mFatalFlag=true`

public:

    //! @brief Constructor.
    //! @param rMessage Error message.
    //! @param rFatalFlag Flag to decide if the error is fatal or not.
    //!                   (Exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rMessage, bool rFatalFlag = true) :
        mMessage(rMessage), mFatalFlag(rFatalFlag)
    {}

    //! @brief Constructor.
    //! @param rCaller Name of the method that throws.
    //! @param rMessage Error message.
    //! @param rFatalFlag Flag to decide if the error is fatal or not.
    //!                   (Exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rCaller, const std::string& rMessage, bool rFatalFlag = true) :
        mMessage(std::string("[") + rCaller +"]\n" + rMessage), mFatalFlag(rFatalFlag)
    {}

    //! @brief Constructor.
    //! @param rCaller Name of the method that throws.
    //! @param rMessage Error message
    //  Overload of const char*, otherwise std::string would be converted to bool and the first ctor is called.
    //! @param rFatalFlag Flag to decide if the error is fatal or not.
    //!                   (Exit the program or able to continue in a consistent way)
    explicit Exception(const std::string& rCaller, const char* rMessage, bool rFatalFlag = true) :
        mMessage(std::string("[") + rCaller +"]\n" + rMessage), mFatalFlag(rFatalFlag)
    {}

    //! @brief Destructor.
    virtual ~Exception() throw() {}

    Exception(const Exception& ) = default;

    //! @brief Return error message of the exception.
    //! @return Error message.
    virtual std::string ErrorMessage() const throw()
    {
        return mMessage;
    }

    const char* what() const noexcept override
    {
        // return ErrorMessage().c_str();
        //
        // For uncaught exceptions mMessage is out of scope. Thus
        // the c_str-pointer to its data points nowhere.
        // A dynamic allocation solves this problem.
        std::string* error = new std::string(ErrorMessage());
        return error->c_str();
    }

    //! @brief Add a message to the exception (to be able to rethrow the exception afterwards).
    //! @param message_ Message to add.
    void AddMessage(const std::string &rMessage)
    {
        mMessage += "\n" + rMessage;
    }

    //! @brief Add a message to the exception (to be able to rethrow the exception afterwards).
    //! @param rCaller Name of the method that throws.
    //! @param message_ Message to add.
    void AddMessage(const std::string& rCaller, const std::string& rMessage)
    {
        mMessage += std::string("\n[") + rCaller + "]\n" + rMessage;
    }

    //! @brief Check whether the exception is fatal/will result in inconsistency of the program.
    bool IsFatal() const
    {
        return mFatalFlag;
    }

    //! @brief Clone the exception.
    //! @return A copy of the exception.
    virtual Exception* Clone()
    {
        return new Exception(*this);
    }
};
} //namespace NuTo
