#pragma once

#include <string>
#include <stdexcept>

namespace NuTo
{
//! @brief Base class for all exceptions thrown in NuTo.
class Exception : public std::exception
{
protected:
    std::string mMessage; //!< error message

public:
    //! @brief Constructor.
    //! @param message Error message.
    explicit Exception(const std::string& message)
        : mMessage(message)
    {
    }

    //! @brief Constructor.
    //! @param caller Name of the method that throws.
    //! @param message Error message.
    explicit Exception(const std::string& caller, const std::string& message)
        : mMessage(std::string("[") + caller + "]\n" + message)
    {
    }

    //! @brief Destructor.
    virtual ~Exception()
    {
    }

    //! @brief Return error message of the exception.
    //! @return Error message.
    virtual std::string ErrorMessage() const
    {
        return mMessage;
    }

    Exception(const Exception&) = default;

    //! @brief Clone the exception.
    //! @return A copy of the exception.
    virtual Exception* Clone()
    {
        return new Exception(*this);
    }

    const char* what() const noexcept override
    {
        // you can't do `return ErrorMessage().c_str();`
        // For uncaught exceptions mMessage is out of scope. Thus
        // the c_str-pointer to its data points nowhere.
        // A dynamic allocation solves this problem.ü
        std::string* error = new std::string(ErrorMessage());
        return error->c_str();
    }
};
} // namespace NuTo
