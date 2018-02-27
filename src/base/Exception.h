#pragma once

#include <string>
#include <stdexcept>

namespace NuTo
{
//! @brief Base class for all exceptions thrown in NuTo.
class Exception : public std::runtime_error
{
public:
    //! @brief Constructor.
    //! @param message Error message.
    explicit Exception(const std::string& message)
        : std::runtime_error(message)
    {
    }

    //! @brief Constructor.
    //! @param caller Name of the method that throws.
    //! @param message Error message.
    explicit Exception(const std::string& caller, const std::string& message)
        : std::runtime_error(std::string("[") + caller + "]\n" + message)
    {
    }
};
} // namespace NuTo
