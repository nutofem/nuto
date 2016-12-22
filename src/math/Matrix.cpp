// $Id$

#include "math/Matrix.h"
namespace NuTo
{
//! @brief ... convert the double entries of the matrix to strings (for file export and std::cout) - special formating
//! @param value ... entry of the matrix
//! @return string ... representation of the entry

template <>
std::string NuTo::Matrix<double>::Convert2String(double value)
{
    std::ostringstream format_message;
    format_message << std::setprecision(12) << std::setw(15) << std::showpoint << value;
    return format_message.str();
}

template <>
std::string NuTo::Matrix<int>::Convert2String(int value)
{
    std::ostringstream format_message;
    format_message << value;
    return format_message.str();
}

template <>
std::string NuTo::Matrix<double>::Convert2String(double value, bool scientific, int precision, int width) const
{
    std::ostringstream format_message;
    format_message << std::setprecision(precision);
    format_message << std::setw(width);
    format_message << std::right;
    if (scientific)
    {
        format_message << std::scientific;
    }
    else
    {
        format_message << std::fixed;
    }
    format_message << value;
    return format_message.str();
}

template <>
std::string NuTo::Matrix<int>::Convert2String(int value, bool scientific, int precision, int width) const
{
    std::ostringstream format_message;
    format_message.setf ( std::ios_base::right);
    format_message.width(width);
    format_message << value;
    return format_message.str();
}
}
