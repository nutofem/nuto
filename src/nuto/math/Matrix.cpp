// $Id$

#include "nuto/math/Matrix.h"
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
std::string NuTo::Matrix<double>::Convert2String(double value, bool scientific, int precision, int width)
{
    std::ostringstream format_message;
    if (scientific)
    {
        format_message.setf ( std::ios_base::right, std::ios::scientific);
    }
    else
    {
        format_message.setf ( std::ios_base::right, std::ios::fixed);
    }
    format_message.precision(precision);
    format_message.width(width);
    format_message << value;
    return format_message.str();
}

template <>
std::string NuTo::Matrix<int>::Convert2String(int value, bool scientific, int precision, int width)
{
    std::ostringstream format_message;
    format_message.setf ( std::ios_base::right);
    format_message.width(width);
    format_message << value;
    return format_message.str();
}
}
