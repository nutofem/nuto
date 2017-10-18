#pragma once

#include <map>
#include <string>
namespace NuTo
{
namespace Plot
{
enum class eLineType
{
    LINES,
    LINESPOINTS,
    POINTS
};

const std::map<eLineType, std::string> GetLineTypeMap();
std::string LineTypeTypeToString(eLineType rOutput);
eLineType LineTypeToEnum(std::string rOutput);

} // namespace Plot
} // namespace NuTo
