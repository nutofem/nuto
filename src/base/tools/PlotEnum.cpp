#include "base/tools/PlotEnum.h"

#include "base/Exception.h"
#include <algorithm>

const std::map<NuTo::Plot::eLineType, std::string> NuTo::Plot::GetLineTypeMap()
{
    const std::map<eLineType, std::string> map = {
            {eLineType::LINES, "LINES"}, {eLineType::LINESPOINTS, "LINESPOINTS"}, {eLineType::POINTS, "POINTS"}};
    return map;
}

std::string NuTo::Plot::LineTypeTypeToString(NuTo::Plot::eLineType rOutput)
{
    try
    {
        return GetLineTypeMap().find(rOutput)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

NuTo::Plot::eLineType NuTo::Plot::LineTypeToEnum(std::string rOutput)
{
    std::transform(rOutput.begin(), rOutput.end(), rOutput.begin(), ::toupper);

    for (auto entry : GetLineTypeMap())
        if (entry.second == rOutput)
            return entry.first;

    throw NuTo::Exception(__PRETTY_FUNCTION__, " Enum undefined or not implemented.");
}
