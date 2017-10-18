#include "base/tools/PlotData.h"

#include "base/tools/PlotEnum.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

NuTo::Plot::PlotData::PlotData(std::string title, std::array<unsigned char, 3> lineColor, eLineType lineType)
    : mLineType(lineType)
    , mTitle(title)
    , mLineColor(lineColor)
{
}

std::string NuTo::Plot::PlotData::GetSetupString() const
{
    std::string setupString;
    std::string lineType = LineTypeTypeToString(mLineType);
    std::transform(lineType.begin(), lineType.end(), lineType.begin(), ::tolower);

    setupString.append("w ");
    setupString.append(lineType);
    setupString.append(" title '");
    setupString.append(mTitle);
    setupString.append("' lt rgb '");
    setupString.append(GetColorCodeString());
    setupString.append("'");

    return setupString;
}

std::string NuTo::Plot::PlotData::GetColorCodeString() const
{
    // https://stackoverflow.com/questions/7639656/getting-a-buffer-into-a-stringstream-in-hex-representation/7639754#7639754
    std::stringstream ss;
    ss << "#" << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < mLineColor.size(); ++i)
        ss << std::setw(2) << static_cast<unsigned int>(mLineColor[i]);
    return ss.str();
}
