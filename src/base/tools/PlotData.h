#pragma once

#include <array>
#include <string>


namespace NuTo
{
namespace Plot
{

enum class eLineType;

class PlotData
{
public:
    // Ctor / Dtor
    // -----------
    PlotData() = default;
    PlotData(const PlotData&) = default;
    PlotData(PlotData&&) = default;
    ~PlotData() = default;

    PlotData(std::string title, std::array<unsigned char, 3> lineColor, eLineType lineType);

    // assignment operator
    // -------------------
    PlotData& operator=(const PlotData&) = default;
    PlotData& operator=(PlotData&&) = default;


    std::string GetSetupString() const;

private:
    std::string GetColorCodeString() const;


    eLineType mLineType;
    std::string mTitle = "";
    std::array<unsigned char, 3> mLineColor = {{0, 0, 0}};
};
} // namespace Plot
} // namespace NuTo
