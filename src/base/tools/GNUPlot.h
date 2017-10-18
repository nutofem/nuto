#pragma once

#include "base/tools/PlotEnum.h"
#include "eigen3/Eigen/Core"
#include <cstdio>
#include <string>
#include <vector>

namespace NuTo
{
namespace Plot
{

class PlotData;

class GNUPlot
{
public:
    // Ctor / Dtor
    // -----------
    GNUPlot(bool keepOpen = false);
    GNUPlot(const GNUPlot&) = delete;
    GNUPlot(GNUPlot&&) = delete;
    ~GNUPlot();

    // assignment operator
    // -------------------
    GNUPlot& operator=(const GNUPlot&) = delete;
    GNUPlot& operator=(GNUPlot&&) = delete;

    void AddPlot(Eigen::VectorXd x, Eigen::VectorXd y, std::array<unsigned char, 3> lineColor = {0, 0, 0},
                 eLineType lineType = eLineType::LINES, std::string title = "");

    void Show() const;
    void Clear();

private:
    std::FILE* mStream;

    std::vector<PlotData> mSetupCurrentPlots;
};
} // namespace Plot
} // namespace NuTo
