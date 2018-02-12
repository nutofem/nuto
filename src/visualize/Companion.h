#pragma once

#include <map>
#include "visualize/Visualizer.h"

namespace NuTo
{
namespace Visualize
{
class Companion
{

public:
    Companion() = default;
    Companion(std::string resultDir);

    void SetResultDirectory(std::string resultDir);

    void AddVisualizer(std::string name, Visualizer&& visualizer);
    void AddVisualizer(std::string name, Group<CellInterface> cells, const HandlerInterface& handler);

    void AddDof(std::string name, DofType dof);

    void AddCellFunction(std::string name,
                         std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> cellFunction,
                         std::string cellFunctionName);

    void AddPointFunction(std::string name, std::function<Eigen::VectorXd(Eigen::VectorXd)> pointFunction,
                          std::string pointFunctionName);

    void Plot(double t, bool asBinary = true);

private:
    struct VisualizationInfo
    {
        Visualizer mVisualizer;
        std::vector<DofType> mDofs;
        std::vector<std::pair<std::function<Eigen::VectorXd(const CellData&, const CellIpData&)>, std::string>>
                mCellFunctions;
        std::vector<std::pair<std::function<Eigen::VectorXd(Eigen::VectorXd)>, std::string>> mPointFunctions;
    };

    std::map<std::string, VisualizationInfo> mVisualize;

    std::string mResultDir;

    int mStep = 0;

    void ThrowOnUnknownName(std::string name);
};
} /* Visualize */
} /* NuTo */
