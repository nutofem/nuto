#pragma once

#include <map>
#include "visualize/Visualizer.h"

namespace NuTo
{
namespace Visualize
{
class PostProcess
{

public:
    PostProcess() = default;
    PostProcess(std::string resultDir);

    void ResultDirectory(std::string resultDir);
    std::string ResultDirectory() const;

    void DefineVisualizer(std::string name, Visualizer&& visualizer);
    void DefineVisualizer(std::string name, Group<CellInterface> cells, const HandlerInterface& handler);

    void Add(std::string name, DofType dof);

    void Add(std::string name, std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> cellFunction,
             std::string cellFunctionName);

    void Add(std::string name, std::function<Eigen::VectorXd(Eigen::VectorXd)> pointFunction,
             std::string pointFunctionName);

    void Plot(double t, bool asBinary = true);

private:
    struct VisualizationInfo
    {
        VisualizationInfo(Visualizer&& visu)
            : mVisualizer(std::move(visu))
        {
        }

        Visualizer mVisualizer;
        std::vector<DofType> mDofs;
        std::vector<std::pair<std::function<Eigen::VectorXd(const CellData&, const CellIpData&)>, std::string>>
                mCellFunctions;
        std::vector<std::pair<std::function<Eigen::VectorXd(Eigen::VectorXd)>, std::string>> mPointFunctions;
    };

    std::map<std::string, VisualizationInfo> mVisualize;

    std::string mResultDir = "not set";

    int mStep = 0;

    void ThrowOnUnknownName(std::string name);
};
} /* Visualize */
} /* NuTo */
