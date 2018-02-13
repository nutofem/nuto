#pragma once

#include <map>
#include "visualize/Visualizer.h"

namespace NuTo
{
namespace Visualize
{

//! Allows you do define multiple visualizers with various visualization objects and writes them at once (including a
//! *pvd file)
class PostProcess
{
public:
    //! default ctor
    PostProcess() = default;

    //! Ctor that sets the result directory.
    //! @param resultDir Result directory.
    PostProcess(std::string resultDir);

    //! Setter for the result directory.
    //! @param resultDir Result directory.
    void ResultDirectory(std::string resultDir);

    //! Getter for the result directory.
    //! @return Result directory.
    std::string ResultDirectory() const;

    //! Moves a visualizer into the PostProcess.
    //! @param name Unique (new) visualizer name.
    //! @param visualizer (Unnamed) visualizer.
    void DefineVisualizer(std::string name, Visualizer&& visualizer);

    //! Constructs a visualizer into the PostProcess.
    //! @param name Unique (new) visualizer name.
    //! @param cells Group of cells you want to visualize.
    //! @param handler Implementation of the HandlerInterface.
    void DefineVisualizer(std::string name, Group<CellInterface> cells, const HandlerInterface& handler);

    //! Adds a dof visualization.
    //! @param name Name of an existing visualizer.
    //! @param dof DofType to visualize.
    void Add(std::string name, DofType dof);

    //! Add a cell data function that should be visualized.
    //! @param name Name of an existint visualizer.
    //! @param f Function that is passed to the cell for evaluation.
    //! @param cellFunctionName Name to be used in the resulting output file for the data array.
    void Add(std::string name, std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> cellFunction,
             std::string cellFunctionName);

    //! Visualize a function y = f(x) over a collection of cells
    //! @param name Name of an existint visualizer.
    //! @param f Function taking the coordinates as an Eigen vector and returning an Eigen vector
    //! @param pointFunctionName Name to be used in the resulting output file for the data array.
    void Add(std::string name, std::function<Eigen::VectorXd(Eigen::VectorXd)> pointFunction,
             std::string pointFunctionName);

    //! Writes one result file per defined visualizer and adds this result file to [visualizerName].pvd
    //! @param t Global time for the pvd file
    //! @param asBinary True for output as binary vtu file. False for output as ASCII.
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
