#include "visualize/PostProcess.h"
#include <boost/filesystem.hpp>

using namespace NuTo;
using namespace NuTo::Visualize;

PostProcess::PostProcess(std::string resultDir)
{
    ResultDirectory(resultDir);
}

void PostProcess::ResultDirectory(std::string resultDir)
{
    namespace fs = boost::filesystem;

    fs::path p(resultDir);

    if (fs::exists(p) and fs::is_directory(p))
        fs::remove_all(p);

    if (fs::exists(p) and fs::is_regular_file(p))
        throw Exception(__PRETTY_FUNCTION__, "Directory name already exists as a file.");

    if (not fs::exists(p))
        fs::create_directory(p);

    mResultDir = p.string();
}

std::string PostProcess::ResultDirectory() const
{
    return mResultDir;
}

void PostProcess::ThrowOnUnknownName(std::string name)
{
    if (mVisualize.find(name) == mVisualize.end())
        throw Exception(__PRETTY_FUNCTION__, "You have to call DefineVisualizer for " + name + " first.");
}

void PostProcess::DefineVisualizer(std::string name, Group<CellInterface> cells, const HandlerInterface& handler)
{
    DefineVisualizer(name, Visualizer(cells, handler));
}

void PostProcess::DefineVisualizer(std::string name, Visualizer&& visualizer)
{
    if (mVisualize.find(name) != mVisualize.end())
        throw Exception(__PRETTY_FUNCTION__, "You already defined " + name + ". Pick a new name for a new Visualizer.");
    mVisualize.emplace(name, std::move(visualizer));
}

void PostProcess::Add(std::string name, DofType dof)
{
    ThrowOnUnknownName(name);
    mVisualize.at(name).mDofs.push_back(dof);
}

void PostProcess::Add(std::string name, std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> cellFunction,
                      std::string cellFunctionName)
{
    ThrowOnUnknownName(name);
    mVisualize.at(name).mCellFunctions.push_back({cellFunction, cellFunctionName});
}

void PostProcess::Add(std::string name, std::function<Eigen::VectorXd(Eigen::VectorXd)> pointFunction,
                      std::string pointFunctionName)
{
    ThrowOnUnknownName(name);
    mVisualize.at(name).mPointFunctions.push_back({pointFunction, pointFunctionName});
}

std::string FormatTime(double t)
{
    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << t;
    return timeFormatted.str();
}

//! Adds the vtuFileName to the pvdFileName with the time step t
//! @param pvdFileName full pvd file name including ".pvd"
//! @param vtuFileName full vtu file name including ".vtu" that is added to the pvdFileName
//! @param double t corresponding time
//! @remark This might as well be a general purpose function/class to be used outside of this class, YAGNI for now.
void AddToPvdFile(std::string pvdFileName, std::string vtuFileName, double t)
{
    boost::filesystem::path p(pvdFileName);
    bool isNew = not boost::filesystem::exists(p);

    std::fstream file;
    if (isNew)
        file.open(pvdFileName, std::fstream::out);
    else
        file.open(pvdFileName, std::fstream::out | std::fstream::in | std::ios_base::ate);

    if (not file.is_open())
        throw Exception(__PRETTY_FUNCTION__, "Error opening file " + pvdFileName);

    std::stringstream endOfXML;
    endOfXML << "</Collection>\n";
    endOfXML << "</VTKFile>\n";
    if (isNew)
    {
        // header /////////////////////////////////////////////////////////////////
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"Collection\">\n";
        file << "<Collection>\n";
    }
    else
    {
        // delete the last part of the xml file
        file.seekp(-endOfXML.str().length(), std::ios_base::end);
    }
    file << "<DataSet timestep=\"" << FormatTime(t) << "\" file=\"" << vtuFileName << "\"/>\n";
    file << endOfXML.str();
}

void PostProcess::Plot(double t, bool asBinary)
{
    if (mResultDir == "not set")
        throw Exception(__PRETTY_FUNCTION__, "You have to set the result directory by calling ResultDirectory(name).");

    for (auto& visuInfo : mVisualize)
    {
        auto& info = visuInfo.second;
        Visualizer& visu = visuInfo.second.mVisualizer;

        for (auto dof : info.mDofs)
            visu.DofValues(dof);

        for (auto cellFunction : info.mCellFunctions)
            visu.CellData(cellFunction.first, cellFunction.second);

        for (auto pointFunction : info.mPointFunctions)
            visu.PointData(pointFunction.first, pointFunction.second);

        boost::filesystem::path resultFile(mResultDir);

        std::string vtuFileName = visuInfo.first + std::to_string(mStep) + ".vtu";
        boost::filesystem::path vtuFile = resultFile / vtuFileName;
        boost::filesystem::path pvdFile = resultFile / std::string(visuInfo.first + ".pvd");

        visu.WriteVtuFile(vtuFile.string(), asBinary);

        AddToPvdFile(pvdFile.string(), vtuFileName, t);
    }
    mStep++;
}
