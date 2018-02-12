#include "visualize/Companion.h"
#include <boost/filesystem.hpp>

using namespace NuTo;
using namespace NuTo::Visualize;

Companion::Companion(std::string resultDir)
{
    SetResultDirectory(resultDir);
}

void Companion::SetResultDirectory(std::string resultDir)
{
    namespace fs = boost::filesystem;

    fs::path p(resultDir);

    if (fs::exists(p) and fs::is_directory(p))
        fs::remove_all(p);

    if (fs::exists(p) and fs::is_regular_file(p))
        throw Exception(__PRETTY_FUNCTION__, "Directory name already exists as a file.");

    if (not fs::exists(p))
        fs::create_directory(p);
}

void Companion::ThrowOnUnknownName(std::string name)
{
    if (mVisualize.find(name) == mVisualize.end())
        throw Exception(__PRETTY_FUNCTION__, "You have to call AddVisualizer for " + name + " first.");
}

void Companion::AddVisualizer(std::string name, Group<CellInterface> cells, const HandlerInterface& handler)
{
    mVisualize[name].mVisualizer = Visualizer(cells, handler);
}

void Companion::AddVisualizer(std::string name, Visualizer&& visualizer)
{
    mVisualize[name].mVisualizer = std::move(visualizer);
}

void Companion::AddDof(std::string name, DofType dof)
{
    ThrowOnUnknownName(name);
    mVisualize[name].mDofs.push_back(dof);
}

void Companion::AddCellFunction(std::string name,
                                std::function<Eigen::VectorXd(const CellData&, const CellIpData&)> cellFunction,
                                std::string cellFunctionName)
{
    ThrowOnUnknownName(name);
    mVisualize[name].mCellFunctions.push_back({cellFunction, cellFunctionName});
}

void Companion::AddPointFunction(std::string name, std::function<Eigen::VectorXd(Eigen::VectorXd)> pointFunction,
                                 std::string pointFunctionName)
{
    ThrowOnUnknownName(name);
    mVisualize[name].mPointFunctions.push_back({pointFunction, pointFunctionName});
}

std::string FormatTime(double t)
{
    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << t;
    return timeFormatted.str();
}

void WritePvdFile(std::string pvdFileName, std::string vtuFileName, double t, bool isNew)
{

    std::fstream file;
    if (isNew)
        file.open(pvdFileName, std::fstream::out);
    else
        file.open(pvdFileName, std::fstream::out | std::fstream::in | std::ios_base::ate);

    if (not file.is_open())
        throw Exception(__PRETTY_FUNCTION__, "Error opening file " + pvdFileName);

    std::stringstream endOfXML;
    endOfXML << "</Collection>" << std::endl;
    endOfXML << "</VTKFile>" << std::endl;
    if (isNew)
    {
        // header /////////////////////////////////////////////////////////////////
        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"Collection\">" << std::endl;
        file << "<Collection>" << std::endl;
    }
    else
    {
        // delete the last part of the xml file
        file.seekp(-endOfXML.str().length(), std::ios_base::end);
    }
    file << "<DataSet timestep=\"" << FormatTime(t) << "\" file=\"" << vtuFileName << "\"/>\n";
    file << endOfXML.str();
}

void Companion::Plot(double t, bool asBinary)
{
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
        boost::filesystem::path vtuFile = resultFile;
        vtuFile /= visuInfo.first + std::to_string(mStep) + ".vtu";

        boost::filesystem::path pvdFile = resultFile;
        pvdFile /= visuInfo.first + ".pvd";

        visu.WriteVtuFile(vtuFile.string(), asBinary);

        WritePvdFile(pvdFile.string(), vtuFile.string(), t, mStep == 0);
    }
    mStep++;
}
