#include "visualize/UnstructuredGrid.h"
#include <algorithm>
#include "visualize/VisualizeException.h"
#include "visualize/XMLWriter.h"


using namespace NuTo::Visualize;

void UnstructuredGrid::ExportVtuDataFile(const std::string& filename, bool asBinary) const
{
    XMLWriter::Export(filename, *this, asBinary);
}

int UnstructuredGrid::AddPoint(Eigen::Vector3d coordinates)
{
    mPoints.push_back(Point(coordinates, mPointDataNames.size()));
    return mPoints.size() - 1;
}

int UnstructuredGrid::AddCell(std::vector<int> pointIds, eCellTypes cellType)
{
    CheckPoints(pointIds);
    mCells.push_back(Cell(pointIds, cellType, mCellDataNames.size()));
    return mCells.size() - 1;
}



void UnstructuredGrid::CheckPoints(std::vector<int> pointIds) const
{
    for (auto pointId : pointIds)
        if (pointId >= static_cast<int>(mPoints.size()))
            throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "Point id not defined.");
}

void UnstructuredGrid::DefinePointData(std::string name)
{
    if (std::find(mPointDataNames.begin(), mPointDataNames.end(), name) != mPointDataNames.end())
        throw VisualizeException(__PRETTY_FUNCTION__, "data name already exist for point data.");

    if (not mPoints.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding points");

    mPointDataNames.push_back(name);
}

void UnstructuredGrid::DefineCellData(std::string name)
{
    if (std::find(mCellDataNames.begin(), mCellDataNames.end(), name) != mCellDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data name already exist for point data.");
    
    if (not mCells.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding cells");

    mCellDataNames.push_back(name);
}

void UnstructuredGrid::SetPointData(int pointIndex, const std::string& name, double data)
{
    SetPointData(pointIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void UnstructuredGrid::SetPointData(int pointIndex, const std::string& name, Eigen::VectorXd data)
{
    mPoints[pointIndex].SetData(GetPointDataIndex(name), data); 
}

void UnstructuredGrid::SetCellData(int cellIndex, const std::string& name, double data)
{
    SetCellData(cellIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void UnstructuredGrid::SetCellData(int cellIndex, const std::string& name, Eigen::VectorXd data)
{
    mCells[cellIndex].SetData(GetCellDataIndex(name), data); 
}

int UnstructuredGrid::GetPointDataIndex(const std::string& name) const
{
    auto it = std::find(mPointDataNames.begin(), mPointDataNames.end(), name);
    if (it == mPointDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mPointDataNames.begin(), it);
}

int UnstructuredGrid::GetCellDataIndex(const std::string& name) const
{
    auto it = std::find(mCellDataNames.begin(), mCellDataNames.end(), name);
    if (it == mCellDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mCellDataNames.begin(), it);
}
