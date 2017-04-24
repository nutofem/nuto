// $Id$
#include "visualize/VisualizeUnstructuredGrid.h"

#include <algorithm>

#include "visualize/Point.h"
#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

int CellTypeToVtk(NuTo::eCellTypes cellType)
{
    switch (cellType)
    {
    case NuTo::eCellTypes::HEXAHEDRON:
        return 12;
    case NuTo::eCellTypes::LINE:
        return 3;
    case NuTo::eCellTypes::PYRAMID:
        return 14;
    case NuTo::eCellTypes::QUAD:
        return 9;
    case NuTo::eCellTypes::TETRAEDER:
        return 10;
    case NuTo::eCellTypes::TRIANGLE:
        return 5;
    case NuTo::eCellTypes::VERTEX:
        return 1;
    }
}

void NuTo::VisualizeUnstructuredGrid::ExportVtuDataFile(const std::string& rFilename) const
{
}

int NuTo::VisualizeUnstructuredGrid::AddPoint(Eigen::Vector3d coordinates)
{
    mPoints.push_back(Point(coordinates, mPointDataNames.size()));
    return mPoints.size() - 1;
}

int NuTo::VisualizeUnstructuredGrid::AddCell(std::vector<int> pointIds, eCellTypes cellType)
{
    mCells.push_back(CellBase(pointIds, cellType, mCellDataNames.size()));
    return mCells.size() - 1;
}



void NuTo::VisualizeUnstructuredGrid::CheckPoints(std::vector<int> pointIds) const
{
    for (auto pointId : pointIds)
        if (pointId >= mPoints.size())
            throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "Point id not defined.");
}

void NuTo::VisualizeUnstructuredGrid::DefinePointData(std::string name)
{
    if (std::find(mPointDataNames.begin(), mPointDataNames.end(), name) != mPointDataNames.end())
        throw VisualizeException(__PRETTY_FUNCTION__, "data name already exist for point data.");

    if (not mPoints.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding points");

    mPointDataNames.push_back(name);
}

void NuTo::VisualizeUnstructuredGrid::DefineCellData(std::string name)
{
    if (std::find(mCellDataNames.begin(), mCellDataNames.end(), name) != mCellDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data name already exist for point data.");
    
    if (not mCells.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding cells");

    mCellDataNames.push_back(name);
}

void NuTo::VisualizeUnstructuredGrid::SetPointData(int pointIndex, const std::string& name, double data)
{
    SetPointData(pointIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void NuTo::VisualizeUnstructuredGrid::SetPointData(int pointIndex, const std::string& name, Eigen::VectorXd data)
{
    mPoints[pointIndex].SetData(GetPointDataIndex(name), data); 
}

void NuTo::VisualizeUnstructuredGrid::SetCellData(int cellIndex, const std::string& name, double data)
{
    SetCellData(cellIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void NuTo::VisualizeUnstructuredGrid::SetCellData(int cellIndex, const std::string& name, Eigen::VectorXd data)
{
    mCells[cellIndex].SetData(GetCellDataIndex(name), data); 
}

int NuTo::VisualizeUnstructuredGrid::GetPointDataIndex(const std::string& name) const
{
    auto it = std::find(mPointDataNames.begin(), mPointDataNames.end(), name);
    if (it == mPointDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mPointDataNames.begin(), it);
}

int NuTo::VisualizeUnstructuredGrid::GetCellDataIndex(const std::string& name) const
{
    auto it = std::find(mCellDataNames.begin(), mCellDataNames.end(), name);
    if (it == mCellDataNames.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mCellDataNames.begin(), it);
}
