// $Id$
#include "visualize/VisualizeUnstructuredGrid.h"

#include <algorithm>

#include "visualize/Point.h"
#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkLine.h>
#include <vtkVertex.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>


vtkSmartPointer<vtkCell> CellToVtkCell(const NuTo::CellBase& c)
{
    vtkSmartPointer<vtkCell> cell;
    switch (c.GetCellType())
    {
    case NuTo::eCellTypes::VERTEX:
        cell = vtkSmartPointer<vtkVertex>::New();
        break;
    case NuTo::eCellTypes::HEXAHEDRON:
        cell = vtkSmartPointer<vtkHexahedron>::New();
        break;
    case NuTo::eCellTypes::LINE:
        cell = vtkSmartPointer<vtkLine>::New();
        break;
    case NuTo::eCellTypes::QUAD:
        cell = vtkSmartPointer<vtkQuad>::New();
        break;
    case NuTo::eCellTypes::TETRAEDER:
        cell = vtkSmartPointer<vtkTetra>::New();
        break;
    case NuTo::eCellTypes::TRIANGLE:
        cell = vtkSmartPointer<vtkTriangle>::New();
        break;
    }
    for (auto pointId : c.GetPointIds())
        cell->GetPointIds()->InsertNextId(pointId);
    return cell;
}

void NuTo::VisualizeUnstructuredGrid::ExportVtuDataFile(const std::string& rFilename) const
{
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    
    // geometry
    
    auto points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& point : mPoints)
        points->InsertNextPoint(point.GetCoordinates().data());
    grid->SetPoints(points);
    
    //grid->Allocate(mCells.size());
    for (const auto& cell : mCells)
    {
        auto vtkcell = CellToVtkCell(cell);
        grid->InsertNextCell(vtkcell->GetCellType(), vtkcell->GetPointIds());
    }
   
    // data

    for (const auto& pointData : mPointDataNames)
    {
        const int id = GetPointDataIndex(pointData);
        const int numComponents = mPoints[0].GetData(id).rows();
        const int numPoints = mPoints.size();
        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(pointData.c_str());
        dataArray->SetNumberOfComponents(numComponents);
        dataArray->SetNumberOfTuples(numPoints);
       
        for (int i = 0; i < numPoints; ++i)
            dataArray->SetTuple(i, mPoints[i].GetData(id).data());
        
        grid->GetPointData()->AddArray(dataArray);
    }
   
    grid->GetCellData()->vtkFieldData::Initialize();
    for (const auto& cellData : mCellDataNames)
    {
        const int id = GetCellDataIndex(cellData);
        const int numComponents = mCells[0].GetData(id).rows();
        const int numCells = mCells.size();
        auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
        dataArray->SetName(cellData.c_str());
        dataArray->SetNumberOfComponents(numComponents);
        dataArray->SetNumberOfTuples(numCells);
        
        for (int i = 0; i < numCells; ++i)
            dataArray->SetTuple(i, mCells[i].GetData(id).data());
        
        grid->GetCellData()->AddArray(dataArray);
    }

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(rFilename.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    writer->Write();
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
