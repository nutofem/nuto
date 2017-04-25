// $Id$
#include "visualize/VisualizeUnstructuredGrid.h"

#include <algorithm>

#include "visualize/Point.h"
#include "visualize/CellBase.h"
#include "visualize/VisualizeException.h"

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPoints.h>


int ToVtkCellType(NuTo::eCellTypes type)
{
    switch (type)
    {
    case NuTo::eCellTypes::VERTEX:
        return VTK_VERTEX;
    case NuTo::eCellTypes::HEXAHEDRON:
        return VTK_HEXAHEDRON;
    case NuTo::eCellTypes::LINE:
        return VTK_LINE;
    case NuTo::eCellTypes::QUAD:
        return VTK_QUAD;
    case NuTo::eCellTypes::TETRAEDER:
        return VTK_TETRA;
    case NuTo::eCellTypes::TRIANGLE:
        return VTK_TRIANGLE;
    case NuTo::eCellTypes::POLYGON:
        return VTK_POLYGON;
    }
}

Eigen::VectorXd TransformData(Eigen::VectorXd nuto)
{
    //                     0  1  2  3  4  5 
    // NuTo voigt format: xx yy zz yz xz xy
    // VTK voigt format:  xx yy zz xy yz xz
    
    if (nuto.rows() != 6)
        return nuto;
    
    Eigen::VectorXd vtk(6);
    vtk.segment(0,3) = nuto.segment(0,3);
    vtk[3] = nuto[5];
    vtk[4] = nuto[3];
    vtk[5] = nuto[4];
    return vtk;
}

void NuTo::VisualizeUnstructuredGrid::ExportVtuDataFile(const std::string& rFilename) const
{
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // geometry

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(mPoints.size());
    for (int i = 0; i < mPoints.size(); ++i)
        points->SetPoint(i, mPoints[i].GetCoordinates().data());
    grid->SetPoints(points);

    //grid->Allocate(mCells.size());
    for (const auto& cell : mCells)
    {
        auto vtkcellType = ToVtkCellType(cell.GetCellType());
        auto ids = vtkSmartPointer<vtkIdList>::New();
        const int numIds = cell.GetPointIds().size();
        ids->SetNumberOfIds(numIds);
        for (int i = 0; i < numIds; ++i)
            ids->SetId(i, cell.GetPointIds()[i]);
        grid->InsertNextCell(vtkcellType, ids);
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

    //grid->GetCellData()->vtkFieldData::Initialize();
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
            dataArray->SetTuple(i, TransformData(mCells[i].GetData(id)).data());

        grid->GetCellData()->AddArray(dataArray);
    }

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(rFilename.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    //writer->SetDataModeToBinary();
    writer->Write();
}

int NuTo::VisualizeUnstructuredGrid::AddPoint(Eigen::Vector3d coordinates)
{
    mPoints.push_back(Point(coordinates, mPointDataNames.size()));
    return mPoints.size() - 1;
}

int NuTo::VisualizeUnstructuredGrid::AddCell(std::vector<int> pointIds, eCellTypes cellType)
{
    CheckPoints(pointIds);
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
