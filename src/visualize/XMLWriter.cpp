#include "visualize/XMLWriter.h"
#include "visualize/UnstructuredGrid.h"

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

vtkSmartPointer<vtkPoints> DefineVtkPoints(const std::vector<NuTo::Visualize::Point>& points)
{
    auto p = vtkSmartPointer<vtkPoints>::New();
    p->SetNumberOfPoints(points.size());
    for (int i = 0; i < points.size(); ++i)
        p->SetPoint(i, points[i].GetCoordinates().data());
    return p;
}

void DefineVtkCells(const std::vector<NuTo::Visualize::Cell>& cells, vtkUnstructuredGrid& grid)
{
    for (const auto& cell : cells)
    {
        auto vtkcellType = ToVtkCellType(cell.GetCellType());
        auto ids = vtkSmartPointer<vtkIdList>::New();
        const int numIds = cell.GetPointIds().size();
        ids->SetNumberOfIds(numIds);
        for (int i = 0; i < numIds; ++i)
            ids->SetId(i, cell.GetPointIds()[i]);
        grid.InsertNextCell(vtkcellType, ids);
    }
}

Eigen::VectorXd TransformData(Eigen::VectorXd data)
{
    //                     0  1  2  3  4  5
    // NuTo voigt format: xx yy zz yz xz xy
    std::swap(data[3], data[5]); 
    // swap to:           xx yy zz xy xz yz
    std::swap(data[4], data[5]); 
    // VTK voigt format:  xx yy zz xy yz xz
    return data;
}

template <typename Ts>
vtkSmartPointer<vtkDoubleArray> DefineDataArray(const Ts data, const std::string& name, int index)
{
    const int num = data.size();
    const int numComponents = data[0].GetData(index).rows();

    auto dataArray = vtkSmartPointer<vtkDoubleArray>::New();
    dataArray->SetName(name.c_str());
    dataArray->SetNumberOfComponents(numComponents);
    dataArray->SetNumberOfTuples(num);

    if (numComponents == 6)
        // write with transformation
        for (int i = 0; i < num; ++i)
            dataArray->SetTuple(i, TransformData(data[i].GetData(index)).data());
    else
        // write without transformation
        for (int i = 0; i < num; ++i)
            dataArray->SetTuple(i, data[i].GetData(index).data());

    return dataArray;
}


void NuTo::Visualize::XMLWriter::Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool asBinary)
{
    const auto& points = unstructuredGrid.mPoints;
    const auto& cells = unstructuredGrid.mCells;
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // define geometry
    grid->SetPoints(DefineVtkPoints(points));
    DefineVtkCells(cells, *grid);

   
    // define data
    for (const auto& name : unstructuredGrid.mPointDataNames)
        grid->GetPointData()->AddArray(DefineDataArray(points, name, unstructuredGrid.GetPointDataIndex(name)));

    for (const auto& name : unstructuredGrid.mCellDataNames)
        grid->GetCellData()->AddArray(DefineDataArray(cells, name, unstructuredGrid.GetCellDataIndex(name)));

    
    // write to file
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);

    if (asBinary)
        writer->SetDataModeToBinary();
    else
        writer->SetDataModeToAscii();

    writer->Write();
}
