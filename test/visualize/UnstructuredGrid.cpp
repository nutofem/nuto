#include "BoostUnitTest.h"
#include "visualize/UnstructuredGrid.h"
#include "visualize/XMLWriter.h"
#include "visualize/VisualizeException.h"

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>

BOOST_AUTO_TEST_CASE(DefinitionOrder)
{
    NuTo::Visualize::UnstructuredGrid visu;
    visu.DefinePointData("Stuff");
    visu.DefineCellData("Stuff");
    int pointId = visu.AddPoint(Eigen::Vector3d::Zero());
    int cellId = visu.AddCell({pointId}, NuTo::eCellTypes::VERTEX);

    // add data to undefined data field
    BOOST_CHECK_THROW(visu.SetPointData(pointId, "OtherStuff", 42.), NuTo::VisualizeException);
    BOOST_CHECK_THROW(visu.SetCellData(cellId, "OtherStuff", 42.), NuTo::VisualizeException);

    // definine data after adding points/cells
    BOOST_CHECK_THROW(visu.DefinePointData("OtherStuff"), NuTo::VisualizeException);
    BOOST_CHECK_THROW(visu.DefineCellData("OtherStuff"), NuTo::VisualizeException);
}

BOOST_AUTO_TEST_CASE(Export)
{
    NuTo::Visualize::UnstructuredGrid visu;
    visu.DefinePointData("Vector");
    visu.DefineCellData("Tensor");
    std::vector<int> pointIds;
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(0., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(1., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(1., 1., 0.)));

    for (auto pointId : pointIds)
        visu.SetPointData(pointId, "Vector", Eigen::Vector3d::Random());
    
    int cellId = visu.AddCell(pointIds, NuTo::eCellTypes::TRIANGLE);
    Eigen::VectorXd voigt(6);
    voigt[0] = 11;
    voigt[1] = 22;
    voigt[2] = 33;
    voigt[3] = 23;
    voigt[4] = 13;
    voigt[5] = 12;
    visu.SetCellData(cellId, "Tensor", voigt); 

    auto file = "VisualizeUnstructuredGridTest.vtu";
    NuTo::Visualize::XMLWriter::Export(file, visu, true);
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file);
    reader->Update();

    vtkUnstructuredGrid& grid = *reader->GetOutput();
    BOOST_CHECK_EQUAL(grid.GetNumberOfPoints(), 3);
    BOOST_CHECK_EQUAL(grid.GetNumberOfCells(), 1);

    vtkDataArray& dataArray = *grid.GetCellData()->GetArray(0);
    BOOST_CHECK_EQUAL(dataArray.GetName(), "Tensor");
    BOOST_CHECK_EQUAL(dataArray.GetNumberOfTuples(), 1);
    BOOST_CHECK_EQUAL(dataArray.GetNumberOfComponents(), 6);

    std::vector<double> data(dataArray.GetTuple(0), dataArray.GetTuple(0) + dataArray.GetNumberOfComponents());
    std::vector<double> vtkVoigt = {11, 22, 33, 12, 23, 13};
    BoostUnitTest::CheckVector(data, vtkVoigt, 6);
}
