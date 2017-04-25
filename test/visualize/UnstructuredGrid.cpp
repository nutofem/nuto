#include "BoostUnitTest.h"
#include "visualize/UnstructuredGrid.h"
#include "visualize/VisualizeException.h"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

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
    voigt[3] = 32;
    voigt[4] = 13;
    voigt[5] = 12;
    visu.SetCellData(cellId, "Tensor", voigt); 

    auto file = "VisualizeUnstructuredGridTest.vtu";
    visu.ExportVtuDataFile(file);
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file);
    reader->Update();
}
