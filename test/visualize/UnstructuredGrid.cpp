#include <iostream>
#include "BoostUnitTest.h"
#include "visualize/UnstructuredGrid.h"
#include "base/Exception.h"

BOOST_AUTO_TEST_CASE(DefinitionOrder)
{
    NuTo::Visualize::UnstructuredGrid visu;
    visu.DefinePointData("Stuff");
    visu.DefineCellData("Stuff");
    int pointId = visu.AddPoint(Eigen::Vector3d::Zero());
    int cellId = visu.AddCell({pointId}, NuTo::eCellTypes::VERTEX);

    // add data to undefined data field
    BOOST_CHECK_THROW(visu.SetPointData(pointId, "OtherStuff", 42.), NuTo::Exception);
    BOOST_CHECK_THROW(visu.SetCellData(cellId, "OtherStuff", 42.), NuTo::Exception);

    // definine data after adding points/cells
    BOOST_CHECK_THROW(visu.DefinePointData("OtherStuff"), NuTo::Exception);
    BOOST_CHECK_THROW(visu.DefineCellData("OtherStuff"), NuTo::Exception);
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
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(0., 1., 0.)));

    for (auto pointId : pointIds)
        visu.SetPointData(pointId, "Vector", Eigen::Vector3d(11, 22, 33));

    int cellId0 = visu.AddCell({0, 1, 2}, NuTo::eCellTypes::TRIANGLE);
    int cellId1 = visu.AddCell({0, 2, 3}, NuTo::eCellTypes::TRIANGLE);
    Eigen::VectorXd voigt(6);
    voigt[0] = 11;
    voigt[1] = 22;
    voigt[2] = 33;
    voigt[3] = 23;
    voigt[4] = 13;
    voigt[5] = 12;
    visu.SetCellData(cellId0, "Tensor", voigt);
    visu.SetCellData(cellId1, "Tensor", voigt);

    auto file = "VisualizeUnstructuredGridTestBinary.vtu";
    visu.ExportVtuDataFile(file, true);

    auto asciiFile = "VisualizeUnstructuredGridTestAscii.vtu";
    visu.ExportVtuDataFile(asciiFile, false);
}
