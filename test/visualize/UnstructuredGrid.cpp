#include <iostream>
#include "BoostUnitTest.h"
#include "nuto/visualize/UnstructuredGrid.h"
#include "nuto/base/Exception.h"

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

int AddLine2ndOrder(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(0., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(1., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(0.5, 0.5, 0.))); // middle Node

    visu.SetPointData(pointIds[0], "Scalar", 2.);
    visu.SetPointData(pointIds[1], "Scalar", 2.);
    visu.SetPointData(pointIds[2], "Scalar", 3.);

    return visu.AddCell(pointIds, NuTo::eCellTypes::LINE2NDORDER);
}

int AddTriangle2ndOrder(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(2., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(3., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(3., 3., 0.)));
    // 2nd order
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(2.5, 1.97, 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(3.01, 2.5, 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(2.44, 2.56, 0.)));

    for (size_t i = 0; i < 3; i++)
        visu.SetPointData(pointIds[i], "Scalar", 2.);
    for (size_t i = 3; i < 6; i++)
        visu.SetPointData(pointIds[i], "Scalar", 3.);

    return visu.AddCell(pointIds, NuTo::eCellTypes::TRIANGLE2NDORDER);
}

int AddQuad2ndOrder(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(4., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(5., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(5., 3., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(4., 3., 0.)));
    // 2nd order
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(4.5, 1.97, 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(5.01, 2.5, 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(4.5, 3.02, 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(4.02, 2.5, 0.)));

    for (size_t i = 0; i < 4; i++)
        visu.SetPointData(pointIds[i], "Scalar", 2.);
    for (size_t i = 4; i < 8; i++)
        visu.SetPointData(pointIds[i], "Scalar", 3.);

    return visu.AddCell(pointIds, NuTo::eCellTypes::QUAD2NDORDER);
}

int AddPyramid(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(6., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(7., 0., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(7., 1., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(6., 1., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(6.5, 0.5, 1.)));

    for (auto pointId : pointIds)
    {
        visu.SetPointData(pointId, "Scalar", 2.);
    }
    return visu.AddCell(pointIds, NuTo::eCellTypes::PYRAMID);
}

int AddTet2ndOrder(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(8., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9., 3., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9., 3., 1.)));
    // 2nd order
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(8.5, 2., 0.05))); // between 0,1
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9., 2.5, 0.05))); // between 1,2
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(8.5, 2.5, 0.05))); // between 2,0
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(8.5, 2.5, 0.55))); // between 0,3
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9.0, 2.5, 0.45))); // between 1,3
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(9.0, 3.0, 0.55))); // between 2,3

    for (size_t i = 0; i < 4; i++)
        visu.SetPointData(pointIds[i], "Scalar", 2.);
    for (size_t i = 4; i < 10; i++)
        visu.SetPointData(pointIds[i], "Scalar", 3.);

    return visu.AddCell(pointIds, NuTo::eCellTypes::TETRAEDER2NDORDER);
}

int AddHex2ndOrder(NuTo::Visualize::UnstructuredGrid& visu)
{
    visu.DefinePointData("Scalar");

    std::vector<int> pointIds;

    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 2., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 3., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 3., 0.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 2., 1.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 2., 1.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 3., 1.)));
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 3., 1.)));
    // 2nd order
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10.5, 2., 0.2))); // 8 - 0,1
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 2.5, 0.2))); // 9 - 1,2
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10.5, 3., 0.2))); // 10 - 2,3
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 2.5, 0.2))); // 11 - 3,0
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10.5, 2., 1.2))); // 12 - 4,5
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 2.5, 1.2))); // 13 - 5,6
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10.5, 3., 1.2))); // 14 - 6,7
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 2.5, 1.2))); // 15 - 7,4
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 2., 0.5))); // 16 - 0,4
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 2., 0.5))); // 17 - 1,13
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(11., 3., 0.5))); // 18 - 2,6
    pointIds.push_back(visu.AddPoint(Eigen::Vector3d(10., 3., 0.5))); // 19 - 3,7

    for (size_t i = 0; i < 8; i++)
        visu.SetPointData(pointIds[i], "Scalar", 2.);
    for (size_t i = 8; i < 20; i++)
        visu.SetPointData(pointIds[i], "Scalar", 3.);

    return visu.AddCell(pointIds, NuTo::eCellTypes::HEXAHEDRON2NDORDER);
}

BOOST_AUTO_TEST_CASE(Export2ndOrderElements)
{
    NuTo::Visualize::UnstructuredGrid visu;

    AddLine2ndOrder(visu);
    AddTriangle2ndOrder(visu);
    AddQuad2ndOrder(visu);
    AddTet2ndOrder(visu);
    AddHex2ndOrder(visu);
    AddPyramid(visu); // this is not a second order element but also new ...

    auto asciiFile = "Export2ndOrderElementsAscii.vtu";
    visu.ExportVtuDataFile(asciiFile, false);
}
