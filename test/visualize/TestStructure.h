#pragma once
#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "math/shapes/Quadrilateral.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace NuTo
{
namespace Test
{

class VisualizeTestStructure
{
    fakeit::Mock<NuTo::CellInterface> cell0;
    fakeit::Mock<NuTo::CellInterface> cell1;
    Quadrilateral myQuad;

public:
    NuTo::Group<NuTo::CellInterface> Cells()
    {
        Eigen::VectorXd bla = Eigen::Vector3d::Constant(1.0);
        fakeit::When(ConstOverloadedMethod(cell0, Interpolate, Eigen::VectorXd(Eigen::VectorXd))).AlwaysReturn(bla);
        fakeit::When(ConstOverloadedMethod(cell1, Interpolate, Eigen::VectorXd(Eigen::VectorXd))).AlwaysReturn(bla);
        fakeit::When(Method(cell0, Eval)).AlwaysReturn({bla, bla, bla, bla});
        fakeit::When(Method(cell1, Eval)).AlwaysReturn({bla, bla, bla, bla});
        fakeit::When(Method(cell0, GetShape)).AlwaysReturn(myQuad);
        fakeit::When(Method(cell1, GetShape)).AlwaysReturn(myQuad);
        return {cell0.get(), cell1.get()};
    }
};

} /* Test */
} /* NuTo */

namespace pt = boost::property_tree;

struct UnstructuredGridCheck
{
public:
    static void CheckNum(std::string filename, int expectedNumPoints, int expectedNumCells)
    {
        pt::ptree tree;
        pt::read_xml(filename, tree);
        int numOfPoints = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfPoints");
        int numOfCells = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfCells");
        BOOST_CHECK_EQUAL(numOfPoints, expectedNumPoints);
        BOOST_CHECK_EQUAL(numOfCells, expectedNumCells);
    }
};
