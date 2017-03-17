/*
 * ElementUniaxial3D.cpp
 *
 *  Created on: 20 May 2015
 *      Author: Thomas Titscher
 */

#define BOOST_TEST_MODULE LinearInterpolation
#define BOOST_TEST_DYN_LINK

#include "ElementUniaxialTest.h"

using namespace NuTo::Interpolation;

void Run(eShapeType rShapeType, eTypeOrder rTypeOrder, std::vector<int> rDiv)
{
    NuToTest::ElementUniaxialTest test;

    NuTo::Structure s(3);
    s.SetShowTime(false);

    auto meshInfo = NuTo::MeshGenerator::Grid(s, {test.lX, test.lY, test.lZ}, rDiv, rShapeType);
    s.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);
    s.ElementTotalConvertToInterpolationType();

    test.Run(s, "Uniaxial3D");
}


BOOST_AUTO_TEST_CASE(Brick)
{
    Run(eShapeType::BRICK3D, eTypeOrder::EQUIDISTANT1, {2, 2, 2});
    Run(eShapeType::BRICK3D, eTypeOrder::EQUIDISTANT2, {1, 2, 3});
    Run(eShapeType::BRICK3D, eTypeOrder::LOBATTO2, {1, 1, 1});
    Run(eShapeType::BRICK3D, eTypeOrder::LOBATTO3, {1, 1, 1});
    Run(eShapeType::BRICK3D, eTypeOrder::LOBATTO4, {1, 1, 1});
}

BOOST_AUTO_TEST_CASE(Tetrahedron)
{
    Run(eShapeType::TETRAHEDRON3D, eTypeOrder::EQUIDISTANT1, {2, 2, 2});
    Run(eShapeType::TETRAHEDRON3D, eTypeOrder::EQUIDISTANT2, {1, 2, 3});
}

BOOST_AUTO_TEST_CASE(Prism)
{
    Run(eShapeType::PRISM3D, eTypeOrder::EQUIDISTANT1, {2, 2, 2});
    Run(eShapeType::PRISM3D, eTypeOrder::EQUIDISTANT2, {1, 2, 3});
}
