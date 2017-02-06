/*
 * ElementUniaxial2D.cpp
 *
 *  Created on: 20 May 2015
 *      Author: Thomas Titscher
 */

#define BOOST_TEST_MODULE LinearInterpolation
#define BOOST_TEST_DYN_LINK

#include "ElementUniaxialTest.h"

using namespace NuTo::Interpolation;

void Run(eShapeType rShapeType, eTypeOrder rTypeOrder, std::array<int, 2> rDiv)
{
    NuToTest::ElementUniaxialTest test;

    NuTo::Structure s(2);
    s.SetShowTime(false);

    auto meshInfo = NuTo::MeshGenerator::Grid<2>(s, {test.lX, test.lY}, rDiv, rShapeType);
    s.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);
    s.ElementTotalConvertToInterpolationType();

    int section = s.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    s.SectionSetThickness(section, test.lZ);
    s.ElementTotalSetSection(section);

    test.Run(s, "Uniaxial2D");
}


BOOST_AUTO_TEST_CASE(Quad)
{
    Run(eShapeType::QUAD2D, eTypeOrder::EQUIDISTANT1, {3, 2});
    Run(eShapeType::QUAD2D, eTypeOrder::EQUIDISTANT2, {3, 2});
    Run(eShapeType::QUAD2D, eTypeOrder::LOBATTO2, {1, 1});
    Run(eShapeType::QUAD2D, eTypeOrder::LOBATTO3, {1, 1});
    Run(eShapeType::QUAD2D, eTypeOrder::LOBATTO4, {1, 1});
}

BOOST_AUTO_TEST_CASE(Triangle)
{
    Run(eShapeType::TRIANGLE2D, eTypeOrder::EQUIDISTANT1, {2, 2});
    Run(eShapeType::TRIANGLE2D, eTypeOrder::EQUIDISTANT2, {2, 2});
    Run(eShapeType::TRIANGLE2D, eTypeOrder::EQUIDISTANT3, {2, 2});
    Run(eShapeType::TRIANGLE2D, eTypeOrder::EQUIDISTANT4, {2, 2});
}
