/*
 * ElementUniaxial1D.cpp
 *
 *  Created on: 20 May 2015
 *      Author: Thomas Titscher
 */

#define BOOST_TEST_MODULE LinearInterpolation
#define BOOST_TEST_DYN_LINK

#include "ElementUniaxialTest.h"

using namespace NuTo::Interpolation;

void Run(eShapeType rShapeType, eTypeOrder rTypeOrder, int rDiv)
{
    NuToTest::ElementUniaxialTest test;

    NuTo::Structure s(1);
    s.SetShowTime(false);

    auto meshInfo = NuTo::MeshGenerator::Grid<1>(s, {test.lX}, {rDiv}, rShapeType);
    s.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);
    s.ElementTotalConvertToInterpolationType();

    int section = s.SectionCreate(NuTo::eSectionType::TRUSS);
    s.SectionSetArea(section, test.lY * test.lZ);

    s.ElementTotalSetSection(section);

    test.Run(s, "Uniaxial1D");
}


BOOST_AUTO_TEST_CASE(Truss)
{
    Run(eShapeType::TRUSS1D, eTypeOrder::EQUIDISTANT1, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::EQUIDISTANT2, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::EQUIDISTANT3, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::EQUIDISTANT4, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::LOBATTO2, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::LOBATTO3, 3);
    Run(eShapeType::TRUSS1D, eTypeOrder::LOBATTO4, 3);
}
