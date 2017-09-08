#include "BoostUnitTest.h"
#include <fakeit.hpp>

#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/ElementInterpolationFEM.h"
#include "mechanics/cell/CellData.h"

BOOST_AUTO_TEST_CASE(CacheN)
{
    fakeit::Mock<NuTo::InterpolationSimple> interpolation;
    Method(interpolation, GetShapeFunctions) = Eigen::Vector2d({42, 6174});
    Method(interpolation, GetDofDimension) = 1;
    Method(interpolation, GetNumNodes) = 2;

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        interpolation.get().GetN(Eigen::Vector2d(0, 0));
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(interpolation, GetShapeFunctions)).Exactly(1));
}

BOOST_AUTO_TEST_CASE(CacheNodeValues)
{
    NuTo::NodeSimple nc0 = NuTo::NodeSimple(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nc1 = NuTo::NodeSimple(Eigen::Vector2d({0, 1}));
    NuTo::NodeSimple nc2 = NuTo::NodeSimple(Eigen::Vector2d({1, 0}));
    fakeit::Mock<NuTo::InterpolationSimple> interpolationCoord;
    NuTo::ElementInterpolationFEM eCoord({&nc0, &nc1, &nc2}, interpolationCoord.get());

    fakeit::Mock<NuTo::ElementInterpolationBase> element;
    Method(element, ExtractNodeValues) = Eigen::Vector2d({42, 6174});
    NuTo::DofType dof("dof", 1, 0);
    NuTo::DofContainer<NuTo::ElementInterpolationBase*> elements;
    elements[dof] = &element.get();

    NuTo::PDE_Element PDE_element(eCoord,elements);
    NuTo::CellData cell(PDE_element);

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        cell.ExtractNodeValues(dof);
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(element, ExtractNodeValues)).Exactly(1));
}
