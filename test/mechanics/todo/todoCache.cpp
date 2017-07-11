#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/InterpolationSimple.h"
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
    fakeit::Mock<NuTo::ElementSimple> element;
    Method(element, ExtractNodeValues) = Eigen::Vector2d({42, 6174});
    NuTo::DofType dof("dof", 1, 0);
    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    elements[dof] = &element.get();

    NuTo::CellData cell(elements);

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        cell.GetNodeValues(dof);
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(element, ExtractNodeValues)).Exactly(1));
}
