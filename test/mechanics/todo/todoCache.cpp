#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/interpolation/InterpolationSimple.h"
#include "mechanics/cell/CellData.h"

BOOST_AUTO_TEST_CASE(CacheNodeValues)
{
    fakeit::Mock<NuTo::CellInterpolationBase> element;
    Method(element, ExtractNodeValues) = Eigen::Vector2d({42, 6174});
    NuTo::DofType dof("dof", 1, 0);
    NuTo::DofContainer<NuTo::CellInterpolationBase*> elements;
    elements[dof] = &element.get();

    NuTo::CellData cell(elements);

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        cell.GetNodeValues(dof);
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(element, ExtractNodeValues)).Exactly(1));
}
