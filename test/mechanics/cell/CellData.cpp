#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellData.h"

BOOST_AUTO_TEST_CASE(CacheNodeValues)
{
    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, ExtractNodeValues) = Eigen::Vector2d({42, 6174});
    NuTo::DofType dof("dof", 1, 0);
    NuTo::ElementCollection elements(mockElement.get());
    elements.AddDofElement(dof, mockElement.get());

    NuTo::CellData cell(elements);

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        auto nodeValues = cell.GetNodeValues(dof);
        BoostUnitTest::CheckEigenMatrix(nodeValues, Eigen::Vector2d(42, 6174));
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockElement, ExtractNodeValues)).Exactly(1));
}
