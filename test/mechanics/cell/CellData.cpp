#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellData.h"
#include "mechanics/elements/Element.h"

BOOST_AUTO_TEST_CASE(CacheNodeValues)
{
    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, ExtractNodeValues) = Eigen::Vector2d({42, 6174});
    NuTo::DofType dof("dof", 1, 0);
    NuTo::DofContainer<const NuTo::ElementInterface*> elements;
    elements[dof] = &mockElement.get();
    NuTo::Element element(mockElement.get(), elements); 

    NuTo::CellData cell(element);

    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        auto nodeValues = cell.GetNodeValues(dof);
        BoostUnitTest::CheckEigenMatrix(nodeValues, Eigen::Vector2d(42, 6174));
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockElement, ExtractNodeValues)).Exactly(1));
}
