#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "nuto/mechanics/cell/CellData.h"

BOOST_AUTO_TEST_CASE(CacheNodeValues)
{
    constexpr int cellID = 1337;
    fakeit::Mock<NuTo::ElementInterface> mockElement;
    Method(mockElement, ExtractNodeValues).Using(0) = Eigen::Vector2d({42, 6174});
    Method(mockElement, ExtractNodeValues).Using(1) = Eigen::Vector2d({.42, .6174});

    fakeit::Mock<NuTo::ElementCollection> elements;
    Method(elements, DofElement) = mockElement.get();

    NuTo::CellData cell(elements.get(), cellID);
    BOOST_CHECK(cell.GetCellId() == cellID);

    NuTo::DofType dof("dof", 1);


    constexpr int numRuns = 10;
    for (int iRun = 0; iRun < numRuns; ++iRun)
    {
        auto nodeValues0 = cell.GetNodeValues(dof);
        BoostUnitTest::CheckEigenMatrix(nodeValues0, Eigen::Vector2d(42, 6174));
        auto nodeValues1 = cell.GetNodeValues(dof, 1);
        BoostUnitTest::CheckEigenMatrix(nodeValues1, Eigen::Vector2d(.42, .6174));
    }
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockElement, ExtractNodeValues).Using(0)).Exactly(1));
    BOOST_CHECK_NO_THROW(fakeit::Verify(Method(mockElement, ExtractNodeValues).Using(1)).Exactly(1));
}
