#include "BoostUnitTest.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"

#include <set>
#include <functional>

using namespace NuTo;

constexpr int cellId = 354;

class CustomIntegrand
{
    std::set<int> mIpIds;
    void CheckIds(int cId, int ipId)
    {
        BOOST_CHECK_EQUAL(cId, cellId);

        // Each IP Id should only be passed once per cell and values have to be 0 and 1 (2 ips)
        BOOST_CHECK(ipId == 0 || ipId == 1);
        BOOST_CHECK(mIpIds.find(ipId) == mIpIds.end());
        mIpIds.insert(ipId);
    }

public:
    void ClearIpIds()
    {
        mIpIds.clear();
    }

    DofVector<double> Vector(const CellData& cellData, const CellIpData& cellIpData)
    {
        CheckIds(cellData.GetCellId(), cellIpData.GetIpId());
        return DofVector<double>();
    }

    DofMatrix<double> Matrix(const CellData& cellData, const CellIpData& cellIpData)
    {
        CheckIds(cellData.GetCellId(), cellIpData.GetIpId());
        return DofMatrix<double>();
    }
};

BOOST_AUTO_TEST_CASE(Pass_Data_To_Integrand)
{
    InterpolationTrussLinear interpolation(1);
    NodeSimple n0(0);
    NodeSimple n1(42);
    ElementFem coordinateElement({n0, n1}, interpolation);
    ElementCollectionFem element(coordinateElement);

    CustomIntegrand integrand;
    using namespace std::placeholders;
    auto vectorF = std::bind(&CustomIntegrand::Vector, integrand, _1, _2);
    auto matrixF = std::bind(&CustomIntegrand::Matrix, integrand, _1, _2);

    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);
    Cell cell(element, integrationType, cellId);
    cell.Integrate(vectorF);

    integrand.ClearIpIds();
    cell.Integrate(matrixF);

    integrand.ClearIpIds();
    cell.Apply(matrixF);
}
