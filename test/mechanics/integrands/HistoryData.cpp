#include "BoostUnitTest.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementCollection.h"

#include <set>
#include <functional>

using namespace NuTo;

constexpr int cellId = 354;

class CustomIntegrand
{
    std::set<int> mIpIds;
    void CheckIds(CellIds ids)
    {
        int cId = ids.cellId;
        int ipId = ids.ipId;
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

    DofVector<double> Vector(const CellIpData& cellIpData)
    {
        CheckIds(cellIpData.Ids());
        return DofVector<double>();
    }

    DofMatrix<double> Matrix(const CellIpData& cellIpData)
    {
        CheckIds(cellIpData.Ids());
        return DofMatrix<double>();
    }
};

BOOST_AUTO_TEST_CASE(Pass_Data_To_Integrand)
{
    InterpolationTrussLinear interpolation;
    NodeCoordinates n0(0);
    NodeCoordinates n1(42);
    CoordinateElementFem coordinateElement({n0, n1}, interpolation);
    ElementCollectionFem element(coordinateElement);

    CustomIntegrand integrand;
    using namespace std::placeholders;
    auto vectorF = std::bind(&CustomIntegrand::Vector, integrand, _1);
    auto matrixF = std::bind(&CustomIntegrand::Matrix, integrand, _1);

    IntegrationTypeTensorProduct<1> integrationType(2, eIntegrationMethod::GAUSS);
    Cell cell(element, integrationType, cellId);
    cell.Integrate(vectorF);

    integrand.ClearIpIds();
    cell.Integrate(matrixF);

    integrand.ClearIpIds();
    cell.Apply(matrixF);
}
