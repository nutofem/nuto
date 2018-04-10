#include "BoostUnitTest.h"
#include <fakeit.hpp>

#include "nuto/math/shapes/Triangle.h"

#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"

double VolumeF(const NuTo::CellIpData&)
{
    return 1.;
}

BOOST_AUTO_TEST_CASE(CellLetsSee)
{
    const double lx = 42;
    const double ly = 12;
    // const double lz = 1; // requires something like "Section" ...
    const double E = 6174;

    NuTo::InterpolationQuadLinear interpolationCoordinates;
    NuTo::NodeCoordinates nCoord0(Eigen::Vector2d({0, 0}));
    NuTo::NodeCoordinates nCoord1(Eigen::Vector2d({lx, 0}));
    NuTo::NodeCoordinates nCoord2(Eigen::Vector2d({lx, ly}));
    NuTo::NodeCoordinates nCoord3(Eigen::Vector2d({0, ly}));
    NuTo::CoordinateElementFem coordinateElement({nCoord0, nCoord1, nCoord2, nCoord3}, interpolationCoordinates);

    NuTo::InterpolationQuadLinear interpolationDisplacements;
    NuTo::NodeSimple nDispl0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl1(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl2(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl3(Eigen::Vector2d({0, 0}));
    NuTo::DofElementFem displacementElement({nDispl0, nDispl1, nDispl2, nDispl3}, interpolationDisplacements);

    NuTo::ElementCollectionFem elements(coordinateElement);
    NuTo::DofType dofDispl("Displacements", 2);
    elements.AddDofElement(dofDispl, displacementElement);

    fakeit::Mock<NuTo::IntegrationTypeBase> intType;
    Method(intType, GetNumIntegrationPoints) = 4;
    Method(intType, GetIntegrationPointWeight) = 1;
    constexpr double a = 0.577350269189626;
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(0)).AlwaysReturn(Eigen::Vector2d({-a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(1)).AlwaysReturn(Eigen::Vector2d({a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(2)).AlwaysReturn(Eigen::Vector2d({a, a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(3)).AlwaysReturn(Eigen::Vector2d({-a, a}));

    auto quad = NuTo::Quadrilateral();
    fakeit::When(Method(intType, GetShape)).AlwaysReturn(quad);

    NuTo::Laws::LinearElastic<2> law(E, 0.0, NuTo::ePlaneState::PLANE_STRAIN);
    using namespace NuTo::Integrands;
    MomentumBalance<2> integrand({dofDispl}, law);
    // bind the functions Gradient and Hessian
    auto GradientF = [&](const NuTo::CellIpData& cellIpData) {
        return integrand.Gradient(cellIpData, /*deltaT = */ 0);

    };
    auto Hessian0F = [&](const NuTo::CellIpData& cellIpData) {
        return integrand.Hessian0(cellIpData, /*deltaT = */ 0);
    };


    NuTo::Cell cell(elements, intType.get(), 1337);
    BOOST_CHECK(cell.Id() == 1337);

    BoostUnitTest::CheckVector(cell.Integrate(GradientF)[dofDispl], Eigen::VectorXd::Zero(8), 8);

    const double ux = 0.4;
    nDispl1.SetValue(0, ux);
    nDispl2.SetValue(0, ux);

    double area = ly;
    double intForce = E * ux / lx * area / 2.;

    BoostUnitTest::CheckVector(cell.Integrate(GradientF)[dofDispl],
                               std::vector<double>({-intForce, 0, intForce, 0, intForce, 0, -intForce, 0}), 8);

    const double uy = 0.2;
    nDispl1.SetValue(0, 0);
    nDispl2.SetValue(0, 0);
    nDispl2.SetValue(1, uy);
    nDispl3.SetValue(1, uy);

    area = lx;
    intForce = E * uy / ly * area / 2;
    BoostUnitTest::CheckVector(cell.Integrate(GradientF)[dofDispl],
                               std::vector<double>({0, -intForce, 0, -intForce, 0, intForce, 0, intForce}), 8);
    {
        // check hessian0
        auto hessian = cell.Integrate(Hessian0F)(dofDispl, dofDispl);
        auto gradient = cell.Integrate(GradientF)[dofDispl];
        auto u = displacementElement.ExtractNodeValues();

        BoostUnitTest::CheckEigenMatrix(gradient, hessian * u);
    }

    BOOST_CHECK_CLOSE(cell.Integrate(VolumeF), lx * ly, 1.e-10);
}

BOOST_AUTO_TEST_CASE(CellShapeMismatch)
{
    fakeit::Mock<NuTo::ElementCollection> elemCollection;
    NuTo::Triangle triangle = NuTo::Triangle();
    Method(elemCollection, GetShape) = triangle;

    fakeit::Mock<NuTo::IntegrationTypeBase> intType;
    const NuTo::Shape& quad = NuTo::Quadrilateral();
    fakeit::When(Method(intType, GetShape)).AlwaysReturn(quad);

    BOOST_CHECK_THROW(NuTo::Cell(elemCollection.get(), intType.get(), 42), NuTo::Exception);
}
