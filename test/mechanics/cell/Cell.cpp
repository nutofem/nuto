#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/Cell.h"
#include "mechanics/elements/ElementFem.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/integrands/MomentumBalance.h"


struct Volume : NuTo::ScalarOperation
{
    double operator()(NuTo::Integrand::Base&, const NuTo::CellData&, const NuTo::CellIpData&) const override
    {
        return 1.;
    }
};

BOOST_AUTO_TEST_CASE(CellLetsSee)
{
    const double lx = 42;
    const double ly = 12;
    // const double lz = 1; // requires something like "Section" ...
    const double E = 6174;

    NuTo::InterpolationQuadLinear interpolationCoordinates(2);
    NuTo::NodeSimple nCoord0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nCoord1(Eigen::Vector2d({lx, 0}));
    NuTo::NodeSimple nCoord2(Eigen::Vector2d({lx, ly}));
    NuTo::NodeSimple nCoord3(Eigen::Vector2d({0, ly}));
    NuTo::ElementFem coordinateElement({&nCoord0, &nCoord1, &nCoord2, &nCoord3}, interpolationCoordinates);

    NuTo::InterpolationQuadLinear interpolationDisplacements(2);
    NuTo::NodeSimple nDispl0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl1(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl2(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl3(Eigen::Vector2d({0, 0}));
    NuTo::ElementFem displacementElement({&nDispl0, &nDispl1, &nDispl2, &nDispl3}, interpolationDisplacements);

    NuTo::DofType dofDispl("Displacements", 2, 0);
    NuTo::DofContainer<const NuTo::ElementInterface*> elements;
    elements[dofDispl] = &displacementElement;

    NuTo::Element element(coordinateElement, elements);

    fakeit::Mock<NuTo::IntegrationTypeBase> intType;
    Method(intType, GetNumIntegrationPoints) = 4;
    Method(intType, GetIntegrationPointWeight) = 1;
    constexpr double a = 0.577350269189626;
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(0)).AlwaysReturn(Eigen::Vector2d({-a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(1)).AlwaysReturn(Eigen::Vector2d({a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(2)).AlwaysReturn(Eigen::Vector2d({a, a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(3)).AlwaysReturn(Eigen::Vector2d({-a, a}));

    NuTo::Laws::LinearElastic<2> law(E, 0.0, NuTo::ePlaneState::PLANE_STRAIN);
    using namespace NuTo::Integrand;
    TimeDependent::MomentumBalance<2> integrand({dofDispl}, law);

    NuTo::Cell cell(element, intType.get(), integrand);

    BoostUnitTest::CheckVector(cell(TimeDependent::Gradient())[dofDispl], Eigen::VectorXd::Zero(8), 8);

    const double ux = 0.4;
    nDispl1.SetValue(0, ux);
    nDispl2.SetValue(0, ux);

    double area = ly;
    double intForce = E * ux / lx * area / 2.;

    BoostUnitTest::CheckVector(cell(TimeDependent::Gradient())[dofDispl],
                               std::vector<double>({-intForce, 0, intForce, 0, intForce, 0, -intForce, 0}), 8);

    const double uy = 0.2;
    nDispl1.SetValue(0, 0);
    nDispl2.SetValue(0, 0);
    nDispl2.SetValue(1, uy);
    nDispl3.SetValue(1, uy);

    area = lx;
    intForce = E * uy / ly * area / 2;
    BoostUnitTest::CheckVector(cell(TimeDependent::Gradient())[dofDispl],
                               std::vector<double>({0, -intForce, 0, -intForce, 0, intForce, 0, intForce}), 8);
    {
        // check hessian0
        auto hessian = cell(TimeDependent::Hessian0())(dofDispl, dofDispl);
        auto gradient = cell(TimeDependent::Gradient())[dofDispl];
        auto u = displacementElement.ExtractNodeValues();

        BoostUnitTest::CheckEigenMatrix(gradient, hessian * u);
    }

    BOOST_CHECK_CLOSE(cell(Volume()), lx*ly, 1.e-10);
}
