#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/Cell.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"

class MockIntegrand : public NuTo::Integrand
{
public:
    std::unique_ptr<NuTo::Integrand> Clone() const override
    {
        return std::make_unique<MockIntegrand>(*this);
    }
    NuTo::DofContainer<Eigen::VectorXd> Gradient(const NuTo::CellData&, const NuTo::CellIPData&) override
    {
        NuTo::DofContainer<Eigen::VectorXd> gradient;
        return gradient;
    }
};

BOOST_AUTO_TEST_CASE(CellLetsSee)
{
    fakeit::Mock<NuTo::Interpolation> mockInterpolation;

    NuTo::NodeSimple n0(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1(Eigen::Vector2d({2, 1}));
    NuTo::NodeSimple n2(Eigen::Vector2d({1, 4}));

    NuTo::ElementSimple coordinateElement({&n0, &n1, &n2}, mockInterpolation.get());

    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    NuTo::IntegrationType2D3NGauss3Ip integrationType;
    MockIntegrand integrand;

    NuTo::Cell cell(coordinateElement, elements, integrationType, integrand);
}
