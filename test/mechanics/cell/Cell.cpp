#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/Cell.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"


class LinearElasticLaw
{
public:
    LinearElasticLaw(double rE, double rNu)
        : mE(rE)
        , mNu(rNu)
    {
    }

    Eigen::Vector3d Stress(const Eigen::VectorXd& rStrain) const
    {
        return C() * rStrain;
    }

    Eigen::Matrix3d C() const
    {

        double factor = mE / (1.0 - (mNu * mNu));
        double C11    = factor;
        double C12    = factor * mNu;
        double C33    = factor * 0.5 * (1.0 - mNu);

        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();

        C(0, 0) = C11;
        C(1, 0) = C12;

        C(0, 1) = C12;
        C(1, 1) = C11;

        C(2, 2) = C33;

        return C;
    }

    double mE;
    double mNu;
};

class MockIntegrand : public NuTo::Integrand
{
public:
    MockIntegrand(std::vector<NuTo::DofType> rDofTypes, const LinearElasticLaw& rLaw)
        : mDofTypes(rDofTypes)
        , mLaw(rLaw)
    {
    }

    std::unique_ptr<NuTo::Integrand> Clone() const override
    {
        return std::make_unique<MockIntegrand>(*this);
    }

    NuTo::DofVector Gradient(const NuTo::CellData& rCellData, const NuTo::CellIPData& rCellIPData) override
    {
        const NuTo::DofType& dof = mDofTypes[0];
        NuTo::BMatrixStrain B    = rCellIPData.GetBMatrixStrain(dof);
        NuTo::NodeValues u       = rCellData.GetNodeValues(dof);
        NuTo::DofVector gradient;

        Eigen::VectorXd strain = B * u;
        gradient[dof]          = B.transpose() * mLaw.Stress(strain);
        return gradient;
    }

private:
    std::vector<NuTo::DofType> mDofTypes;
    const LinearElasticLaw& mLaw;
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
    NuTo::ElementSimple coordinateElement({&nCoord0, &nCoord1, &nCoord2, &nCoord3}, interpolationCoordinates);

    NuTo::InterpolationQuadLinear interpolationDisplacements(2);
    NuTo::NodeSimple nDispl0(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl1(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl2(Eigen::Vector2d({0, 0}));
    NuTo::NodeSimple nDispl3(Eigen::Vector2d({0, 0}));
    NuTo::ElementSimple displacementElement({&nDispl0, &nDispl1, &nDispl2, &nDispl3}, interpolationDisplacements);

    NuTo::DofType dofDispl("Displacements", 2, 0);
    NuTo::DofContainer<NuTo::ElementSimple*> elements;
    elements[dofDispl] = &displacementElement;

    fakeit::Mock<NuTo::IntegrationTypeBase> intType;
    constexpr double a = 0.577350269189626;
    Method(intType, GetNumIntegrationPoints)   = 4;
    Method(intType, GetIntegrationPointWeight) = 1;
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(0)).AlwaysReturn(Eigen::Vector2d({-a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(1)).AlwaysReturn(Eigen::Vector2d({a, -a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(2)).AlwaysReturn(Eigen::Vector2d({a, a}));
    fakeit::When(Method(intType, GetLocalIntegrationPointCoordinates).Using(3)).AlwaysReturn(Eigen::Vector2d({-a, a}));

    LinearElasticLaw law(E, 0.0);
    MockIntegrand integrand({dofDispl}, law);

    NuTo::Cell cell(coordinateElement, elements, intType.get(), integrand);

    auto gradient = cell.Gradient();
    BoostUnitTest::CheckVector(gradient[dofDispl], Eigen::VectorXd::Zero(8), 8);

    const double ux = 0.4;
    nDispl1.SetValue(0, ux);
    nDispl2.SetValue(0, ux);

    double area     = ly;
    double intForce = E * ux / lx * area / 2.;

    std::cout << cell.Gradient()[dofDispl] << std::endl;
    BoostUnitTest::CheckVector(cell.Gradient()[dofDispl],
                               std::vector<double>({-intForce, 0, intForce, 0, intForce, 0, -intForce, 0}), 8);

    const double uy = 0.2;
    nDispl1.SetValue(0, 0);
    nDispl2.SetValue(0, 0);
    nDispl2.SetValue(1, uy);
    nDispl3.SetValue(1, uy);

    area     = lx;
    intForce = E * uy / ly * area / 2;
    BoostUnitTest::CheckVector(cell.Gradient()[dofDispl],
                               std::vector<double>({0, -intForce, 0, -intForce, 0, intForce, 0, intForce}), 8);
}
