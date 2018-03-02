#include "BoostUnitTest.h"

#include "mechanics/integrands/GradientDamage.h"

#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/interpolation/InterpolationTrussLobatto.h"

using namespace NuTo;

void CheckHessian0(ElementCollectionFem& element,
                   Integrands::GradientDamage<1, Constitutive::DamageLawExponential>& gdm, Eigen::VectorXd d,
                   Eigen::VectorXd eeq, double kappa, double delta = 1.e-7, double tolerance = 1.e-4)
{
    gdm.mKappas(0, 0) = kappa;
    for (int i = 0; i < 3; ++i)
        element.DofElement(gdm.mDisp).GetNode(i).SetValue(0, d[i]);
    for (int i = 0; i < 3; ++i)
        element.DofElement(gdm.mEeq).GetNode(i).SetValue(0, eeq[i]);

    CellData cellData(element, 0);
    Eigen::VectorXd ip = Eigen::VectorXd::Ones(0) * 0.3333333333333;
    Jacobian jac(element.CoordinateElement().ExtractNodeValues(),
                 element.CoordinateElement().GetDerivativeShapeFunctions(ip), 1);

    CellIpData cipd(cellData, jac, ip, 0);

    auto hessian0 = gdm.Hessian0(cipd);
    auto gradient = gdm.Gradient(cipd);

    auto hessianDiff = hessian0;
    hessianDiff *= 0;

    BOOST_TEST_MESSAGE("D Kappa D Eeq " << gdm.DkappaDeeq(cipd));

    for (int i = 0; i < 3; ++i)
    {
        CellData cellDataI(element, 0);
        CellIpData cipdI(cellDataI, jac, ip, 0);
        element.DofElement(gdm.mDisp).GetNode(i).SetValue(0, d[i] + delta);

        auto diff = gdm.Gradient(cipdI);
        diff.AddScaled(gradient, -1);
        diff *= 1. / delta;

        element.DofElement(gdm.mDisp).GetNode(i).SetValue(0, d[i]);

        // BOOST_TEST_MESSAGE("Diff\n" << diff);

        hessianDiff(gdm.mDisp, gdm.mDisp).col(i) = diff[gdm.mDisp];
        hessianDiff(gdm.mEeq, gdm.mDisp).col(i) = diff[gdm.mEeq];
    }

    for (int i = 0; i < 3; ++i)
    {
        CellData cellDataI(element, 0);
        CellIpData cipdI(cellDataI, jac, ip, 0);
        element.DofElement(gdm.mEeq).GetNode(i).SetValue(0, eeq[i] + delta);

        auto diff = gdm.Gradient(cipdI);
        diff.AddScaled(gradient, -1);
        diff *= 1. / delta;

        element.DofElement(gdm.mEeq).GetNode(i).SetValue(0, eeq[i]);

        // BOOST_TEST_MESSAGE("Diff\n" << diff);

        hessianDiff(gdm.mDisp, gdm.mEeq).col(i) = diff[gdm.mDisp];
        hessianDiff(gdm.mEeq, gdm.mEeq).col(i) = diff[gdm.mEeq];
    }
    // BOOST_TEST_MESSAGE("Hessian0\n" << hessian0);
    // BOOST_TEST_MESSAGE("Hessian0Diff\n" << hessianDiff);

    auto Check = [&](DofType d0, DofType d1) {
        BOOST_TEST_MESSAGE("" << d0.GetName() << ", " << d1.GetName());
        double absTolerance = std::max(tolerance, tolerance * hessian0(d0, d1).cwiseAbs().maxCoeff());
        BoostUnitTest::CheckEigenMatrix(hessian0(d0, d1), hessianDiff(d0, d1), absTolerance);
    };
    // Check(gdm.mDisp, gdm.mDisp);
    // Check(gdm.mDisp, gdm.mEeq);
    // Check(gdm.mEeq, gdm.mDisp);
    Check(gdm.mEeq, gdm.mEeq);
}

BOOST_AUTO_TEST_CASE(NeumannBc1Din2D)
{
    // coordinate element
    NodeSimple n0(0);
    NodeSimple n1(1);
    NodeSimple n2(2);

    InterpolationTrussLobatto interpolation(2);
    ElementCollectionFem element({{n0, n1, n2}, interpolation});

    // displacement element nodes
    NodeSimple nd0(0);
    NodeSimple nd1(0);
    NodeSimple nd2(0);
    DofType disp("displacements", 1);
    element.AddDofElement(disp, {{nd0, nd1, nd2}, interpolation});

    NodeSimple ne0(0);
    NodeSimple ne1(0);
    NodeSimple ne2(0);
    ScalarDofType eeq("eeq");
    element.AddDofElement(eeq, {{ne0, ne1, ne2}, interpolation});

    double E = 20000;
    double nu = 0.2;
    double k0 = 4. / E;
    double beta = 360;
    double alpha = 0.96;
    double c = 200;
    Laws::LinearElastic<1> elasticLaw(E, nu);
    Constitutive::ModifiedMisesStrainNorm<1> strainNorm(nu, 1);
    Constitutive::DamageLawExponential damageLaw(k0, beta, alpha);

    Integrands::GradientDamage<1, Constitutive::DamageLawExponential> gdm(disp, eeq, c, elasticLaw, damageLaw,
                                                                          strainNorm);
    gdm.mKappas.setZero(1, 1);

    double delta = k0 / 200;
    double tol = 1.e-7;

    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(0, 0, 0), 0, delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(k0, k0, k0), k0, delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(9 * k0, 10 * k0, 11 * k0), 9 * k0,
                  delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(k0, k0, k0), k0, delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(9 * k0, 10 * k0, 11 * k0), 9 * k0,
                  delta, tol);
}
