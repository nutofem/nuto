#include "BoostUnitTest.h"

#include "nuto/mechanics/integrands/GradientDamage.h"

#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"

using namespace NuTo;

void CheckHessian0(ElementCollectionFem& element, Integrands::GradientDamage<1>& gdm, Eigen::VectorXd d,
                   Eigen::VectorXd eeq, double kappa, double delta = 1.e-7, double tolerance = 1.e-4)
{
    BOOST_TEST_MESSAGE("Checking for d = " << d.transpose() << " eeq = " << eeq.transpose() << " kappa = " << kappa);
    gdm.mKappas(0, 0) = kappa;
    for (int i = 0; i < 3; ++i)
        element.DofElement(gdm.mDisp).GetNode(i).SetValue(0, d[i]);
    for (int i = 0; i < 3; ++i)
        element.DofElement(gdm.mEeq).GetNode(i).SetValue(0, eeq[i]);

    CellData cellData(element, 0);
    Eigen::VectorXd ip = Eigen::Vector2d::Ones() * 0.67124;
    Jacobian jac(element.CoordinateElement().ExtractNodeValues(),
                 element.CoordinateElement().GetDerivativeShapeFunctions(ip));

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

        hessianDiff(gdm.mDisp, gdm.mEeq).col(i) = diff[gdm.mDisp];
        hessianDiff(gdm.mEeq, gdm.mEeq).col(i) = diff[gdm.mEeq];
    }

    auto Check = [&](DofType d0, DofType d1, double tol) {
        BOOST_TEST_MESSAGE("" << d0.GetName() << ", " << d1.GetName());
        double absTolerance = std::max(tol, tol * hessian0(d0, d1).cwiseAbs().maxCoeff());
        BoostUnitTest::CheckEigenMatrix(hessian0(d0, d1), hessianDiff(d0, d1), absTolerance);
    };
    Check(gdm.mDisp, gdm.mDisp, tolerance);
    Check(gdm.mDisp, gdm.mEeq, tolerance);
    Check(gdm.mEeq, gdm.mDisp, tolerance * 100); // this involves d(||e||) / d (e) which is known to be less accurate
    Check(gdm.mEeq, gdm.mEeq, tolerance);
}

BOOST_AUTO_TEST_CASE(GradientDamage1D)
{
    // coordinate element
    CoordinateNode n0(0);
    CoordinateNode n1(1);
    CoordinateNode n2(2);

    InterpolationTrussLobatto interpolation(2);
    ElementCollectionFem element({{n0, n1, n2}, interpolation});

    // displacement element nodes
    DofNode nd0(0);
    DofNode nd1(0);
    DofNode nd2(0);
    DofType disp("displacements", 1);
    element.AddDofElement(disp, {{nd0, nd1, nd2}, interpolation});

    DofNode ne0(0);
    DofNode ne1(0);
    DofNode ne2(0);
    ScalarDofType eeq("eeq");
    element.AddDofElement(eeq, {{ne0, ne1, ne2}, interpolation});

    Material::Softening m;
    m.E = 20000;
    m.nu = 0.;
    m.ft = 4;
    m.fc = 40;
    m.gf = m.ft / 360;
    m.fMin = 0.04 * m.ft;
    m.c = 0.5;

    Integrands::GradientDamage<1> gdm(disp, eeq, m, Laws::eDamageApplication::UNILATERAL);
    gdm.mKappas.setZero(1, 1);

    double k0 = m.ft / m.E;
    double delta = k0 / 10000;
    double tol = 1.e-4;

    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(0, 0, 0), 0, delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(k0, k0, k0), 2 * k0, delta, tol);
    CheckHessian0(element, gdm, Eigen::Vector3d(0.0, 0.01, 0.02), Eigen::Vector3d(9 * k0, 10 * k0, 11 * k0), 9 * k0,
                  delta, tol);
}
