#include <iostream>
#include <fstream>

#include "nuto/math/EigenIO.h"
#include "nuto/math/EigenCompanion.h"
#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/MeshFemDofConvert.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLobatto.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/tools/CellStorage.h"
#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/AdaptiveSolve.h"
#include "nuto/mechanics/tools/GlobalFractureEnergyIntegrator.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"

using namespace NuTo;

/* MATERIAL */
double E = 30000;
double nu = 0.2;
double ft = 4;
double fc = 40;
double c = 1.00;
double k0 = ft / E;

//! Solves a 1D tensile test. The load displacement curve is integrated to obtain the global fracture energy. The
//! localization is triggered by predamaging two elements in the middle of the structure. The fracture energy
//! contribution of this imperfection is subtracted.
//! @tparam TGdm gradient damage model type
//! @param gdm gradient damage model integrand
//! @param L length of the truss structure
//! @param nElements number of elements
//! @param boundaryDisplacement boundary displacement
//! @remark This method may fail due to various reasons:
//!   - boundaryDisplacement not big enough (structure does not fully unload)
//!   - combination of c, gf (local fracture energy parameter) and L causes a snap-back
//!   - ???
template <typename TGdm>
double GlobalFractureEnergy(TGdm& gdm, double L = 50, int nElements = 200, double boundaryDisplacement = 0.1)
{
    DofType d = gdm.mDisp;
    ScalarDofType eeq = gdm.mEeq;

    MeshFem mesh = UnitMeshFem::Transform(UnitMeshFem::CreateLines(nElements),
                                          [&](Eigen::VectorXd x) { return Eigen::VectorXd::Constant(1, x[0] * L); });

    InterpolationTrussLobatto interpolationD(2);
    AddDofInterpolation(&mesh, d, interpolationD);
    AddDofInterpolation(&mesh, eeq);

    Constraint::Constraints constraints;
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d), {eDirection::X}));
    constraints.Add(d, Constraint::Component(mesh.NodesAtAxis(eDirection::X, d, L), {eDirection::X},
                                             Constraint::RhsRamp(1, boundaryDisplacement)));

    IntegrationTypeTensorProduct<1> integration(3, eIntegrationMethod::GAUSS);
    const int nIp = integration.GetNumIntegrationPoints();
    CellStorage cellStorage;
    auto cells = cellStorage.AddCells(mesh.ElementsTotal(), integration);

    gdm.mKappas.setZero(cells.Size(), nIp);
    gdm.mKappas.row(cells.Size() / 2) = Eigen::VectorXd::Constant(nIp, 3 * k0);
    gdm.mKappas.row(cells.Size() / 2 + 1) = Eigen::VectorXd::Constant(nIp, 3 * k0);

    TimeDependentProblem equations(&mesh);
    equations.AddGradientFunction(cells, TimeDependentProblem::Bind(gdm, &TGdm::Gradient));
    equations.AddHessian0Function(cells, TimeDependentProblem::Bind(gdm, &TGdm::Hessian0));
    equations.AddUpdateFunction(cells, TimeDependentProblem::Bind(gdm, &TGdm::Update));

    QuasistaticSolver problem(equations, {d, eeq});
    problem.SetQuiet();
    problem.mTolerance = 1.e-6;
    problem.SetConstraints(constraints);

    int dofLeft = mesh.NodeAtCoordinate(EigenCompanion::ToEigen(L), d).GetDofNumber(0);

    auto loadDispFileName = "CalibrationLD.dat";
    std::ofstream loadDisplacement(loadDispFileName);
    auto doStep = [&](double t) { return problem.DoStep(t, "MumpsLU"); };
    auto postProcessF = [&](double) { problem.WriteTimeDofResidual(loadDisplacement, d, {dofLeft}); };

    AdaptiveSolve adaptiveSolve(doStep, postProcessF);
    adaptiveSolve.SetQuiet();
    adaptiveSolve.dt = 0.01;
    adaptiveSolve.dtMin = 1.e-10;
    adaptiveSolve.dtMax = 0.01;
    adaptiveSolve.Solve(1.);

    loadDisplacement.close();
    Eigen::VectorXd loads = EigenIO::ReadFromFile(loadDispFileName).col(2);
    Eigen::VectorXd disps = EigenIO::ReadFromFile(loadDispFileName).col(1);

    Tools::GlobalFractureEnergyIntegrator gfIntegrator(loads, disps);
    double GfTotal = gfIntegrator.IntegrateSofteningCurve(1., 0.1);

    double lPreDamage = 2. * nElements / L;
    double preDamageIntegral = 0;
    double deltaK = k0 / 1000.;

    // integrate from 0 to 3*k0 with increasing damage
    for (double k = 0; k <= 3 * k0; k += deltaK)
        preDamageIntegral += (1 - gdm.mDamageLaw.Damage(k)) * E * k * deltaK;

    // subtract the contribution from 0 to 3*k0 with constant damage of omega(3*k0)
    double omega = gdm.mDamageLaw.Damage(3 * k0);
    for (double k = 0; k <= 3 * k0; k += deltaK)
        preDamageIntegral -= (1 - omega) * E * k * deltaK;

    return GfTotal - lPreDamage * preDamageIntegral;
}

//! Finds the root of @fFunction by using the secant method (no derivatives required)
//! @param fFunction function to find the root for
//! @param x0 first guess
//! @param x1 second guess
//! @param tolerance algorithm stopts if abs(fFunction) < tolerance
//! @param maxIterations algorithm stops after this number of steps
double SecantMethod(std::function<double(double)> fFunction, double x0, double x1, double tolerance = 1.e-4,
                    int maxIterations = 20)
{
    std::vector<double> x;
    std::vector<double> f;

    auto Info = [&](int n) {
        std::cout << "Secant step " << n << ": x_n = " << x[n] << ", f(x_n) = " << f[n] << ".\n";
    };

    x.push_back(x0);
    f.push_back(fFunction(x0));
    Info(0);

    x.push_back(x1);
    f.push_back(fFunction(x1));
    Info(1);

    for (int n = 2; n < maxIterations; ++n)
    {
        double xn = (x[n - 2] * f[n - 1] - x[n - 1] * f[n - 2]) / (f[n - 1] - f[n - 2]);
        x.push_back(xn);
        f.push_back(fFunction(xn));
        Info(n);

        if (std::abs(f[n]) < tolerance)
            return xn;
    }
    throw std::runtime_error("No convergence!");
}

int main()
{
    const double GlobalFractureEnergyParameter = 0.1;

    auto f = [&](double gf) {

        DofType d("Displacements", 1);
        ScalarDofType eeq("NonlocalEquivalentStrains");

        Laws::LinearElasticDamage<1> elasticLaw(E, nu);
        Constitutive::DamageLawExponential dmg(k0, ft / gf, 0.9999999);
        Constitutive::ModifiedMisesStrainNorm<1> strainNorm(nu, fc / ft);

        NonlocalInteraction::Decreasing interaction(0.1, 5);
        using Gdm = Integrands::GradientDamage<1, Constitutive::DamageLawExponential, NonlocalInteraction::Decreasing>;
        Gdm gdm(d, eeq, c, elasticLaw, dmg, strainNorm, interaction);
        return GlobalFractureEnergy(gdm) - GlobalFractureEnergyParameter;
    };

    double gf_guess1 = 0.015;
    double gf_guess2 = 0.01;

    double gf = SecantMethod(f, gf_guess1, gf_guess2);
    std::cout << gf << std::endl;
}
