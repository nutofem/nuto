#include <iostream>
#include <fstream>
#include <boost/math/tools/roots.hpp>

#include "nuto/math/EigenIO.h"
#include "nuto/math/EigenCompanion.h"
#include "nuto/mechanics/integrands/GradientDamage.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
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
double GlobalFractureEnergy(TGdm& gdm, Material::Softening material, double L = 50, int nElements = 200,
                            double boundaryDisplacement = 0.2)
{
    DofType d = gdm.mDisp;
    ScalarDofType eeq = gdm.mEeq;

    GeometryMeshFem geoMesh = UnitMeshFem::Transform(UnitMeshFem::CreateLines(nElements), [&](Eigen::VectorXd x) {
        return Eigen::VectorXd::Constant(1, x[0] * L);
    });
    MeshFem mesh(geoMesh);

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

    double k0 = material.ft / material.E;
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
    double crossSection = 1.;
    double GfTotal = gfIntegrator.IntegrateSofteningCurve(crossSection, 1.e-3 * material.ft * crossSection);
    // an exception is thrown if the forces in the load-displacement-curve do _not_ drop below this 1.e-5 * ft

    double lPreDamage = 2. * nElements / L;
    double preDamageIntegral = 0;
    double deltaK = k0 / 1000.;

    // integrate from 0 to 3*k0 with increasing damage
    for (double k = 0; k <= 3 * k0; k += deltaK)
        preDamageIntegral += (1 - gdm.mDamageLaw.Damage(k)) * material.E * k * deltaK;

    // subtract the contribution from 0 to 3*k0 with constant damage of omega(3*k0)
    double omega = gdm.mDamageLaw.Damage(3 * k0);
    for (double k = 0; k <= 3 * k0; k += deltaK)
        preDamageIntegral -= (1 - omega) * material.E * k * deltaK;

    return GfTotal - lPreDamage * preDamageIntegral;
}

double FindRootWithoutDerivative(std::function<double(double)> f, double guess, int significantBits = 10,
                                 long unsigned maxIter = 20, double factor = 2)
{
    std::pair<double, double> result = boost::math::tools::bracket_and_solve_root(
            f, guess, factor, true, boost::math::tools::eps_tolerance<double>(significantBits), maxIter);
    // the root is somewhere between result.first and result.second
    return 0.5 * (result.first + result.second);
}

int main()
{
    const double GlobalFractureEnergyParameter = 0.11;
    Material::Softening material = Material::DefaultConcrete();
    material.fMin = 1.e-6;

    auto f = [&](double gf) {
        std::cout << "Calculating for local gf = " << gf << " ... ";
        material.gf = gf;
        std::cout << std::flush;
        DofType d("Displacements", 1);
        ScalarDofType eeq("NonlocalEquivalentStrains");

        NonlocalInteraction::Decreasing interaction(0.1, 5);
        using Gdm = Integrands::GradientDamage<1, NonlocalInteraction::Decreasing>;
        Gdm gdm(d, eeq, material, Laws::eDamageApplication::FULL, interaction);
        double Gf = GlobalFractureEnergy(gdm, material, 100, 200, 1);
        std::cout << "gives global Gf = " << Gf << ".\n";
        return Gf - GlobalFractureEnergyParameter;
    };

    double gfGuess = 0.05;
    double gfCalibrated = FindRootWithoutDerivative(f, gfGuess);
    std::cout << "The fracture energy parameter gf (for GF = " << GlobalFractureEnergyParameter << ") is "
              << gfCalibrated << ".\n";
}
