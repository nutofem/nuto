/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */
#include <iostream>
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/integrands/MomentumBalance.h"
#include "nuto/mechanics/constitutive/LocalIsotropicDamage.h"
#include "nuto/mechanics/constitutive/damageLaws/DamageLawExponential.h"

#include "nuto/mechanics/tools/QuasistaticSolver.h"
#include "nuto/mechanics/tools/GlobalFractureEnergyIntegrator.h"
/*               3   2
 *   /|          /  /           \
 *   /|         /e0/          ---\
 *   /|        /  /           ---/
 *   /|      _/  /              /
 *           0   1
 *
 *  + interface angle & thickness
 *  - angle 90 = vertical
 *  - (0,0) is at the middle of the structure
 *
 */

void CheckFractureEnergy2D(int angleDegree, double interfaceThickness)
{
    const double ly2 = 2.; // half of Length_y

    const double angleRad = M_PI / 180. * angleDegree;

    const double projectedThickness = interfaceThickness / std::sin(angleRad);
    const double xInterfaceOffset = ly2 / std::tan(angleRad);
    const double interfaceLength = 2 * ly2 / std::sin(angleRad);

    std::cout << "projectedThickness " << projectedThickness << '\n';
    std::cout << "xInterfaceOffset " << xInterfaceOffset << '\n';
    std::cout << "interfaceLength " << interfaceLength << '\n';

    // lower nodes
    NuTo::CoordinateNode n0(Eigen::Vector2d({-xInterfaceOffset - projectedThickness / 2., -ly2}));
    NuTo::CoordinateNode n1(Eigen::Vector2d({-xInterfaceOffset + projectedThickness / 2., -ly2}));
    // upper nodes
    NuTo::CoordinateNode n2(Eigen::Vector2d({+xInterfaceOffset + projectedThickness / 2., ly2}));
    NuTo::CoordinateNode n3(Eigen::Vector2d({+xInterfaceOffset - projectedThickness / 2., ly2}));

    NuTo::DofNode nd0(Eigen::Vector2d::Zero());
    NuTo::DofNode nd1(Eigen::Vector2d::Zero());
    NuTo::DofNode nd2(Eigen::Vector2d::Zero());
    NuTo::DofNode nd3(Eigen::Vector2d::Zero());
    NuTo::InterpolationQuadLinear interpolation;

    NuTo::CoordinateElementFem cElm{{n0, n1, n2, n3}, interpolation};
    NuTo::ElementCollectionFem element(cElm);
    NuTo::DofType d("Displ", 2);
    element.AddDofElement(d, {{nd0, nd1, nd2, nd3}, interpolation});

    auto material = NuTo::Material::DefaultConcrete();
    double Gf = 0.1;
    material.nu = 0.;
    material.gf = Gf / interfaceThickness;
    material.fMin = 0;
    NuTo::Laws::LocalIsotropicDamage<2> law(material);
    NuTo::Integrands::MomentumBalance<2> momentum(d, law);

    NuTo::DofNumbering::Build({nd0, nd1, nd2, nd3}, d, {}); // numbering without constraints

    NuTo::IntegrationTypeTensorProduct<2> integration(2, NuTo::eIntegrationMethod::GAUSS);

    law.mEvolution.mKappas.setZero(1, integration.GetNumIntegrationPoints());

    NuTo::Cell cell(element, integration, 0);
    auto Gradient = [&](const auto& data) { return momentum.Gradient(data, 0); };

    int dofBC1 = nd1.GetDofNumber(0);
    int dofBC2 = nd2.GetDofNumber(0);

    int numLoadSteps = 200;
    double bcEnd = 0.4;

    Eigen::VectorXd displ(numLoadSteps + 1);
    Eigen::VectorXd force(numLoadSteps + 1);

    for (int i = 0; i < numLoadSteps + 1; ++i)
    {
        double bc = bcEnd * i / (numLoadSteps);
        displ[i] = bc;
        nd1.SetValue(0, bc);
        nd2.SetValue(0, bc);
        auto internalForces = cell.Integrate(Gradient);
        force[i] = (internalForces[d][dofBC1] + internalForces[d][dofBC2]);
    }

    //    std::cout << force << std::endl;

    NuTo::Tools::GlobalFractureEnergyIntegrator integrator(force, displ);
    double crackArea = interfaceLength;
    double globalFractureEnergy = integrator.IntegrateSofteningCurve(crackArea, 0.01);
    double error = std::abs(Gf - globalFractureEnergy);
    double tolerance = Gf / 10.;

    std::cout << "angle: " << angleDegree << "\t thickness: " << interfaceThickness << "\t GF: " << globalFractureEnergy
              << "\t Error: " << error << std::endl;
    if (error > tolerance)
    {
        throw;
    }
}


int main()
{
    CheckFractureEnergy2D(90, .1);
    CheckFractureEnergy2D(90, .01);
    CheckFractureEnergy2D(90, .001);

    CheckFractureEnergy2D(75, .001);

    return EXIT_SUCCESS;
}
