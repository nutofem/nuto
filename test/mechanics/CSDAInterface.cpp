/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */

#include <cmath>
#include <boost/filesystem/operations.hpp>
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
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

int FindLocalElementIndex(int rGlobalDof, Eigen::VectorXi rGlobalDofs)
{
    for (int i = 0; i < rGlobalDofs.rows(); ++i)
    {
        if (rGlobalDofs[i] == rGlobalDof)
            return i;
    }
    throw;
}

void CheckFractureEnergy2D(int rAngleDegree, double rInterfaceThickness)
{
    NuTo::Structure s(2);
    s.SetShowTime(false);
    s.GetLogger().SetQuiet(true);

    const double ly2 = 2.;  // half of Length_y
    const double lz = 6.;

    const double angleRad = M_PI / 180. * rAngleDegree;

    const double projectedThickness = rInterfaceThickness / std::sin(angleRad);
    const double xInterfaceOffset =  ly2 / std::tan(angleRad);
    const double interfaceLength = 2 * ly2 / std::sin(angleRad);

    s.GetLogger() << "projectedThickness " << projectedThickness << '\n';
    s.GetLogger() << "xInterfaceOffset " << xInterfaceOffset << '\n';
    s.GetLogger() << "interfaceLength " << interfaceLength << '\n';

    // lower nodes
    s.NodeCreate(0, Eigen::Vector2d({-xInterfaceOffset - projectedThickness / 2., -ly2}));
    s.NodeCreate(1, Eigen::Vector2d({-xInterfaceOffset + projectedThickness / 2., -ly2}));

    // upper nodes
    s.NodeCreate(2, Eigen::Vector2d({+xInterfaceOffset + projectedThickness / 2.,  ly2}));
    s.NodeCreate(3, Eigen::Vector2d({+xInterfaceOffset - projectedThickness / 2.,  ly2}));


    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


    Eigen::VectorXi ids(4);
    ids << 0, 1, 2, 3;
    s.ElementCreate(0, it, ids);

    s.ElementTotalConvertToInterpolationType();

    int mySection = s.SectionCreate("Plane_Stress");
    s.SectionSetThickness(mySection, lz);
    s.ElementTotalSetSection(mySection);
    s.ElementTotalConvertToInterpolationType();

    using namespace NuTo::Constitutive;
    s.ConstitutiveLawCreate(0, eConstitutiveType::LOCAL_DAMAGE_MODEL);

    constexpr double fractureEnergy         = 0.1;

    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / rInterfaceThickness);
    s.ConstitutiveLawSetDamageLaw(0, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    s.ElementSetConstitutiveLaw(0, 0);

    s.NodeBuildGlobalDofs();
    int dofBC1 = s.NodeGetNodePtr(1)->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0);
    int dofBC2 = s.NodeGetNodePtr(2)->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0);

    int localDofIndex1 = FindLocalElementIndex(dofBC1, s.ElementBuildGlobalDofsRow(0)[NuTo::Node::eDof::DISPLACEMENTS]);
    int localDofIndex2 = FindLocalElementIndex(dofBC2, s.ElementBuildGlobalDofsRow(0)[NuTo::Node::eDof::DISPLACEMENTS]);

    int numLoadSteps = 200;
    double bcEnd = 0.4;

    Eigen::VectorXd displ(numLoadSteps+1);
    Eigen::VectorXd force(numLoadSteps+1);

    auto globalDofs = s.NodeExtractDofValues(0);
    auto& globalDisplacementDofs = globalDofs.J[NuTo::Node::eDof::DISPLACEMENTS];

    for (int i = 0; i < numLoadSteps+1; ++i)
    {
        double bc = bcEnd * i / (numLoadSteps);
        displ[i] = bc;
        globalDisplacementDofs[dofBC1] = bc;
        globalDisplacementDofs[dofBC2] = bc;

        s.StructureBase::NodeMergeDofValues(globalDofs);
        auto internalForces = s.ElementBuildInternalGradient(0)[NuTo::Node::eDof::DISPLACEMENTS];
        force[i] = (internalForces[localDofIndex1] + internalForces[localDofIndex2]);
    }

//    std::cout << force << std::endl;

    NuTo::Tools::GlobalFractureEnergyIntegrator integrator(force, displ);
    double crackArea = lz * interfaceLength;
    double globalFractureEnergy = integrator.IntegrateSofteningCurve(crackArea, 0.01);
    double error = std::abs(fractureEnergy - globalFractureEnergy);
    double tolerance = fractureEnergy / 10.;

    std::cout << "angle: " << rAngleDegree << "\t thickness: " << rInterfaceThickness << "\t GF: " << globalFractureEnergy << "\t Error: " << error << std::endl;
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