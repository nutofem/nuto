/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */

#include <cmath>
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/tools/GlobalFractureEnergyIntegrator.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "mechanics/mesh/MeshCompanion.h"

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


    std::vector<int> ids({0, 1, 2, 3});
    s.ElementCreate(0, it, ids);

    s.ElementTotalConvertToInterpolationType();

    auto mySection = NuTo::SectionPlane::Create(lz, false);
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

        s.NodeMergeDofValues(globalDofs);
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


void CSDA2D()
{

    /*        \/
     * 3    2  7     6
     *
     *
     * 0    1  4     5
     *
     */


    constexpr double thickness2 = 0.1;
    constexpr double lx2 = 10;
    constexpr double ly = 5;
    constexpr double lz = 2;



    NuTo::Structure s(2);

    s.NodeCreate(0, Eigen::Vector2d({-lx2, 0}));
    s.NodeCreate(1, Eigen::Vector2d({-thickness2, 0}));
    s.NodeCreate(2, Eigen::Vector2d({-thickness2, ly}));
    s.NodeCreate(3, Eigen::Vector2d({-lx2, ly}));

    s.NodeCreate(4, Eigen::Vector2d({thickness2, 0}));
    s.NodeCreate(5, Eigen::Vector2d({lx2, 0}));
    s.NodeCreate(6, Eigen::Vector2d({lx2, ly}));
    s.NodeCreate(7, Eigen::Vector2d({thickness2, ly}));


    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    int it2 = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementCreate(0, it, {0, 1, 2, 3});
    s.ElementCreate(1, it, {4, 5, 6, 7});

     s.ElementCreate(2, it, {1,4,7,2});
//    s.ElementCreate(2, it2, {1, 4, 7});
//    s.ElementCreate(3, it2, {1, 7, 2});

    using namespace NuTo::Constitutive;
    int LIN = 0;
    int CSDA = 1;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO,       0.0);

    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy         = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / (thickness2*2.));
    s.ConstitutiveLawSetDamageLaw(CSDA, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    s.ElementSetConstitutiveLaw(0, LIN);
    s.ElementSetConstitutiveLaw(1, LIN);
    s.ElementSetConstitutiveLaw(2, CSDA);
//    s.ElementSetConstitutiveLaw(3, CSDA);

    auto mySection = NuTo::SectionPlane::Create(lz, true);
    s.ElementTotalSetSection(mySection);
    s.ElementTotalConvertToInterpolationType();


    const auto& nodeFixXY = s.NodeGetAtCoordinate(Eigen::Vector2d({-lx2, 0}));
    const auto& nodeFixY = s.NodeGetAtCoordinate(Eigen::Vector2d({lx2, 0}));
    const auto& nodeBC = s.NodeGetAtCoordinate(Eigen::Vector2d({thickness2, ly}));

    double deltaD = -.5;
    using namespace NuTo::Constraint;
    using NuTo::Node::eDof;
    s.Constraints().Add(eDof::DISPLACEMENTS, Component(nodeFixXY, {NuTo::eDirection::X, NuTo::eDirection::Y}));
    s.Constraints().Add(eDof::DISPLACEMENTS, Component(nodeFixY, {NuTo::eDirection::Y}));
    s.Constraints().Add(eDof::DISPLACEMENTS, Direction(nodeBC, Eigen::Vector2d::UnitY(), RhsRamp(1, deltaD)));
    
    s.NodeBuildGlobalDofs();
    std::cout << s.GetNumTotalActiveDofs() << std::endl;
    std::cout << s.GetNumTotalDependentDofs() << std::endl;

    s.AddVisualizationComponent(s.GroupGetElementsTotal(), NuTo::eVisualizeWhat::DISPLACEMENTS);

    NuTo::NewmarkDirect newmark(&s);

    s.SetShowTime(false);
    newmark.SetShowTime(false);

    newmark.SetTimeStep(0.1);
    newmark.SetMinTimeStep(0.001);
    newmark.SetMaxTimeStep(0.1);
    newmark.SetToleranceForce(1e-6);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetPerformLineSearch(true);
    newmark.SetMaxNumIterations(20);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("./CSDA2D", deleteDirectory);
    newmark.Solve(1);
}



void PrismCreate(NuTo::Interpolation::eTypeOrder rCoordinateInterpolation)
{
    constexpr double thickness = .1;
    constexpr double lx = 10;
    constexpr double ly = 2;
    constexpr double lz = 5;

    NuTo::Structure s(3);

    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, rCoordinateInterpolation);

    int gMatrix = s.GroupCreate(NuTo::eGroupId::Elements);
    int gAggreg = s.GroupCreate(NuTo::eGroupId::Elements);

    if (rCoordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
    {
        s.NodeCreate(0, Eigen::Vector3d({-lx, 0, 0}));
        s.NodeCreate(1, Eigen::Vector3d({0,-ly, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, 0, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, ly, 0}));
        s.NodeCreate(4, Eigen::Vector3d({0, -ly/2., lz/2.}));
        s.NodeCreate(5, Eigen::Vector3d({0, ly/2., lz/2.}));
        s.NodeCreate(6, Eigen::Vector3d({0, 0, lz}));
        s.NodeCreate(7, Eigen::Vector3d({lx, 0, 0}));

        s.ElementCreate(1, it, {0, 1, 2, 4});
        s.ElementCreate(2, it, {0, 2, 3, 5});
        s.ElementCreate(3, it, {0, 4, 2, 5});
        s.ElementCreate(4, it, {0, 4, 5, 6});
        s.ElementCreate(5, it, {7, 1, 2, 4});
        s.ElementCreate(6, it, {7, 2, 3, 5});
        s.ElementCreate(7, it, {7, 4, 2, 5});
        s.ElementCreate(8, it, {7, 4, 5, 6});


        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gMatrix, 2);
        s.GroupAddElement(gMatrix, 3);
        s.GroupAddElement(gMatrix, 4);

        s.GroupAddElement(gAggreg, 5);
        s.GroupAddElement(gAggreg, 6);
        s.GroupAddElement(gAggreg, 7);
        s.GroupAddElement(gAggreg, 8);

    }
    else
    {
        Eigen::Vector3d(0.5, 0.0, 0.0);
        Eigen::Vector3d(0.5, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.0, 0.5);
        Eigen::Vector3d(0.0, 0.5, 0.5);
        Eigen::Vector3d(0.5, 0.0, 0.5);



        s.NodeCreate(0, Eigen::Vector3d({0,   -ly/2, 0}));
        s.NodeCreate(1, Eigen::Vector3d({lx,      0, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0,    ly/2, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0,       0, lz}));

        s.NodeCreate(4, Eigen::Vector3d({lx/2.,-ly/4, 0}));
        s.NodeCreate(5, Eigen::Vector3d({lx/2., ly/4, 0}));
        s.NodeCreate(6, Eigen::Vector3d({0,        0, 0}));

        s.NodeCreate(7, Eigen::Vector3d({0,    -ly/4, lz/2}));
        s.NodeCreate(8, Eigen::Vector3d({0,     ly/4, lz/2}));
        s.NodeCreate(9, Eigen::Vector3d({lx/2.,   0, lz/2}));


        s.NodeCreate(10, Eigen::Vector3d({-lx,      0, 0}));
        s.NodeCreate(11, Eigen::Vector3d({-lx/2.,-ly/4, 0}));
        s.NodeCreate(12, Eigen::Vector3d({-lx/2., ly/4, 0}));
        s.NodeCreate(13, Eigen::Vector3d({-lx/2.,   0, lz/2}));



        s.ElementCreate(1, it, {0,  1, 2, 3,  4,  5, 6, 7, 8,  9});
        s.ElementCreate(2, it, {0, 10, 2, 3, 11, 12, 6, 7, 8, 13});
//            s.ElementCreate(2, it, {1, 2, 3, 7});

        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gAggreg, 2);
    }
    auto prism = NuTo::MeshCompanion::ElementPrismsCreate(s, gMatrix, gAggreg, thickness);

    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(prism.second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementTotalConvertToInterpolationType();
}



void CSDA3D()
{

    NuTo::Structure s(3);

    constexpr double thickness = .1;
    constexpr double lx = 10;
    constexpr double ly = 2;
    constexpr double lz = 5;

    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    int it2 = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::PRISM3D);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.NodeCreate(0, Eigen::Vector3d({-lx, 0, 0}));
    s.NodeCreate(1, Eigen::Vector3d({-thickness/2,-ly, 0}));
    s.NodeCreate(2, Eigen::Vector3d({-thickness/2, ly, 0}));
    s.NodeCreate(3, Eigen::Vector3d({-thickness/2, 0, lz}));

    s.NodeCreate(4, Eigen::Vector3d({thickness/2,-ly, 0}));
    s.NodeCreate(5, Eigen::Vector3d({thickness/2, ly, 0}));
    s.NodeCreate(6, Eigen::Vector3d({thickness/2, 0, lz}));

    s.NodeCreate(7, Eigen::Vector3d({lx, 0, 0}));

    s.ElementCreate(0, it, {0, 1, 2, 3});
    s.ElementCreate(1, it, {4, 5, 6, 7});

    s.ElementCreate(2, it2, {1, 2, 3, 4, 5, 6});


    using namespace NuTo::Constitutive;
    int LIN = 0;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO,       0.0);

    int CSDA = 1;
    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy         = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS,       200.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / thickness);
    s.ConstitutiveLawSetDamageLaw(CSDA, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);


    s.ElementSetConstitutiveLaw(0, LIN);
    s.ElementSetConstitutiveLaw(1, LIN);
    s.ElementSetConstitutiveLaw(2, CSDA);

    s.ElementTotalConvertToInterpolationType();

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;

    s.ElementInfo(10);
    s.NodeInfo(10);

    const auto& nodeFixXYZ = s.NodeGetAtCoordinate(Eigen::Vector3d({-lx, 0, 0}));
    const auto& nodeFixYZ = s.NodeGetAtCoordinate(Eigen::Vector3d({lx, 0, 0}));
    
    int groupNodeFixZ = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeRadiusRange(groupNodeFixZ, Eigen::Vector3d({0, 0, lz}), 0, 2*thickness);
    auto groupZ = s.GroupGetGroupPtr(groupNodeFixZ)->AsGroupNode();

    using namespace NuTo::Constraint;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, Component(nodeFixXYZ, {NuTo::eDirection::X, NuTo::eDirection::Y, NuTo::eDirection::Z}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, Component(nodeFixYZ, {NuTo::eDirection::Y, NuTo::eDirection::Z}));
    double deltaD = -.5;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, Direction(*groupZ, Eigen::Vector3d::UnitZ(), RhsRamp(1, deltaD))) ;

    s.AddVisualizationComponent(s.GroupGetElementsTotal(), NuTo::eVisualizeWhat::DISPLACEMENTS);

    s.NodeBuildGlobalDofs();
    std::cout << s.GetNumTotalActiveDofs() << std::endl;
    std::cout << s.GetNumTotalDependentDofs() << std::endl;

    NuTo::NewmarkDirect newmark(&s);

    s.SetShowTime(false);
    newmark.SetShowTime(false);

    newmark.SetTimeStep(0.1);
    newmark.SetMinTimeStep(1.e-12);
    newmark.SetMaxTimeStep(0.1);
    newmark.SetToleranceForce(1.e-06);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetPerformLineSearch(true);
    newmark.SetMaxNumIterations(100);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("./CSDA3D", deleteDirectory);
    newmark.Solve(1);
}

int main()
{

    CSDA2D();
    CSDA3D();

    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    CheckFractureEnergy2D(90, .1);
    CheckFractureEnergy2D(90, .01);
    CheckFractureEnergy2D(90, .001);

    CheckFractureEnergy2D(75, .001);

    return EXIT_SUCCESS;
}
