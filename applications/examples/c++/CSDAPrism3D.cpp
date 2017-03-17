/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */

#include <boost/filesystem/operations.hpp>
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"

#include "visualize/VisualizeEnum.h"

#include "mechanics/mesh/MeshCompanion.h"


int main()
{

    NuTo::Structure s(3);
    auto ids = s.ImportFromGmsh("CSDAMesh.msh");

    int gMatrix = ids[0].first;
    int gAggreg = ids[1].first;

    using namespace NuTo::Constitutive;
    int LIN = 0;
    int CSDA = 1;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO,       0.0);

    const double thickness = 0.1;

    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy         = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS,       200.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / thickness);
    s.ConstitutiveLawSetDamageLaw(CSDA, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;

    auto prism = NuTo::MeshCompanion::ElementPrismsCreate(s, gMatrix, gAggreg, thickness);
    s.ElementTotalSetConstitutiveLaw(LIN);
    s.ElementGroupSetConstitutiveLaw(prism.first, CSDA);

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;


    s.InterpolationTypeAdd(ids[0].second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(ids[1].second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(prism.second,  NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementTotalConvertToInterpolationType();

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;

    s.ElementInfo(10);
    s.NodeInfo(10);

    int nodeFixXYZ = s.NodeGetIdAtCoordinate(Eigen::Vector3d({-5, 0, 0}), 1.e-5);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitX(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitY(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitZ(), 0);

    int nodeFixYZ = s.NodeGetIdAtCoordinate(Eigen::Vector3d({5, 0, 0}), 1.e-5);
    s.ConstraintLinearSetDisplacementNode(nodeFixYZ, Eigen::Vector3d::UnitY(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixYZ, Eigen::Vector3d::UnitZ(), 0);

    int groupNodeFixZ = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeRadiusRange(groupNodeFixZ, Eigen::Vector3d({0, 0, 2}), 0, 2*thickness);

    s.ConstraintLinearSetDisplacementNodeGroup(groupNodeFixZ, Eigen::Vector3d::UnitY(), 0);

    int BC = s.ConstraintLinearSetDisplacementNodeGroup(groupNodeFixZ, Eigen::Vector3d::UnitZ(), 0);

    s.NodeBuildGlobalDofs();
    std::cout << s.GetNumTotalActiveDofs() << std::endl;
    std::cout << s.GetNumTotalDependentDofs() << std::endl;

    s.AddVisualizationComponent(s.GroupGetElementsTotal(), NuTo::eVisualizeWhat::DISPLACEMENTS);

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
    newmark.SetResultDirectory("./CSDAPrism3DOut", deleteDirectory);

    double deltaD = .5;
    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, 1, -deltaD;

    newmark.AddTimeDependentConstraint(BC, dispRHS);
    newmark.Solve(1);

    return EXIT_SUCCESS;
}
