/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/constitutive/LocalIsotropicDamage.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"

#include "mechanics/tools/QuasistaticSolver.h"
#include "mechanics/tools/GlobalFractureEnergyIntegrator.h"

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
    NuTo::NodeSimple n0(Eigen::Vector2d({-xInterfaceOffset - projectedThickness / 2., -ly2}));
    NuTo::NodeSimple n1(Eigen::Vector2d({-xInterfaceOffset + projectedThickness / 2., -ly2}));
    // upper nodes
    NuTo::NodeSimple n2(Eigen::Vector2d({+xInterfaceOffset + projectedThickness / 2., ly2}));
    NuTo::NodeSimple n3(Eigen::Vector2d({+xInterfaceOffset - projectedThickness / 2., ly2}));

    NuTo::NodeSimple nd0(Eigen::Vector2d::Zero());
    NuTo::NodeSimple nd1(Eigen::Vector2d::Zero());
    NuTo::NodeSimple nd2(Eigen::Vector2d::Zero());
    NuTo::NodeSimple nd3(Eigen::Vector2d::Zero());
    NuTo::InterpolationQuadLinear interpolation;

    NuTo::ElementCollectionFem element({{n0, n1, n2, n3}, interpolation});
    NuTo::DofType d("Displ", 2);
    element.AddDofElement(d, {{nd0, nd1, nd2, nd3}, interpolation});

    NuTo::Laws::LinearElasticDamage<2> elasticDamage(20000., 0.);
    double k0 = 4. / 20000;
    double Gf = 0.1;
    double gf = 4. * interfaceThickness / Gf;
    NuTo::Constitutive::DamageLawExponential damageLaw(k0, gf, 1.);
    using Law = NuTo::Laws::LocalIsotropicDamage<2, NuTo::Constitutive::DamageLawExponential,
                                                 NuTo::Laws::EvolutionImplicit<2>>;
    Law law(elasticDamage, damageLaw, {{0., 10}});

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


/*
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
        s.NodeCreate(1, Eigen::Vector3d({0, -ly, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, 0, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, ly, 0}));
        s.NodeCreate(4, Eigen::Vector3d({0, -ly / 2., lz / 2.}));
        s.NodeCreate(5, Eigen::Vector3d({0, ly / 2., lz / 2.}));
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


        s.NodeCreate(0, Eigen::Vector3d({0, -ly / 2, 0}));
        s.NodeCreate(1, Eigen::Vector3d({lx, 0, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, ly / 2, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, 0, lz}));

        s.NodeCreate(4, Eigen::Vector3d({lx / 2., -ly / 4, 0}));
        s.NodeCreate(5, Eigen::Vector3d({lx / 2., ly / 4, 0}));
        s.NodeCreate(6, Eigen::Vector3d({0, 0, 0}));

        s.NodeCreate(7, Eigen::Vector3d({0, -ly / 4, lz / 2}));
        s.NodeCreate(8, Eigen::Vector3d({0, ly / 4, lz / 2}));
        s.NodeCreate(9, Eigen::Vector3d({lx / 2., 0, lz / 2}));


        s.NodeCreate(10, Eigen::Vector3d({-lx, 0, 0}));
        s.NodeCreate(11, Eigen::Vector3d({-lx / 2., -ly / 4, 0}));
        s.NodeCreate(12, Eigen::Vector3d({-lx / 2., ly / 4, 0}));
        s.NodeCreate(13, Eigen::Vector3d({-lx / 2., 0, lz / 2}));


        s.ElementCreate(1, it, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
        s.ElementCreate(2, it, {0, 10, 2, 3, 11, 12, 6, 7, 8, 13});
        //            s.ElementCreate(2, it, {1, 2, 3, 7});

        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gAggreg, 2);
    }
    auto prism = NuTo::MeshCompanion::ElementPrismsCreate(s, gMatrix, gAggreg, thickness);

    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(prism.second, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementTotalConvertToInterpolationType();
}


void CSDA3D(int order)
{

    NuTo::Structure s(3);

    constexpr double thickness = .1;
    constexpr double lx = 10;
    constexpr double ly = 2;
    constexpr double lz = 5;

    NuTo::Interpolation::eTypeOrder displInterpol =
            order == 1 ? NuTo::Interpolation::eTypeOrder::EQUIDISTANT1 : NuTo::Interpolation::eTypeOrder::EQUIDISTANT2;

    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, displInterpol);

    int it2 = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::PRISM3D);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::DISPLACEMENTS, displInterpol);

    s.NodeCreate(0, Eigen::Vector3d({-lx, 0, 0}));
    s.NodeCreate(1, Eigen::Vector3d({-thickness / 2, -ly, 0}));
    s.NodeCreate(2, Eigen::Vector3d({-thickness / 2, ly, 0}));
    s.NodeCreate(3, Eigen::Vector3d({-thickness / 2, 0, lz}));

    s.NodeCreate(4, Eigen::Vector3d({thickness / 2, -ly, 0}));
    s.NodeCreate(5, Eigen::Vector3d({thickness / 2, ly, 0}));
    s.NodeCreate(6, Eigen::Vector3d({thickness / 2, 0, lz}));

    s.NodeCreate(7, Eigen::Vector3d({lx, 0, 0}));

    s.ElementCreate(0, it, {0, 1, 2, 3});
    s.ElementCreate(1, it, {4, 5, 6, 7});

    s.ElementCreate(2, it2, {1, 2, 3, 4, 5, 6});


    using namespace NuTo::Constitutive;
    int LIN = 0;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS, 20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO, 0.0);

    int CSDA = 1;
    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS, 200.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO, 0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetDamageLaw(CSDA, DamageLawExponential::Create(4. / 200., 4. * thickness / fractureEnergy));


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
    s.GroupAddNodeRadiusRange(groupNodeFixZ, Eigen::Vector3d({0, 0, lz}), 0, 2 * thickness);
    auto groupZ = s.GroupGetGroupPtr(groupNodeFixZ)->AsGroupNode();

    using namespace NuTo::Constraint;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        Component(nodeFixXYZ, {NuTo::eDirection::X, NuTo::eDirection::Y, NuTo::eDirection::Z}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        Component(nodeFixYZ, {NuTo::eDirection::Y, NuTo::eDirection::Z}));
    double deltaD = -.5;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        Direction(*groupZ, Eigen::Vector3d::UnitZ(), RhsRamp(1, deltaD)));

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
    newmark.PostProcessing().SetResultDirectory("./CSDA3D_" + std::to_string(order), deleteDirectory);
    newmark.Solve(1);
}
*/
int main()
{

    // CSDA3D(1);
    // CSDA3D(2);

    // PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    // PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    CheckFractureEnergy2D(90, .1);
    CheckFractureEnergy2D(90, .01);
    CheckFractureEnergy2D(90, .001);

    CheckFractureEnergy2D(75, .001);

    return EXIT_SUCCESS;
}
