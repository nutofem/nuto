#include "/usr/lib/openmpi/include/mpi.h"
#include <boost/mpi.hpp>
#include <boost/filesystem.hpp>

#include "mechanics/structures/unstructured/StructureFETI.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/NewmarkFeti.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "base/ErrorEnum.h"
#include "mechanics/sections/SectionEnum.h"



using std::cout;
using std::endl;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eSectionType;

// geometry
constexpr   int         dimension                   = 2;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      youngsModulus               = 4.0e4;
constexpr   double      poissonsRatio               = 0.2;
constexpr   double      tensileStrength             = 3;
constexpr   double      compressiveStrength         = 30;
constexpr   double      fractureEnergy              = 0.1;

// integration
constexpr   double      timeStep                    = 1.0;
constexpr   double      toleranceDisp               = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -0.01;
constexpr   double      maxInterations              = 10;

const auto directionX = Eigen::Matrix<double, dimension, 1>::UnitX();
const auto directionY = Eigen::Matrix<double, dimension, 1>::UnitY();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Feti();

void SingleDomain();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    boost::mpi::environment  env(argc, argv);
    boost::mpi::communicator world;

    Feti();
    cout << "Finished Feti" << endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (world.rank() == 0)
        SingleDomain();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Feti()
{
    boost::mpi::communicator world;
    const int rank = world.rank();

    NuTo::StructureFETI structure(dimension);
    structure.SetShowTime(false);

    std::string meshFile = "feti_rectangle_tri.msh" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile,interpolationTypeId);

    // section
    int sectionId = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, thickness);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    // constraints
    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    // loads
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    // time integration
    NuTo::NewmarkFeti timeIntegration(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + std::string("/resultFeti_") + std::to_string(structure.mRank));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetMaxNumIterations         ( maxInterations            );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);
    timeIntegration.Solve(simulationTime);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SingleDomain()
{
    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);

    std::string meshFile = "feti_rectangle_tri.msh_compare";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalSetInterpolationType(interpolationTypeId);

    structure.ImportFromGmsh(meshFile);
    structure.ElementTotalSetInterpolationType(interpolationTypeId);
    structure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    structure.SetVerboseLevel(10);
    structure.Info();
    structure.ElementTotalConvertToInterpolationType();

    // section
    int sectionId = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, thickness);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    // constraints
    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    // loads
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    // time integration
    NuTo::NewmarkDirect timeIntegration(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + std::string("/resultFeti_reference"));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetMaxNumIterations         ( maxInterations            );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);
    timeIntegration.Solve(simulationTime);


}
