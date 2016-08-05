#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/constitutive/laws/PhaseField.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

using std::cout;
using std::endl;

int main()
{
    try
    {

        constexpr unsigned int dimension = 2;

        constexpr bool performLineSearch = true;
        constexpr bool automaticTimeStepping = true;

        constexpr double youngsModulus = 2.1e5;
        constexpr double poissonsRatio = 0.3;
        constexpr double thickness = 1.0;
        constexpr double lengthScale = 1.0;
        constexpr double fractureEnergy = 2.7;
        constexpr double artificialViscosity = 0.1;

        constexpr double timeStep = 1e-2;
        constexpr double minTimeStep = 1e-5;
        constexpr double maxTimeStep = 1e-1;
        constexpr double toleranceForce = 1e-6;
        constexpr double simulationTime = 1.0;
        constexpr double load = 0.01;

        constexpr double tolerance = 1.0e-6;

        boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/resultPhaseField/");
        boost::filesystem::remove_all(resultPath);
        boost::filesystem::create_directory(resultPath);

        const NuTo::FullVector<double, dimension> directionX = NuTo::FullVector<double, dimension>::UnitX();
        const NuTo::FullVector<double, dimension> directionY = NuTo::FullVector<double, dimension>::UnitY();

        cout << "**********************************************" << endl;
        cout << "**  strucutre                               **" << endl;
        cout << "**********************************************" << endl;

        NuTo::Structure myStructure(dimension);
        myStructure.SetShowTime(false);
        myStructure.SetNumTimeDerivatives(1);

        cout << "**********************************************" << endl;
        cout << "**  integration sheme                       **" << endl;
        cout << "**********************************************" << endl;

        NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
        myIntegrationScheme.SetTimeStep(timeStep);
        myIntegrationScheme.SetMinTimeStep(minTimeStep);
        myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
        myIntegrationScheme.SetToleranceForce(toleranceForce);
        myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
        myIntegrationScheme.SetPerformLineSearch(performLineSearch);
        myIntegrationScheme.SetResultDirectory(resultPath.string(), true);

        cout << "**********************************************" << endl;
        cout << "**  section                                 **" << endl;
        cout << "**********************************************" << endl;

        int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
        myStructure.SectionSetThickness(mySection, thickness);

        cout << "**********************************************" << endl;
        cout << "**  material                                **" << endl;
        cout << "**********************************************" << endl;

        NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus,
                                                                  poissonsRatio,
                                                                  lengthScale,
                                                                  fractureEnergy,
                                                                  artificialViscosity);


        int matrixMaterial = myStructure.AddConstitutiveLaw(phaseField);

        cout << "**********************************************" << endl;
        cout << "**  geometry                                **" << endl;
        cout << "**********************************************" << endl;

        NuTo::FullVector<int, -1> nodeIds(3);
        NuTo::FullVector<double, -1> nodeCoords(2);

        nodeCoords[0] = 0.0;
        nodeCoords[1] = 0.0;
        nodeIds[0] = myStructure.NodeCreate(nodeCoords);

        nodeCoords[0] = 2.0;
        nodeCoords[1] = 0.0;
        nodeIds[1] = myStructure.NodeCreate(nodeCoords);

        nodeCoords[0] = 1.0;
        nodeCoords[1] = 1.0;
        nodeIds[2] = myStructure.NodeCreate(nodeCoords);

        int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);

        myStructure.ElementCreate(myInterpolationType, nodeIds);

        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::CRACKPHASEFIELD, NuTo::Interpolation::EQUIDISTANT1);

        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::STATICDATA);
        myStructure.InterpolationTypeInfo(myInterpolationType);

        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        myStructure.ElementTotalSetSection(mySection);
        myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

        cout << "**********************************************" << endl;
        cout << "**  bc                                      **" << endl;
        cout << "**********************************************" << endl;

        NuTo::FullVector<double, 2> center;

        // bottom left boundary
        center[0] = 0;
        center[1] = 0;
        int grpNodes_bottom_left = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_left, center, 0, tolerance);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, directionX, 0);

        // bottom right boundary
        int grpNodes_bottom_right = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeCoordinateRange(grpNodes_bottom_right, 1, -tolerance, tolerance);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_right, directionY, 0);

        cout << "**********************************************" << endl;
        cout << "**  load                                    **" << endl;
        cout << "**********************************************" << endl;

        // middle top load
        center[0] = 1;
        center[1] = 1;
        int grpNodes_load = myStructure.GroupCreate(NuTo::Groups::Nodes);

        myStructure.GroupAddNodeRadiusRange(grpNodes_load, center, 0, tolerance);

        int loadId = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_load, directionY, 0);

        cout << "**********************************************" << endl;
        cout << "**  visualization                           **" << endl;
        cout << "**********************************************" << endl;

        int groupId = myStructure.GroupGetElementsTotal();
        myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::CRACK_PHASE_FIELD);

        cout << "**********************************************" << endl;
        cout << "**  solver                                  **" << endl;
        cout << "**********************************************" << endl;

        myStructure.NodeBuildGlobalDofs();
        myStructure.CalculateMaximumIndependentSets();

        NuTo::FullMatrix<double, 2, 2> dispRHS;
        dispRHS(0, 0) = 0;
        dispRHS(1, 0) = simulationTime;
        dispRHS(0, 1) = 0;
        dispRHS(1, 1) = load;

        myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);

        myIntegrationScheme.Solve(simulationTime);

    } catch (...)
    {
        cout << "Test failed" << endl;
        return EXIT_FAILURE;
    }

    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}

