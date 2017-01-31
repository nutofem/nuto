#include <stdlib.h>
#include <iostream>
#include <boost/filesystem.hpp>
#include <math.h>
#include <boost/tokenizer.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

#include "base/Logger.h"
#include <math/EigenCompanion.h>
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/EigenSolverArpack.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/VelocityVerlet.h"
#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/RungeKutta3.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/NystroemQinZhu.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/MechanicsEnums.h"

#include "base/Exception.h"

#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRVector2General.h"

// Test of explicit time integration schemes (Kunge Kutta 4, Velocity Verlet and Nystroem Qin Zhu)
// Dynamic simulation of a 1D elastic bar under shock load
int Run(NuTo::Structure& myStructure, int timeIntegrationScheme)
{
    std::string resultDir = boost::filesystem::initial_path().string() + std::string("/ResultsExplicitTimeIntegration");
    // delete result directory
    if (boost::filesystem::exists(resultDir)) // does p actually exist?
    {
        if (boost::filesystem::is_directory(resultDir)) // is p a directory?
        {
            boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    // INPUT
    // geometry parameter
    double mL = 20; // length [m]
    double mA = 25e-4; // cross section [m²]

    // material parameter
    double youngsModulus = 70600e6; // Youngs modulus [N/m2] aluminium
    double poissonRatio = 0.; // Poisson ratio
    double density = 2824; // density [kg/m3]

    double TSchlag = 0.001; // time of shock
    double mFreq = 1; // in this load case only because of time step computation
    double mAmplS = -1e6;

    // mesh parameter
    int numElements = 100; // number of elements
    double timeStep = 1e-6;

    // computation simulation time
    double waveSpeedL =
            sqrt(youngsModulus * (1. - poissonRatio) / (density * (1. - 2. * poissonRatio) * (1. + poissonRatio)));
    double simulationTime(1. * mL / waveSpeedL);
    double timePeriod(1 / mFreq);

    // set number of time derivatives to 2 (nodes have disp, vel and accelerations)
    myStructure.SetNumTimeDerivatives(2);

    // set interpolation order
    int myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    // create nodes
    int theNode(0);
    int numNodesX = 1 * numElements + 1;
    double deltaX = mL / numElements;

    for (int countX = 0; countX < numNodesX; countX++)
    {
        Eigen::VectorXd coordinates(1);
        coordinates[0] = countX * deltaX;
        myStructure.NodeCreate(theNode, coordinates);
        theNode++;
    }

    // create elements
    std::vector<int> nodes(2);
    for (int countX = 0; countX < numElements; countX++)
    {
        nodes[0] = countX;
        nodes[1] = (countX + 1);
        myStructure.ElementCreate(myInterpolationType, nodes);
    }

    // set interpolation type to nodes
    myStructure.ElementTotalConvertToInterpolationType();

    // create section
    int mySection = myStructure.SectionCreate("TRUSS");
    myStructure.SectionSetArea(mySection, mA);
    myStructure.ElementTotalSetSection(mySection);

    // create constitutive law
    int myMaterial = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(
            myMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMaterial, NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);
    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);

    // create NodeGroups
    // locate Load: Node origin
    int grpNodes_Load = myStructure.GroupCreate("Nodes");
    double direction = 0; // X
    double min = 0. - 0.01 * (mL / numElements);
    double max = 0. + 0.01 * (mL / numElements);
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Load, direction, min, max);
    std::vector<int> Group1 = myStructure.GroupGetMemberIds(grpNodes_Load);
    std::cout << "nodes in Group Load\n";
    for (int i : Group1) std::cout << i << '\t';
    std::cout << std::endl;


    // locate results: Node x=mL
    int grpNodes_Disp = myStructure.GroupCreate("Nodes");
    direction = 0; // X
    min = mL - 0.01 * (mL / numElements);
    max = mL + 0.01 * (mL / numElements);
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Disp, direction, min, max);
    std::vector<int> Group2 = myStructure.GroupGetMemberIds(grpNodes_Disp);
    for (int i : Group2) std::cout << i << '\t';
    std::cout << std::endl;


    myStructure.CalculateMaximumIndependentSets();

    // Info
    myStructure.SetShowTime(false);
    myStructure.Info();

    // chose TimeIntegration
    NuTo::TimeIntegrationBase* myIntegrationScheme(0);
    std::cout << "timeIntegrationScheme " << timeIntegrationScheme << std::endl;
    switch (timeIntegrationScheme)
    {
    case 0:
    {
        NuTo::RungeKutta4* tmp = new NuTo::RungeKutta4(&myStructure);
        myIntegrationScheme = tmp;
        break;
    }
    case 1:
    {
        NuTo::VelocityVerlet* tmp = new NuTo::VelocityVerlet(&myStructure);
        myIntegrationScheme = tmp;
        break;
    }
    case 2:
    {
        NuTo::NystroemQinZhu* tmp = new NuTo::NystroemQinZhu(&myStructure);
        myIntegrationScheme = tmp;
        break;
    }
    }

    // time step (explicit method)
    double minTimeStepAccuracy(timePeriod / 100); // to cover the external wave
    // calculate critical time step
    double criticalTimeStep(minTimeStepAccuracy);
    if (myIntegrationScheme->HasCriticalTimeStep())
    {
        double criticalTimeStepMethod = myIntegrationScheme->CalculateCriticalTimeStep();
        if (criticalTimeStepMethod < minTimeStepAccuracy)
        {
            criticalTimeStep = criticalTimeStepMethod;
        }
    }
    if (timeStep == 0 || timeStep > criticalTimeStep)
    {
        timeStep = criticalTimeStep;
    }
    myIntegrationScheme->SetTimeStep(timeStep);

    // set unit load
    myStructure.SetNumLoadCases(1);
    Eigen::VectorXd directionL(1);
    directionL[0] = 1;
    myStructure.LoadCreateNodeGroupForce(0, grpNodes_Load, directionL, 1);

    // set load
    Eigen::Matrix<double, 4, 2> forceRHS;
    forceRHS(0, 0) = 0;
    forceRHS(0, 1) = mAmplS;
    forceRHS(1, 0) = TSchlag;
    forceRHS(1, 1) = mAmplS;
    forceRHS(2, 0) = TSchlag + timeStep;
    forceRHS(2, 1) = 0.0;
    forceRHS(3, 0) = simulationTime + 1;
    forceRHS(3, 1) = 0.0;

    // apply load
    myIntegrationScheme->SetTimeDependentLoadCase(0, forceRHS);

    // fixed displacement at origin
    Eigen::VectorXd directionD(1);
    directionD(0) = 1;
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Disp, directionD, 0.0);

    // set output
    myStructure.SetShowTime(false);
    myStructure.SetNumProcessors(1);

    myIntegrationScheme->AddResultTime("Time");
    std::vector<int> GroupNodesLoad = myStructure.GroupGetMemberIds(grpNodes_Load);
    int RightNode = GroupNodesLoad[0];
    myIntegrationScheme->AddResultNodeDisplacements("DisplacementsNodeRight", RightNode);
    int plotElement = myStructure.GetNumElements()/2;
    myIntegrationScheme->AddResultElementIpData("StressCenterElement",plotElement,NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    // only plot at every 5%
    myIntegrationScheme->SetMinTimeStepPlot(simulationTime * 0.05);

    // set result directory
    bool deleteDirectory(false);
    myIntegrationScheme->SetResultDirectory(resultDir, deleteDirectory);

    myStructure.NodeBuildGlobalDofs();

    // solve (perform explicit time integration)
    myIntegrationScheme->Solve(simulationTime);

    // analytical solution
    double uanal = mAmplS * waveSpeedL * TSchlag / (youngsModulus * mA);

    // read in result file
    boost::filesystem::path resultFile = resultDir;
    resultFile /= std::string("DisplacementsNodeRight.dat");

    Eigen::MatrixXd result = NuTo::EigenCompanion::ReadFromFile(resultFile.string());
    int EndR = result.rows();
    std::cout << "difference " << fabs(result(EndR - 1) - uanal) << "\n";

    if (fabs(result(EndR - 1) - uanal) > 1e-4)
    {
        std::cout << "difference " << fabs(result(EndR - 1) - uanal) << "\n";
        std::cout << "[ExplicitTimeIntegration] result is not correct." << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        return EXIT_SUCCESS;
    }
}

int main()
{
    try
    {
        {
            NuTo::Structure myStructure(1);
            int error1 = Run(myStructure, 0);
            return error1;
        }
        {
            NuTo::Structure myStructure(1);
            int error2 = Run(myStructure, 1);
            return error2;
        }
        {
            NuTo::Structure myStructure(1);
            int error3 = Run(myStructure, 2);
            return error3;
        }
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Error executing ExplicitTimeIntegration " << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "Error executing ExplicitTimeIntegration " << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
