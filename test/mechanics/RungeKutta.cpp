#include <stdlib.h>
#include <sstream>
#include "boost/filesystem.hpp"

#include "math/FullMatrix.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/RungeKutta3.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/RungeKuttaCashKarp.h"
#include "mechanics/timeIntegration/RungeKuttaDormandPrince.h"

#include "mechanics/tools/MeshGenerator.h"

/*
 *  TT:
 *  This is just a test to see if the time integration scheme compiles/runs. I highly doubt the
 *  correctness of the test setup and its results.
 */

void Run(NuTo::Structure& myStructure, NuTo::RungeKuttaBase& rTimeIntegrationScheme)
{
    std::string resultDir = boost::filesystem::initial_path().string()+std::string("/ResultsRungeKutta");
    //delete result directory
    if (boost::filesystem::exists(resultDir))// does p actually exist?
    {
        if (boost::filesystem::is_directory(resultDir))      // is p a directory?
        {
            boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    double mL = 100;
    double mH = 10;
    int numElementsX = 5;
    int numElementsY = 1;
    double density = 1;//N and mm^3


    //set number of time derivatives to 2 (nodes have disp, vel and accelerations)
    myStructure.SetNumTimeDerivatives(2);

    int myInterpolationType = myStructure.InterpolationTypeCreate("Quad2D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    //section
    double thickness(1);
    int mySection = myStructure.SectionCreate("Plane_Stress");
    myStructure.SectionSetThickness(mySection, thickness);
    myStructure.ElementTotalSetSection(mySection);

    //constitutive
    int myMatLattice = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLattice);

    NuTo::MeshGenerator::MeshRectangularPlane(
            myStructure,
            mySection,
            myMatLattice,
            myInterpolationType,
            std::array<int,2>(
                    {   numElementsX,numElementsY}),
            std::array<double,2>(
                    {   mL, mH}));

    myStructure.ElementTotalConvertToInterpolationType();

    //create node groups bottom boundary
    int nOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 2>(
                    {   0.,0.}), 1.e-6);

    //right boundary
    int grpNodes_Right = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Right, 0, mL, mL);

    int grpNodes_Left = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Right, 0, 0., 0.);

    //calculate maximum independent sets
    myStructure.CalculateMaximumIndependentSets();

    //set constraints
    //directionX
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> directionX(2, 1);
    directionX.SetValue(0, 0, 1.0);
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> directionY(2, 1);
    directionY.SetValue(1, 0, 1.0);

    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Left, directionX, 0);

    myStructure.LoadCreateNodeGroupForce(0, grpNodes_Right, directionX, 1000);

    myStructure.Info();
    NuTo::RungeKutta4 myIntegrationScheme(&myStructure);

    double simulationTime(1);

    //set a sinusoidal quarter wave
    double period(5);
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> forceRHS(101, 2);
    for (int count = 0; count < forceRHS.GetNumRows() - 1; count++)
    {
        double t = ((double) count) / ((double) forceRHS.GetNumRows() - 2.) * 0.25 * period;
        forceRHS(count, 0) = t;
        forceRHS(count, 1) = 0.1 * sin(t / period * 2. * M_PI);
        //loadRHSFactor(count,0) = 0;
    }
    forceRHS(forceRHS.GetNumRows() - 1, 0) = forceRHS(forceRHS.GetNumRows() - 2, 0) + 1;
    forceRHS(forceRHS.GetNumRows() - 1, 1) = forceRHS(forceRHS.GetNumRows() - 2, 1);

//        forceRHS.Info();

    rTimeIntegrationScheme.SetTimeDependentLoadCase(0, forceRHS);

    rTimeIntegrationScheme.SetMaxTimeStep(10);
    rTimeIntegrationScheme.SetMinTimeStep(0.0001 * myIntegrationScheme.GetMaxTimeStep());

    //set output during the simulation to false
    myStructure.SetShowTime(false);

    rTimeIntegrationScheme.AddResultTime("Time");
    rTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right", grpNodes_Right);

    //set result directory
    bool deleteResultDirectoryFirst(true);
    rTimeIntegrationScheme.SetResultDirectory(resultDir, deleteResultDirectoryFirst);

    //solve (perform Newton raphson iteration
    rTimeIntegrationScheme.Solve(simulationTime);

}

int main()
{
    try
    {
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKutta2 rk(&mStructure);
            Run(mStructure, rk);
        }
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKutta3 rk(&mStructure);
            Run(mStructure, rk);
        }
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKutta38 rk(&mStructure);
            Run(mStructure, rk);
        }
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKutta4 rk(&mStructure);
            Run(mStructure, rk);
        }
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKuttaCashKarp rk(&mStructure);
            Run(mStructure, rk);
        }
        {
            NuTo::Structure mStructure(2);
            NuTo::RungeKuttaDormandPrince rk(&mStructure);
            Run(mStructure, rk);
        }

    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "Error executing RungeKutta " << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "Error executing RungeKutta " << std::endl;
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
