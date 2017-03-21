#include <sstream>
#include "boost/filesystem.hpp"

#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/RungeKutta3.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/RungeKuttaCashKarp.h"
#include "mechanics/timeIntegration/RungeKuttaDormandPrince.h"

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/sections/SectionPlane.h"
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

    //section
    double thickness(1);
    auto mySection = NuTo::SectionPlane::Create(thickness, false);
    myStructure.ElementTotalSetSection(mySection);

    //constitutive
    int myMatLattice = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice, NuTo::Constitutive::eConstitutiveParameter::DENSITY, density);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLattice);



    auto meshInfo = NuTo::MeshGenerator::Grid(myStructure, {mL, mH},
                                              {numElementsX,numElementsY}, NuTo::Interpolation::eShapeType::QUAD2D);

    myStructure.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    myStructure.ElementGroupSetSection(meshInfo.first, mySection);
    myStructure.ElementGroupSetConstitutiveLaw(meshInfo.first, myMatLattice);


    myStructure.ElementTotalConvertToInterpolationType();

    //create node groups bottom boundary
    int nOrigin = myStructure.NodeGetIdAtCoordinate(Eigen::Vector2d::Zero(), 1.e-6);

    //right boundary
    int grpNodes_Right = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Right, 0, mL, mL);

    int grpNodes_Left = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNodeCoordinateRange(grpNodes_Right, 0, 0., 0.);

    //calculate maximum independent sets
    myStructure.CalculateMaximumIndependentSets();

    //set constraints
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, Eigen::Vector2d::UnitY(), 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Left, Eigen::Vector2d::UnitX(), 0);

    myStructure.LoadCreateNodeGroupForce(0, grpNodes_Right, Eigen::Vector2d::UnitX(), 1000);

    myStructure.Info();
    NuTo::RungeKutta4 myIntegrationScheme(&myStructure);

    double simulationTime(1);

    //set a sinusoidal quarter wave
    double period(5);
    Eigen::MatrixXd forceRHS(101, 2);
    for (int count = 0; count < forceRHS.rows() - 1; count++)
    {
        double t = ((double) count) / ((double) forceRHS.rows() - 2.) * 0.25 * period;
        forceRHS(count, 0) = t;
        forceRHS(count, 1) = 0.1 * sin(t / period * 2. * M_PI);
        //loadRHSFactor(count,0) = 0;
    }
    forceRHS(forceRHS.rows() - 1, 0) = forceRHS(forceRHS.rows() - 2, 0) + 1;
    forceRHS(forceRHS.rows() - 1, 1) = forceRHS(forceRHS.rows() - 2, 1);

//        forceRHS.Info();

    rTimeIntegrationScheme.SetTimeDependentLoadCase(0, forceRHS);

    rTimeIntegrationScheme.SetMaxTimeStep(10);
    rTimeIntegrationScheme.SetMinTimeStep(0.0001 * myIntegrationScheme.GetMaxTimeStep());

    //set output during the simulation to false
    myStructure.SetShowTime(false);

    rTimeIntegrationScheme.AddResultTime("Time");
    rTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right", grpNodes_Right);
    rTimeIntegrationScheme.AddResultElementIpData("StressElement",1,NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS); //only for test issues

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
