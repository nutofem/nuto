#include <stdlib.h>
#include <sstream>
#include "boost/filesystem.hpp"

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#define createResult true

int main(int argc,char *argv[])
{
try
{
    std::string resultDir;
	if (argc==1)
    {
#ifdef DEBUG
        resultDir = std::string("./ResultsNewmarkPlane2D4N");
#else
        resultDir = boost::filesystem::initial_path().string()+std::string("/ResultsNewmarkPlane2D4N");
#endif
    }
    else
    {
    	if (argc!=2)
    	{
    	    std::cout << "Error executing NewmarkPlane2D4N, expecting one argument, got "<< argc-1 << std::endl;
    		exit(-1);
    	}
    	resultDir = argv[1];
    }

    //delete result directory
    if (boost::filesystem::exists(resultDir))    // does p actually exist?
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
    double elementLength = 1;
    double elementHeight = 10;
    double density = 1;       //N and mm^3

    //create structure
    NuTo::Structure myStructure(2);

	//create nodes
    int numNodesX = (int)(floor(mL/elementLength))+1;
    int numNodesY = (int)(floor(mH/elementHeight))+1;
    double deltaX = mL/(numNodesX-1);
    double deltaY = mH/(numNodesY-1);

    int nodeNum(0);
    for (int countY=0; countY<numNodesY; countY++)
    {
        for (int countX=0; countX<numNodesX; countX++)
        {
        	NuTo::FullVector<double,Eigen::Dynamic> coordinates(2);
        	coordinates(0) = countX*deltaX;
        	coordinates(1) = countY*deltaY;
        	myStructure.NodeCreate(nodeNum,std::string("DISPLACEMENTS"),coordinates,2);
        	nodeNum++;
        }
    }

	//create elements
    int numElementsX = numNodesX-1;
    int numElementsY = numNodesY-1;
    for (int countY=0; countY<numElementsY; countY++)
    {
        for (int countX=0; countX<numElementsX; countX++)
        {
        	NuTo::FullVector<int,Eigen::Dynamic> nodes(4);
        	nodes(0) = countX  +  countY   *numNodesX;
        	nodes(1) = countX+1+  countY   *numNodesX;
        	nodes(2) = countX+1+ (countY+1)*numNodesX;
        	nodes(3) = countX  + (countY+1)*numNodesX;
        	myStructure.ElementCreate(std::string("PLANE2D4N"), nodes, std::string("ConstitutiveLawIp"), std::string("StaticData"));
        }
    }

    //section
	double thickness(1);
    int mySection = myStructure.SectionCreate("Plane_Stress");
    myStructure.SectionSetThickness(mySection,thickness);
    myStructure.ElementTotalSetSection(mySection);

	//constitutive
    int myMatLattice = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLattice,30000);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLattice,0.0);
    myStructure.ConstitutiveLawSetDensity(myMatLattice,density);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLattice);

    //wavespeed
    //double nu(myStructure.ConstitutiveLawGetPoissonsRatio(myMatLattice));
    //double E(myStructure.ConstitutiveLawGetYoungsModulus(myMatLattice));
    //double rho(myStructure.ConstitutiveLawGetDensity(myMatLattice));
    //std::cout << "P-wave speed " << sqrt(E*(1.-nu)/((1.+nu)*(1.-2*nu))) << "\n";
    //std::cout << "nu E " << nu << " " << E << "\n";

    //double cs(sqrt(E/(2.*(1.+nu))/rho));
    //double kappaPlaneStrain(3.-4.*nu); //plane strain
    //double kappaPlaneStress((3.-nu)/(1.+nu)); //plane stress
    //double cp(cs*sqrt((kappaPlaneStress+1.)/(kappaPlaneStress-1.)));

    //std::cout << "S-wave speed" << cs << "\n";
    //std::cout << "P-wave speed plane strain" << cs*sqrt((kappaPlaneStrain+1.)/(kappaPlaneStrain-1.)) << "\n";
    //std::cout << "P-wave speed plane stress" << cp << "\n";
    //std::cout << "time to arrive at other end" << mL/(cp) << "\n";

    //create node groups bottom boundary
	int grpNodes_Bottom = myStructure.GroupCreate("Nodes");
	int direction=1;
	double min=0-0.01*elementHeight;
	double max=0+0.01*elementHeight;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Bottom,direction,min,max);

	//top boundary
	int grpNodes_Top = myStructure.GroupCreate("Nodes");
	direction=1;
	min=mH-0.01*elementHeight;
	max=mH+0.01*elementHeight;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Top,direction,min,max);

	//left boundary
	int grpNodes_Left = myStructure.GroupCreate("Nodes");
	direction=0;
	min=0;
	max=0;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Left,direction,min,max);

	//right boundary
	int grpNodes_Right = myStructure.GroupCreate("Nodes");
	direction=0;
	min=mL-0.01*elementLength;
	max=mL+0.01*elementLength;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Right,direction,min,max);

	//allNodes
	int grpNodes_All = myStructure.GroupCreate("Nodes");
	direction=0;
	min=0-0.01*elementLength;
	max=mL+0.01*elementLength;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_All,direction,min,max);

	//intersect the groups for the left bottom node
	int groupLeftBottomSupport = myStructure.GroupIntersection(grpNodes_Bottom,grpNodes_Left);
	if (myStructure.GroupGetNumMembers(groupLeftBottomSupport)!=1)
	{
		std::cout << "group for bottom left boundary node has " << myStructure.GroupGetNumMembers(groupLeftBottomSupport) << " members."<< std::endl;
	    exit(-1);
	}

	int groupLeftTopSupport = myStructure.GroupIntersection(grpNodes_Top,grpNodes_Left);
	if (myStructure.GroupGetNumMembers(groupLeftTopSupport)!=1)
	{
		std::cout << "group for top left boundary node has " << myStructure.GroupGetNumMembers(groupLeftTopSupport) << " members."<< std::endl;
	    exit(-1);
	}

    //calculate maximum independent sets
	myStructure.CalculateMaximumIndependentSets();

#ifdef ENABLE_VISUALIZE
    //add visualization
    //myStructure.AddVisualizationComponentSection();
    //myStructure.AddVisualizationComponentConstitutive();
	myStructure.AddVisualizationComponentDisplacements();
	myStructure.AddVisualizationComponentEngineeringStrain();
	myStructure.AddVisualizationComponentEngineeringStress();
    //myStructure.AddVisualizationComponentDamage();
#endif

    //set constraints
	//directionX
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionX(2,1);
    DirectionX.SetValue(0,0,1.0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionY(2,1);
    DirectionY.SetValue(1,0,1.0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionXY(2,1);
    DirectionXY.SetValue(0,0,1.0);
    DirectionXY.SetValue(1,0,1.0);

     //export the initial plot
	boost::filesystem::path resultFile;
	resultFile = resultDir;
	resultFile /= std::string("ElementOutput2D_0")+std::string(".vtk");
	myStructure.ExportVtkDataFileElements(resultFile.string());
	resultFile = resultDir;
	resultFile /= std::string("ParticleOutput2D_0")+std::string(".vtk");
	myStructure.ExportVtkDataFileNodes(resultFile.string());

	myStructure.Info();

	//NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> velocity(2,1);
	//velocity(0,0) = 0.1;
	//velocity(1,0) = 0.0;
	//myStructure.NodeGroupSetVelocities(grpNodes_All,velocity);

	//create a unit load on all the nodes on the left side
//	myStructure.LoadCreateNodeGroupForce(grpNodes_Left,DirectionX , 1);
	myStructure.LoadCreateNodeGroupForce(0,groupLeftBottomSupport,DirectionX, 1);

	NuTo::NewmarkDirect myIntegrationScheme;

	myIntegrationScheme.SetDampingCoefficientMass(0.05);
	myIntegrationScheme.SetDynamic(true);

	//set a sinusoidal load in xy-direction(diagonal) on the lower left node
	double simulationTime(10.);
    //step wise loading
/*    double maxAmplitude(0.0);
    double period(20);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> loadRHSFactor(simulationTime/period*10+1,2);
    for (int count=0; count<loadRHSFactor.GetNumRows(); count++)
    {
    	double t = count*simulationTime/(loadRHSFactor.GetNumRows()-1);
    	loadRHSFactor(count,0) = t;
    	loadRHSFactor(count,0) = maxAmplitude*sin(t/period*2.*M_PI);
    	//loadRHSFactor(count,0) = 0;
    }
    //myIntegrationScheme.SetExternalLoads(time, loadRHSFactor);
    //std::cout << "loadRHSFactor " << loadRHSFactor << "\n";
*/
	//set a sinusoidal quarter wave
    double period(5);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispRHS(51,2);
    for (int count=0; count<dispRHS.GetNumRows()-1; count++)
    {
    	double t = ((double)count)/((double)dispRHS.GetNumRows()-2.)*0.25*period ;
    	dispRHS(count,0) = t;
    	dispRHS(count,1) = 0.1*sin(t/period*2.*M_PI);
    	//loadRHSFactor(count,0) = 0;
    }
	dispRHS(dispRHS.GetNumRows()-1,0) = dispRHS(dispRHS.GetNumRows()-2,0)+1;
	dispRHS(dispRHS.GetNumRows()-1,1) = dispRHS(dispRHS.GetNumRows()-2,1);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Left,DirectionX,0);
    int constraintRightDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Right,DirectionX,0);
    myIntegrationScheme.SetTimeDependentConstraint(constraintRightDisp, dispRHS);

    //5 timesteps to capture the quarter wave (5/quarter wave)
    //myIntegrationScheme.SetMaxTimeStep(period/20.);
    myIntegrationScheme.SetMaxTimeStep(10);
    myIntegrationScheme.SetMinTimeStep(0.001*myIntegrationScheme.GetMaxTimeStep());

    //set output during the simulation to false
    myStructure.SetShowTime(false);
    myStructure.SetNumProcessors(8);

    //set output to be calculated at the left and right nodes
    NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> mGroupNodesReactionForces(2,1);
    mGroupNodesReactionForces(0,0) = grpNodes_Left;
    mGroupNodesReactionForces(1,0) = grpNodes_Right;
    myIntegrationScheme.SetGroupNodesReactionForces(mGroupNodesReactionForces);

    //set result directory
    bool deleteResultDirectoryFirst(true);
    myIntegrationScheme.SetResultDirectory(resultDir,deleteResultDirectoryFirst);

    //solve (perform Newton raphson iteration
    myIntegrationScheme.Solve(myStructure, simulationTime);

    //read in the result file
	resultFile = resultDir;
	resultFile /= std::string("resultAllLoadSteps.dat");

	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> result;
    result.ReadFromFile(resultFile.string());
    std::cout << "result" << std::endl;
    result.Info(15,12,true);

	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> resultRef(2,11);
	resultRef(0,7) = -1; //load on fixed node
	resultRef(1,0) = 10;
	resultRef(1,1) = 0.1;
	resultRef(1,2) = 1.500009228786e+01;
	resultRef(1,3) = 6.452563911884e-02;
	resultRef(1,4) = 1.700347727382e-02;
	resultRef(1,5) = 1.508162134529e+01;
	resultRef(1,6) = 0.;
	resultRef(1,7) = -3.001683921119e+02;
	resultRef(1,8) = 0.;
	resultRef(1,9) = 3.016324269061e+02;
	resultRef(1,10) = 0.;

    if ((resultRef-result).cwiseAbs().maxCoeff()>1e-4)
    {
    	std::cout << "difference " << (resultRef-result).cwiseAbs().maxCoeff() << "\n";
        std::cout<< "real result" << std::endl;
    	result.Info();
        std::cout<< "ref result" << std::endl;
        resultRef.Info();
    	std::cout << "difference " << (resultRef-result) << "\n";
        std::cout << "[NewmarkPlane2D4N] result is not correct." << std::endl;
        return -1;
    }
}
catch (NuTo::MechanicsException& e)
{
    std::cout << "Error executing NewmarlPlane2D4N "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
catch (NuTo::Exception& e)
{
    std::cout << "Error executing NewmarlPlane2D4N "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
}


