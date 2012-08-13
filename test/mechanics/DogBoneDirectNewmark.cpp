#include <stdlib.h>
#include <sstream>
#include "boost/filesystem.hpp"

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

void ReplaceStringInFile(std::string rFileName, std::string rStringOld, std::string rStringNew, int rNumber)
{
	std::cout << std::string("changing parameter in file ")+rFileName+std::string(" from ")+ rStringOld+ std::string(" to ") + rStringNew << std::endl;
	std::ifstream iFile(rFileName,std::ios::binary);
	if (iFile.is_open())
	{
		iFile.seekg(0,std::ios_base::end);
		long s=iFile.tellg();
		char *buffer=new char[s];
		iFile.seekg(0);
		iFile.read(buffer,s);
		iFile.close();
		std::string fileStr(buffer,s);
		delete[] buffer;
		size_t off=0;
		int numReplaced(0);
		while ((off=fileStr.find(rStringOld,off))!=std::string::npos)
		{
			fileStr.replace(off,rStringOld.length(),rStringNew);
			off+=rStringNew.length();
			numReplaced++;
		}
		if(numReplaced!=rNumber)
		{
			//check, if new string is already existing
			int numFound(0);
			off = 0;
			while ((off=fileStr.find(rStringNew,off))!=std::string::npos)
			{
				off+=rStringNew.length();
				numFound++;
			}
			if (numFound!=rNumber)
			{
				std::cout << " to be replaced " << numReplaced << " found total after replacement " << numFound << "!=" << rNumber << " (given as input)" << "\n";
				throw NuTo::Exception(std::string("[ReplaceStringInFile] Error changing parameter in file ")+rFileName+std::string(" from ")+
					rStringOld+ std::string(" to ") + rStringNew);
			}
		}
		std::ofstream ofile(rFileName);
		ofile.write(fileStr.c_str(),fileStr.size());
	}
	else
	{
		std::cout << "File " << rFileName << " could not be opened." << "\n";
		exit(0);
	}
}


int main(int argc,char *argv[])
{
try
{
    std::string resultDir;
    std::string execDir;
	if (argc==1)
    {
//#ifdef DEBUG
//    this is only relevant if you run the program via eclipse
//    resultDir = std::string("/home/unger3/develop/nutoBuildMyExampleDebug/test/mechanics/DogBoneDirectNewmarkresults");
//    #else
//#endif
        resultDir = boost::filesystem::initial_path().string()+std::string("/resultsDogBoneDirectNewmark");
    }
    else
    {
    	if (argc!=2)
    	{
    	    std::cout << "Error executing DogboneDirectNewmark, expecting one arguments, got "<< argc-1 << std::endl;
    		exit(-1);
    	}
    	resultDir = argv[1];
    }

	int dimension(2);

    NuTo::Structure myStructure(dimension);

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

    //copy the geo file into the result directory
    boost::filesystem::path resultFile(resultDir);
    resultFile/= std::string("DogBoneDirectNewmarkGeometry.geo");

    boost::filesystem::path srcFile(resultFile.parent_path().parent_path());
    srcFile/= std::string("DogBoneDirectNewmarkGeometry.geo");

    //std::cout << "copy " << srcFile.string() << " to " << resultFile.string() << "\n";
    boost::filesystem::copy_file (srcFile, resultFile);

    int mSeed=2;

    double D(0.05); //dog bone specimen
    std::stringstream streamD;
    streamD << D;
    ReplaceStringInFile(resultFile.string(), std::string("D = 0.5;"), std::string("D = ")+streamD.str()+std::string(";"),1);
    double meshSize(D*0.2); //dog bone specimen
    std::stringstream streamM;
    streamM << meshSize;
    ReplaceStringInFile(resultFile.string(), std::string("lc = 0.005;"), std::string("lc = ")+streamM.str()+std::string(";"),1);

    //mesh
    std::string gmshCommand(std::string("cd ")+ resultDir+std::string(";gmsh -2 -order 1 DogBoneDirectNewmarkGeometry.geo"));
    std::cout << gmshCommand << "\n";
    std::cout<< "in the test file, the gmsh command is not executed. If you want to use it, uncomment the next line and make sure, ./gmsh is added to the path." << "\n";
    boost::filesystem::path mshFile(resultDir);
    mshFile /= std::string("DogBoneDirectNewmarkGeometry.msh");
    if (system(gmshCommand.c_str())!=0)
    {
        boost::filesystem::path srcFileMsh(resultFile.parent_path().parent_path());
        srcFileMsh/= std::string("DogBoneDirectNewmarkGeometry.msh");
        boost::filesystem::copy_file (srcFileMsh, mshFile);
    }

    //import mesh
    NuTo::FullMatrix<int> createdGroupIds;
    myStructure.ImportFromGmsh(mshFile.string(),2,"displacements", "CONSTITUTIVELAWIPNONLOCAL", "STATICDATANONLOCAL", createdGroupIds);

	//section
	double thickness(1);
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection,thickness);
    myStructure.ElementTotalSetSection(mySection);

	//constitutive
    double mDensityConcrete = 2500;       //kg and m^3
    bool isElastic(false);
    int myMat(0);
    if (isElastic)
    {
		myMat = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    }
    else
    {
		myMat = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticityEngineeringStress");
	    double fct(3e6);
		myStructure.ConstitutiveLawSetTensileStrength(myMat,fct);
	    myStructure.ConstitutiveLawSetCompressiveStrength(myMat,10*fct);
		myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMat,fct*12.5);
	    myStructure.ConstitutiveLawSetFractureEnergy(myMat,50);
		double nonlocalRadius(meshSize*1.2);
		myStructure.ConstitutiveLawSetNonlocalRadius(myMat,nonlocalRadius);
    }

	myStructure.ConstitutiveLawSetYoungsModulus(myMat,30e9);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMat,0.2);
	myStructure.ConstitutiveLawSetDensity(myMat,mDensityConcrete);

    myStructure.ElementTotalSetConstitutiveLaw(myMat);

    //Build nonlocal elements
    if (!isElastic)
	    myStructure.BuildNonlocalData(myMat);

    //wavespeed
    double nu(myStructure.ConstitutiveLawGetPoissonsRatio(myMat));
    double E(myStructure.ConstitutiveLawGetYoungsModulus(myMat));
    double rho(myStructure.ConstitutiveLawGetDensity(myMat));

    double cs(sqrt(E/(2.*(1.+nu))/rho));
    double kappaPlaneStrain(3.-4.*nu); //plane strain
    double kappaPlaneStress((3.-nu)/(1.+nu)); //plane stress
    double cp(cs*sqrt((kappaPlaneStress+1.)/(kappaPlaneStress-1.)));

    std::cout << "S-wave speed " << cs << "\n";
    //std::cout << "P-wave speed plane strain" << cs*sqrt((kappaPlaneStrain+1.)/(kappaPlaneStrain-1.)) << "\n";
    std::cout << "P-wave speed plane stress " << cp << "\n";
    std::cout << "P-wave speed solid " << sqrt(E*(1.-nu)/(rho*(1.+nu)*(1.-2.*nu))) << "\n";
    std::cout << "time to arrive at other end " << 1.5*D/cp << "\n";


    //create node groups bottom boundary
	int grpNodes_Bottom = myStructure.GroupCreate("Nodes");
	int direction=1;
	double min=0.;
	double max=0.;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Bottom,direction,min,max);

	//top boundary
	int grpNodes_Top = myStructure.GroupCreate("Nodes");
	direction=1;
	min=1.5*D;
	max=1.5*D;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Top,direction,min,max);

	//left boundary
	int grpNodes_Left = myStructure.GroupCreate("Nodes");
	direction=0;
	min=0.;
	max=0.;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Left,direction,min,max);

	//right boundary
	int grpNodes_Right = myStructure.GroupCreate("Nodes");
	direction=0;
	min=D;
	max=D;
	myStructure.GroupAddNodeCoordinateRange(grpNodes_Right,direction,min,max);

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
	myStructure.AddVisualizationComponentVelocity();
	myStructure.AddVisualizationComponentAcceleration();
	myStructure.AddVisualizationComponentEngineeringStrain();
	myStructure.AddVisualizationComponentEngineeringStress();
	myStructure.AddVisualizationComponentDamage();
#endif

	//directionX
    NuTo::FullMatrix<double> DirectionX(2,1);
    DirectionX.SetValue(0,0,1.0);
    NuTo::FullMatrix<double> DirectionY(2,1);
    DirectionY.SetValue(1,0,1.0);

    int constraintBottomDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Bottom,DirectionY,0);
    //int constraintBottomRot = myStructure.ConstraintLinearSetRotationNodeGroup(grpNodes_Bottom,0);
    int constraintBottomLeftX = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupLeftBottomSupport,DirectionX,0);
    //mConstraintBottomLeftX = myStructure.ConstraintLinearSetDisplacementNodeGroup(mGrpNodes_Bottom,DirectionX,0);
    int constraintTopDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_Top,DirectionY,0);
    //int constraintTopRot = myStructure.ConstraintLinearSetRotationNodeGroup(grpNodes_Top,0);
    //mConstraintTopLeftX = myStructure.ConstraintLinearSetDisplacementNodeGroup(mGroupLeftTopSupport,DirectionX,0);
    //mConstraintTopLeftX = myStructure.ConstraintLinearSetDisplacementNodeGroup(mGrpNodes_Top,DirectionX,0);

    //export the initial plot
    myStructure.ElementTotalUpdateTmpStaticData();
	resultFile = resultDir;
	resultFile /= std::string("ElementOutput2D_0")+std::string(".vtk");
	myStructure.ExportVtkDataFileElements(resultFile.string());
	resultFile = resultDir;
	resultFile /= std::string("NodeOutput2D_0")+std::string(".vtk");
	myStructure.ExportVtkDataFileNodes(resultFile.string());

	myStructure.Info();

	NuTo::NewmarkDirect myIntegrationScheme;

	myIntegrationScheme.SetDampingCoefficientMass(0.0);

	double simulationTime(0.002);
	double finalDisplacement(0.00004*D);
	double relAccelerationTime(0.1); //10% of the simulation time are used to ramp the velocitity linearly from 0 to v_end
	double acceleration(finalDisplacement/((relAccelerationTime-relAccelerationTime*relAccelerationTime)*simulationTime*simulationTime));
	double timeChange(relAccelerationTime*simulationTime); //transition between constant acceleration and constand velocity

	int numLoadSteps(2);

	bool linearAcceleration(true);
	NuTo::FullMatrix<double> dispRHS;
	if (linearAcceleration)
	{
		dispRHS.Resize(10*numLoadSteps+2,2);
		for (int count=0; count<10*numLoadSteps+1; count++)
		{
			dispRHS(count,0) = (double)count/(10.*numLoadSteps)*timeChange;
			dispRHS(count,1) = 0.5 * acceleration * dispRHS(count,0) * dispRHS(count,0);
		}
		dispRHS(10*numLoadSteps+1,0) = simulationTime;
		dispRHS(10*numLoadSteps+1,1) = finalDisplacement;
	}
	else
	{
		dispRHS.Resize(2,2);
		dispRHS(0,0) = 0;
		dispRHS(0,1) = 0; //disp

		dispRHS(1,0) = simulationTime;
		dispRHS(1,1) = finalDisplacement;
	}

    //NuTo::FullMatrix<double> velocity(2,1);
	//velocity(0,0) = 0.0;
	//velocity(1,0) = finalDisplacement/simulationTime;
	//std::cout << "velocity " << velocity.Trans() << "\n";
	//myStructure.NodeGroupSetVelocities(grpNodes_Top,velocity);

    myIntegrationScheme.SetDisplacements(constraintTopDisp, dispRHS);
    myIntegrationScheme.SetMaxTimeStep(simulationTime/numLoadSteps);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetMinTimeStep(1e-5*myIntegrationScheme.GetMaxTimeStep());
    myIntegrationScheme.SetToleranceForce(1e-5);
    myIntegrationScheme.SetMaxNumIterations(10);

    //set output during the simulation to false
    myStructure.SetShowTime(false);
    //myStructure.SetNumProcessors(8);

    //set output to be calculated at the left and right nodes
    NuTo::FullMatrix<int> groupNodesReactionForces(2,1);
    groupNodesReactionForces(0,0) = grpNodes_Top;
    groupNodesReactionForces(1,0) = grpNodes_Bottom;
    myIntegrationScheme.SetGroupNodesReactionForces(groupNodesReactionForces);

    //set result directory
    bool deleteDirectory(false);
    myIntegrationScheme.SetResultDirectory(resultDir,deleteDirectory);

    //solve (perform Newton raphson iteration
    myIntegrationScheme.Solve(myStructure, simulationTime);
    //read in the result file
	resultFile = resultDir;
	resultFile /= std::string("resultAllLoadSteps.dat");

	NuTo::FullMatrix<double> result;
    result.ReadFromFile(resultFile.string());
    //result.Info(15,12,true);

    boost::filesystem::path resultFileRef(resultFile.parent_path().parent_path());
    resultFileRef /= std::string("DogboneDirectNewmarkresultAllLoadStepsRef.dat");
	NuTo::FullMatrix<double> resultRef;
	resultRef.ReadFromFile(resultFileRef.string());
    std::cout << "result \n" << result << "\n";
    std::cout << "resultFileRef \n" << resultRef << "\n";

    if ((resultRef-result).Abs().Max()>1e-4)
    {
    	std::cout << "difference " << (resultRef-result).Abs().Max() << "\n";
    	result.Info();
        std::cout<< "reference results" << std::endl;
        resultRef.Info();
        std::cout << "[DogboneDirectNewmark] result is not correct." << std::endl;
        return -1;
    }

}
catch (NuTo::MechanicsException& e)
{
    std::cout << "Error executing DogboneDirectNewmark "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
catch (NuTo::Exception& e)
{
    std::cout << "Error executing DogboneDirectNewmark "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
}


