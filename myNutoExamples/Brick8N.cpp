// $Id: Brick8N.cpp 649 2013-11-08 11:10:59Z unger3 $

#include <iostream>
#include "nuto/base/Exception.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/timeIntegration/JumpDirect.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include <boost-1_55/boost/filesystem.hpp>

#include "GetMaterialData.h"

#define GetValue(variable) GetParameterValue(#variable, variable)
#define GetName(variable) GetVariableName(#variable)

int main(int argc, char* argv[])
{
try
{
    int readFlag = false;
    if(readFlag)
    {
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> a;
    a.ReadFromFile("stiffnessMatrix.txt");
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVector2(a);
    NuTo::FullVector<double,Eigen::Dynamic> rhsVector;
    rhsVector.ReadFromFile("rhsVector.txt");
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
    //stiffnessMatrix.SetOneBasedIndexing();
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(stiffnessMatrixVector2);
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
    displacementVector.WriteToFile("disp.txt"," ");
    a = rhsVector - stiffnessMatrix * displacementVector;
    std::cout << "residual: " << a.Norm() << std::endl;
    }
    else
    {
        // definitions
//        double YoungsModulus = 20000.;
//        double PoissonsRatio = 0.2;
        double Width = 1000.;
        double Height = 1000.;
        double Length = 1000.;
        int NumElementsX = 3; // by my tests bulo 3 Jörg 6
        int NumElementsY = 3; // Jörg 6
        int NumElementsZ = 3; // Jörg 6
        double Force = 1.;
        bool EnableDisplacementControl = true;
        double BoundaryDisplacement = -2.; //by my tests bulo -3, by Jörg bulo -5

        // create structure
        NuTo::Structure myStructure(3);

        // create material law
    	int Material1 = myStructure.ConstitutiveLawCreate("DamageViscoPlasticityHardeningEngineeringStress");

    	double Density, youngsModulus, nu, TensileStrength, CompressiveStrength,
    		BiaxialCompressiveStrength, Viscosity, ViscosityExponent, DamageDistribution,
    		ViscoplasticYieldSurfaceOffset, FractureEnergy, HardeningValue, HardeningExponent;

//    	GetValue(Density);							// kg/(mm)³
//    	GetValue(youngsModulus);					// MPa
//    	GetValue(nu);								// 1
//    	GetValue(TensileStrength);					// MPa
//    	GetValue(CompressiveStrength);				// MPa
//    	GetValue(BiaxialCompressiveStrength);		// MPa
//    	GetValue(Viscosity);						// 1s
//    	GetValue(ViscosityExponent);				// 1
//    	GetValue(DamageDistribution);				// 1
//    	GetValue(ViscoplasticYieldSurfaceOffset);	// MPa
//    	GetValue(FractureEnergy);					// J/m³ = Pa
//    	GetValue(HardeningValue);					// MPa
//    	GetValue(HardeningExponent);				// 1
    	Density = 7000.e-8; // bulo 2300.e-9
    	youngsModulus = 29000.;
    	nu = 0.2;
    	TensileStrength = 6.74;
    	CompressiveStrength = 65.4;
    	BiaxialCompressiveStrength = 74;
    	Viscosity = 1.8;
    	ViscosityExponent = 2.6;
    	DamageDistribution = 0.5;
    	ViscoplasticYieldSurfaceOffset = -21.2;
    	FractureEnergy = 0.024;
    	HardeningValue = 5.8;
    	HardeningExponent = 1800;


    	myStructure.ConstitutiveLawSetDensity(Material1, Density);
    	myStructure.ConstitutiveLawSetYoungsModulus(Material1, youngsModulus);
    	myStructure.ConstitutiveLawSetPoissonsRatio(Material1, nu);
    	myStructure.ConstitutiveLawSetTensileStrength(Material1, TensileStrength);
    	myStructure.ConstitutiveLawSetCompressiveStrength(Material1, CompressiveStrength);
    	myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(Material1, BiaxialCompressiveStrength);
    	myStructure.ConstitutiveLawSetViscosity(Material1, Viscosity);
    	myStructure.ConstitutiveLawSetViscosityExponent(Material1, ViscosityExponent);
    	myStructure.ConstitutiveLawSetDamageDistribution(Material1, DamageDistribution);
    	myStructure.ConstitutiveLawSetViscoplasticYieldSurfaceOffset(Material1, ViscoplasticYieldSurfaceOffset);
    	myStructure.ConstitutiveLawSetFractureEnergy(Material1, FractureEnergy);
    	myStructure.ConstitutiveLawSetHardeningValue(Material1, HardeningValue);
    	myStructure.ConstitutiveLawSetHardeningExponent(Material1, HardeningExponent);
//		Set fatigue flag
    	myStructure.ConstitutiveLawSetFatigueExtrapolation(Material1,true);

//        int Material1 = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
//        myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus);
//        myStructure.ConstitutiveLawSetPoissonsRatio(Material1, PoissonsRatio);

        int mySection1 = myStructure.SectionCreate("VOLUME");

        // create nodes
        NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(3);
        int node = 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            nodeCoordinates(2) = (double)zCount * Height/(double)NumElementsZ;
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                nodeCoordinates(1) = (double)yCount * Width/(double)NumElementsY;
                for(int xCount = 0; xCount < NumElementsX + 1; xCount++)
                {
                    nodeCoordinates(0) = (double)xCount * Length/(double)NumElementsX;
                    //std::cout << "node: " << node << " coordinates: " << nodeCoordinates.GetValue(0,0) << "," << nodeCoordinates.GetValue(1,0) << "," << nodeCoordinates.GetValue(2,0) << std::endl;
                    myStructure.NodeCreate(node, "displacements", nodeCoordinates,2);
                    node ++;
                }
            }
        }

        // create elements
        NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(8);
        int element = 0;
        for(int zCount = 0; zCount < NumElementsZ; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY; yCount++)
            {
                for(int xCount = 0; xCount < NumElementsX; xCount++)
                {
                    int node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount;
                    elementIncidence(0) = node1;
                    elementIncidence(1) = node1 + 1;
                    elementIncidence(2) = node1 + NumElementsX + 2;
                    elementIncidence(3) = node1 + NumElementsX + 1;
                    elementIncidence(4) = node1 + (NumElementsX + 1) * (NumElementsY + 1);
                    elementIncidence(5) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1;
                    elementIncidence(6) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2;
                    elementIncidence(7) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1;
                    //std::cout << "element: " << element << " incidence: " << std::endl;
                    //elementIncidence.Info();
                    myStructure.ElementCreate(element, "Brick8N", elementIncidence,"CONSTITUTIVELAWIP","STATICDATA");
                    myStructure.ElementSetConstitutiveLaw(element,Material1);
            	    myStructure.ElementSetSection(element,mySection1);
                    element ++;
                }
            }
        }

        // boundary conditions
        NuTo::FullVector<double,Eigen::Dynamic> direction(3);
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1);
                //std::cout << "node: " << node << std::endl;
                myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0);
            }
        }
        direction(0)= 0;
        direction(1)= 0;
        direction(2)= 1;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0);
        direction(0)= 0;
        direction(1)= 1;
        direction(2)= 0;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);

        // create the set of the nodes experiencing the time-dependent active constraints
        int ActiveNodesBC = myStructure.GroupCreate("Nodes");	// myChange

        // apply nodes
        if(EnableDisplacementControl)
        {
            std::cout << "Displacement control" << std::endl;
            // boundary displacments
            direction(0)= 1;
            direction(1)= 0;
            direction(2)= 0;
            for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
            {
                //std::cout << zCount << std::endl;
                for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
                {
                    int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    //std::cout << "node: " << node << std::endl;
// Joerg            myStructure.ConstraintLinearSetDisplacementNode(node, direction, BoundaryDisplacement);
                    myStructure.GroupAddNode(ActiveNodesBC,node);  // myChange
                }
            }
        }
        else
        {
            std::cout << "Load control" << std::endl;
            // apply load to nodes
            direction(0)= 1;
            direction(1)= 0;
            direction(2)= 0;
            for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
            {
                double nodeForce;
                if(zCount == 0 || zCount == NumElementsZ)
                {
                    nodeForce = Force / (4 *NumElementsY * NumElementsZ);
                }
                else
                {
                    nodeForce = Force / (2 *NumElementsY * NumElementsZ);
                }
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX;
                //std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
                myStructure.LoadCreateNodeForce(1, node, direction, nodeForce);
                for(int yCount = 1; yCount < NumElementsY; yCount++)
                {
                    node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    //std::cout << "apply force to node: " << node << " force: " << 2 * nodeForce << std::endl;
                    myStructure.LoadCreateNodeForce(1, node, direction, 2 * nodeForce);
                }
                node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1;
                //std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
                myStructure.LoadCreateNodeForce(1 ,node, direction, nodeForce);
            }
        }

        // pseudo time-dependent analysis
        std::cout << "#ActiveNodeBC = " << myStructure.GroupGetNumMembers(ActiveNodesBC) << std::endl;
        // start analysis
//        NuTo::NewmarkDirect myIntegrationScheme(&myStructure);	//	NewmarkDirect
        NuTo::JumpDirect myIntegrationScheme(&myStructure);			//	JumpDirect
        myIntegrationScheme.SetDynamic(true);
        myIntegrationScheme.SetDampingCoefficientMass(0.0);

        // direction of the time-dependent constraint
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;

        // displacement time dependent factor
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> displacementsBC(2,2);
        displacementsBC(0,0) = 0.;   // t=0
        displacementsBC(0,1) = 0.;   // u(0)=0
        displacementsBC(1,0) = 1.;   // t=1
        displacementsBC(1,1) = 1.;   // u(1)= 1

        displacementsBC.col(1) *= BoundaryDisplacement;

        // apply time-dependent constraints to the node set ActiveNodeBC
        int constraintLoading = myStructure.ConstraintLinearSetDisplacementNodeGroup(ActiveNodesBC,direction,0.);
        myIntegrationScheme.SetTimeDependentConstraint(constraintLoading, displacementsBC);
   		myIntegrationScheme.SetHarmonicIncrementation(30);


        // apply harmonic excitation
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> HarmonicConstraintFactor(1,3);	//	JumpDirect
        HarmonicConstraintFactor(0,0) = 0.1*BoundaryDisplacement;	//	set displacement
        HarmonicConstraintFactor(0,1) = 2.;							//	set frequency
        HarmonicConstraintFactor(0,2) = 2.9;						// 	set number of cycles

        myIntegrationScheme.SetHarmonicConstraint(HarmonicConstraintFactor);		//	JumpDirect
   		myIntegrationScheme.SetHarmonicExtrapolation(true);

        // specify solver parameters
        myIntegrationScheme.SetToleranceForce(1.e-4);
        myStructure.SetShowTime(true);
        myIntegrationScheme.SetAutomaticTimeStepping(true);
        myIntegrationScheme.SetTimeStep(0.05);

        // output and visualization for the results
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.AddVisualizationComponentDamage();
        myStructure.ExportVtkDataFileElements("Brick8N.vtk");

        boost::filesystem::path resultDir = boost::filesystem::current_path();
        std::cout << resultDir << std::endl;
//        resultDir += std::string("/Brick8NOutput");	//	NewmarkDirect
        resultDir += std::string("/Brick8NOutputJump");		//	JumpDirect
        std::cout << resultDir << std::endl;
        boost::filesystem::create_directory(resultDir);
//        boost::filesystem::path resultFile = resultDir;
//        std::cout << resultFile << std::endl;
//        resultFile /= "out1.vtk";

        myIntegrationScheme.AddResultGroupNodeForce("Out_NodeForceBC",ActiveNodesBC);
        myIntegrationScheme.AddResultNodeDisplacements("Out_NodeDisplBC",1);
        myIntegrationScheme.AddResultTime("Out_Time");
        myIntegrationScheme.SetResultDirectory(resultDir.string(),true);

        // solving

        myStructure.CalculateMaximumIndependentSets();
        try{
            myIntegrationScheme.Solve(100); // solve untill time = 1
        }
        	catch(NuTo::Exception& e)
    	{
    		std::cout << e.ErrorMessage() << std::endl;
    	}

        myStructure.Info();


        // Joerg's analysis
//        // start analysis
//        // build global dof numbering
//        myStructure.NodeBuildGlobalDofs();
//
//        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
//        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVector2;
//        NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
//        myStructure.CalculateMaximumIndependentSets();
//        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixVector2, dispForceVector);
//        stiffnessMatrixVector2.RemoveZeroEntries(0,1e-14);
//        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> A(stiffnessMatrix);
//        //A.WriteToFile("stiffnessMatrix.txt"," ");
//        //stiffnessMatrix.Info();
//        //dispForceVector.Info();
//
//        // build global external load vector
//        NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
//        myStructure.BuildGlobalExternalLoadVector(1,extForceVector);
//        //extForceVector.Info();
//
//        // calculate right hand side
//        NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;
//        rhsVector.WriteToFile("rhsVector.txt"," ");
//
//        // solve
//        NuTo::SparseDirectSolverMUMPS mySolver;
//        NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(stiffnessMatrixVector2);
//        stiffnessMatrix.SetOneBasedIndexing();
//        mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
//        displacementVector.WriteToFile("displacementVector.txt"," ");
//
//        // write displacements to node
//        myStructure.NodeMergeActiveDofValues(displacementVector);
//
//        // calculate residual
//        displacementsBC(0,0) = 0.;   // t=0
//        displacementsBC(0,1) = 0.;   // u(0)=0
//        displacementsBC(1,0) = 1.;   // t=2
//        displacementsBC(1,1) = 1.;   // u(2)= 1
//        NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
//        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
//        NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
//        std::cout << "residual: " << residualVector.Norm() << std::endl;
//
//
//
//#ifdef ENABLE_VISUALIZE
//        // visualize results
//        myStructure.AddVisualizationComponentDisplacements();
//        myStructure.AddVisualizationComponentEngineeringStrain();
//        myStructure.AddVisualizationComponentEngineeringStress();
//        myStructure.ExportVtkDataFileElements("Brick8N.vtk");
//#endif // ENABLE_VISUALIZE
    }
}
catch (NuTo::Exception& e)
{
	std::cout << "Error executing Brick8N "<< std::endl;
	std::cout << e.ErrorMessage() << std::endl;
	return(-1);
}
return 0;
}

