#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/nodes/NodeDisplacements3D.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif

#include <iostream>
#include "nuto/base/NuToObject.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/ConjugateGradientStructureGrid.h"

int main()
{
	//double microTomm =0.001;
		double microTomm =1;

	//   int readFlag = false;
    double PoissonsRatio = 0.2;
    //for local base stiffness matrix
    double YoungsModulus = 1.;

    int numVoxel;
    const double* voxelSpacing;
    const double* gridOrigin;
    const int* gridDimension;

    // create structure
    try
    {
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
        NuTo::StructureGrid myGrid(3);

		// read entries
		//myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("../nuto/examples/c++/InputStructureGrid3D");
		//myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("InputCheckerboard");
		myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("InputTest");
		numVoxel=myGrid.GetNumVoxels();
		voxelSpacing=myGrid.GetVoxelSpacing();
		gridOrigin=myGrid.GetGridOrigin();
		gridDimension=myGrid.GetGridDimension();
		assert (numVoxel==gridDimension[0]*gridDimension[1]*gridDimension[2]);
		std::cout<<__FILE__<<"  variab: spac: "<<voxelSpacing[0]<< " gridDim: "<<gridDimension[0]<<std::endl;
		NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> imageValues (numVoxel,1);
		//imageValues.NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>::ImportFromVtkASCIIFile( "../nuto/examples/c++/InputStructureGrid3D");
		//imageValues.NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>::ImportFromVtkASCIIFile( "InputCheckerboard");
		imageValues.NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>::ImportFromVtkASCIIFile( "InputTest");

		std::cout<<"first value "<< imageValues(0,0) << std::endl;
		std::cout<<"numVoxel "<< numVoxel << std::endl;

		//RB
		//double Force = 1.;
		bool EnableDisplacementControl = true;
		double BoundaryDisplacement = 1.0;
		//double BoundaryDisplacement = (gridDimension[0]*voxelSpacing[0]*microTomm)/20.0;
		std::cout<<__FILE__<<"  Boundary Displacement: "<<BoundaryDisplacement<<std::endl;
    	//calculate one element stiffness matrix with E=1
        NuTo::Structure myHelpStruc(3);

        myHelpStruc.SetVerboseLevel(1);
        // create material law
        int myMat=myHelpStruc.ConstitutiveLawCreate("LinearElastic");
        myHelpStruc.ConstitutiveLawSetPoissonsRatio(myMat, PoissonsRatio);
        myHelpStruc.ConstitutiveLawSetYoungsModulus(myMat, YoungsModulus);

        // create nodes
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> nodeCoordinates(3, 1);
        NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> elementIncidence(8,1);

        nodeCoordinates(0, 0) = -voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] * microTomm / 2;
        elementIncidence(0,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] * microTomm / 2;
        elementIncidence(1,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] * microTomm / 2;
        elementIncidence(2,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] * microTomm / 2;
        elementIncidence(3,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] * microTomm / 2;
        elementIncidence(4,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] * microTomm / 2;
        elementIncidence(5,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] * microTomm / 2;
        elementIncidence(6,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] * microTomm / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] * microTomm / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] * microTomm / 2;
        elementIncidence(7,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        // first element create

       // elementIncidence.Info();
        int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
        //myHelpStruc.ElementSetConstitutiveLaw(element,"Material1");
        myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);

        //myHelpStruc.NodeInfo(0);
        myHelpStruc.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrix;
        NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rows;
        NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> coluums;
        myHelpStruc.ElementStiffness(0,stiffnessMatrix,rows,coluums );

        std::cout<<"Element Stiffness created"<<std::endl;
        //stiffnessMatrix.Info(12,6);
  		stiffnessMatrix.WriteToFile("stiffnessMatrixOld.txt"," ");

		//grid structure create
		//material values are smaller than threshold value
		//color values form 0 255
		//thresholdMaterialValues: 180, 188 ->last value with material
		int thresholdMaterialValue=188; //last value for "air", voxel which does not get element
		myGrid.CreateNodeGrid("DISPLACEMENTS",thresholdMaterialValue);

		//myGrid.CreateNodeGrid("COORDINATES");
		std::cout<<"NodeGrid created"<<std::endl;
		int numNodes=myGrid.GetNumNodes();
		std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;


		//set Modul for each color
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> myMapColorModul(256,1);
		for(int count=0;count<thresholdMaterialValue;count++)
			myMapColorModul(count,0)=100000.;
		for(int count=thresholdMaterialValue;count<255;count++)
			myMapColorModul(count,0)=0.;

		//myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");
		//generiert Knoten, Freiheitsgrade noch nicht gesetzt


		myGrid.CreateElementGrid(stiffnessMatrix,myMapColorModul,"VOXEL8N");
		std::cout<<"ElementGrid created with "<<myGrid.GetNumElements()<<" Elements" <<std::endl;

		// boundary conditions
		int NumElementsX = (int) gridDimension[0];
		int NumElementsY = (int) gridDimension[1];
		int NumElementsZ = (int) gridDimension[2];

		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> direction(3,1);
		direction(0,0)= 1;
		direction(1,0)= 0;
		direction(2,0)= 0;
		std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;

		int numGridNodes=(NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);
		for (int count = 0;count<numGridNodes;count+= (NumElementsX + 1))
      	{
			//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint x"<< count <<std::endl;
			//new way of constraint saving
			myGrid.NodeSetConstraintSwitch(count,0,false);
      	}

		direction(0,0)= 0;
		direction(1,0)= 1;
		direction(2,0)= 0;
		// Test Randbedingungen y=0
		// for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);++count)
		for (int count =0;count < numGridNodes;count+= (NumElementsX + 1)*(NumElementsY + 1))
		{
			//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint y "<< count <<std::endl;
			//new way of constraint saving
			myGrid.NodeSetConstraintSwitch(count,1,false);
		}

		direction(0,0)= 0;
		direction(1,0)= 0;
		direction(2,0)= 1;
		// Test Randbedingungen z=0
		//for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);++count)
		for (int count = 0;count<(NumElementsX + 1)*(NumElementsY + 1);count+=(NumElementsX + 1))
		{
			//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint z "<< count <<std::endl;
			//new way of constraint saving
			myGrid.NodeSetConstraintSwitch(count,2,false);
		}

		// apply nodes
		if(EnableDisplacementControl)
		{
			std::cout << "Displacement control" << std::endl;
			// boundary displacments
			direction(0,0)= 1;
			direction(1,0)= 0;
			direction(2,0)= 0;
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> displacements(3,1);
			displacements(0,0)= BoundaryDisplacement;
			displacements(1,0)= 0;
			displacements(2,0)= 0;

			for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
			{
				//std::cout << zCount << std::endl;
				for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
				{
					int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
					// std::cout << "node displacement constraint: " << node << std::endl;
					if(myGrid.NodeGridGetNodePtr(node))
						//stays important for both ways
						myGrid.NodeSetDisplacements(node, displacements);

					//new way of constraint saving
					myGrid.NodeSetConstraintSwitch(node,0,false);
				}
			}
		}
		else
		{
			std::cout <<__FILE__<<" "<<__LINE__<< "Load control" << "not implemented"<<std::endl;
			return 1;
			//! @TODO: Add special configurations if edge node not exist
			// apply load to nodes
			/*
			direction(0,0)= 1;
			direction(1,0)= 0;
			direction(2,0)= 0;
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
				int myNodeNumber;
				int flag=0;
				try
				{
					myNodeNumber=myGrid.NodeGetNodeNumberFromId(node); //node from type NodeGridCoordinates
				}
				catch(NuTo::MechanicsException& e)
				{
					std::cout<<__FILE__<<__LINE__<<"node existiert nicht"<<node<<std::endl;
					flag=1;
				}
				if(flag==0)
				{
					std::cout<<__FILE__<<__LINE__<<"myNodeNumber"<<myNodeNumber<<"  von  "<<myGrid.GetNumNodes()<<std::endl;
					myGrid.LoadCreateNodeForce(myNodeNumber, direction, nodeForce);
				}
				for(int yCount = 1; yCount < NumElementsY; yCount++)
				{
					node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
					std::cout << "apply force to node: " << node << " force: " << 2 * nodeForce << std::endl;

					myGrid.LoadCreateNodeForce(node, direction, 2 * nodeForce);
				}
				node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1;
				//std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
				myGrid.LoadCreateNodeForce(node, direction, nodeForce);
			}
*/

		}
		// start analysis
		std::cout<<__FILE__<<" "<<__LINE__<<"  start analysis"<<std::endl;
		// build global dof numbering
		myGrid.SetVerboseLevel(0);
		//myGrid.NodeBuildGlobalDofs();

		std::cout<<__FILE__<<" "<<__LINE__<<"  glob dofs "<<myGrid.GetNumNodes()*3<<std::endl;
		//std::cout<<__FILE__<<" "<<__LINE__<<" active dofs "<<myGrid.GetNumActiveDofs()<<std::endl;

		//NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> *voxelLocation(myGrid.GetNumElements(),4);
		//voxelLocation=myGrid.GetVoxelNumAndLocMatrix();
		//  std::cout<<__FILE__<<" "<<__LINE__<<" VoxelLocationmatrix: \n"<<voxelLocation<< "\n\n"<<std::endl;
		/*
		int *dofs;

		for (int n=0;n< myGrid.GetNumNodes();++n)
		{
			NuTo::NodeBase* node=myGrid.NodeGetNodePtr(n);
			dofs=node->GetGlobalDofs();
		}
		*/

		myGrid.CalculateVoxelLocations();
		myGrid.SetAllNodeIds();
		myGrid.SetAllElementIds();
		myGrid.SetAllNodeIdsAtNode();
		myGrid.SetAllPartCoefficientMatrix0();
#ifdef SHOW_TIME
	end=clock();
	std::cout<<"[NuTo::StructureGrid3D] structure set " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

		NuTo::CallbackHandlerGrid myCallback;

		std::cout<<__FILE__<<" "<<__LINE__<<"  callback created"<<std::endl;
		NuTo::ConjugateGradientStructureGrid myOptimizer((unsigned int) myGrid.GetNumNodes()*3);
		std::cout<<__FILE__<<" "<<__LINE__<<"  optimizer created"<<std::endl;

		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> startVector(myGrid.GetNumNodes()*3,1);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rDisplacements(3,1);

		for (int myNode=0; myNode<numNodes;++myNode)
		{
			myGrid.NodeGetDisplacements(myNode,rDisplacements);
			startVector.SetBlock(3*myNode,0,rDisplacements);
		}
		std::cout<<__FILE__<<" "<<__LINE__<<"  "<<std::endl;
		myOptimizer.SetVerboseLevel(5);

		myOptimizer.SetParameters(startVector);
		std::cout<<__FILE__<<" "<<__LINE__<<"  Parameters set, Anzahl "<<myOptimizer.GetNumParameters()<<std::endl;
#ifdef ENABLE_MECHANICS
		myOptimizer.SetGridStructure(&myGrid);
 #endif // ENABLE_MECHANICS

		std::cout<<__FILE__<<" "<<__LINE__<<"  Grid set"<<std::endl;

		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> returnVector(myGrid.GetNumNodes()*3,1);
		std::cout<<__FILE__<<" "<<__LINE__<<" startVector filled, last value"<<startVector(myGrid.GetNumNodes()*3-1,0)<<std::endl;
		//set callback routines for the calculation of the objective function, gradient etc
		//this works, because Neural network has been derived from CallbackHandler of the optimization module
		myOptimizer.SetCallback(dynamic_cast<NuTo::CallbackHandler*>(&myCallback));
		myOptimizer.Optimize();
		std::cout<<__FILE__<<" "<<__LINE__<<"  optimiert"<<std::endl;
		return 0;
		returnVector=myOptimizer.GetParameters();
		//
		// myGrid.NodeMergeActiveDofValues(returnVector);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispAll(0,1);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispNode(3,1);
		/*
		for (int count=0; count<myGrid.GetNumNodes(); ++count)
		{
			myGrid.NodeGetDisplacements(count,dispNode);
			dispAll.AppendRows(dispNode);
		}
		*/
		//dispAll=returnVector;
		//dispAll.WriteToFile("displacements.txt"," ");
		returnVector.WriteToFile("displacements.txt"," ");
		//dispAll.Info(12,6);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispAnsys(myGrid.GetNumNodes(),1);
		//dispAnsys.ReadFromFile("disp_cube11_ansys.txt");
		/*
		try
		{
			dispAnsys.ReadFromFile("solu_disp_checkerboard_cube11.txt");
			//NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispAllDiff=dispAll-dispAnsys;
			dispAllDiff.WriteToFile("displDiff.txt"," ");
			Eigen::VectorXd diffVector;
			Eigen::VectorXd dispAllVector;
			diffVector=dispAllDiff.mEigenMatrix;
			dispAllVector=dispAll.mEigenMatrix;
			std::cout<<__FILE__<<" "<<__LINE__<<" diff norm " <<diffVector.norm()<<std::endl;
			std::cout<<__FILE__<<" "<<__LINE__<<" error " <<diffVector.norm()/dispAllVector.norm()*100<<" %"<<std::endl;

			//returnVector.WriteToFile("ActiveDofdisplacements.txt"," ");
		}
		catch(NuTo::MechanicsException& e)
		{}
		*/
		//std::cout<<__FILE__<<__LINE__<<"teste vorhandene Objekte \n"<<
		// myGrid.m
		//<<std::endl;
		return 0;
		/*
#ifdef ENABLE_VISUALIZE
		// visualize element
		myGrid.AddVisualizationComponentDisplacements();
		myGrid.AddVisualizationComponentEngineeringStrain();
		myGrid.AddVisualizationComponentEngineeringStress();
		myGrid.ExportVtkDataFile("Grid3D.vtk");
#endif //ENABLE_VISUALIZE
*/
		std::cout<<"numpar "<<myOptimizer.GetNumParameters()<<std::endl;

		//visualize results
		//myGrid.ExportVtkDataFile("StrukturedGrid3D.vtk","DISPLACEMENTS");


	}
	catch (NuTo::Exception& e)
	{
		std::cout << e.ErrorMessage() << std::endl;
	}
	return 0;
}
