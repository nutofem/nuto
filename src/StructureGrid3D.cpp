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


#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/ConjugateGradientLinear.h"

int main()
{
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
        NuTo::StructureGrid myGrid(3);

        // read entries
        myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("/home/fuhlrott/develop/nuto_work/nuto/examples/c++/InputStructureGrid3D");
        numVoxel=myGrid.GetNumVoxels();
        voxelSpacing=myGrid.GetVoxelSpacing();
        gridOrigin=myGrid.GetGridOrigin();
        gridDimension=myGrid.GetGridDimension();
        std::cout<<__FILE__<<"  variab: spac: "<<voxelSpacing[0]<< " gridDim: "<<gridDimension[0]<<std::endl;
        NuTo::FullMatrix<int> imageValues (numVoxel,1);
        imageValues.NuTo::FullMatrix<int>::ImportFromVtkASCIIFile( "/home/fuhlrott/develop/nuto_work/nuto/examples/c++/InputStructureGrid3D");

        std::cout<<"first value "<< imageValues(0,0) << std::endl;
        std::cout<<"numVoxel"<< numVoxel << std::endl;

        //RB
        //double Force = 1.;
        bool EnableDisplacementControl = true;
        double BoundaryDisplacement = 50;

        //calculate one element stiffness matrix with E=1
        NuTo::Structure myHelpStruc(3);
        // create material law
        int myMat=myHelpStruc.ConstitutiveLawCreate("LinearElastic");
        myHelpStruc.ConstitutiveLawSetPoissonsRatio(myMat, PoissonsRatio);
        myHelpStruc.ConstitutiveLawSetYoungsModulus(myMat, YoungsModulus);

        // create nodes
        NuTo::FullMatrix<double> nodeCoordinates(3, 1);
        NuTo::FullMatrix<int> elementIncidence(8,1);
        int count = 0;

         // for first voxel

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(0,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(1,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(2,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(3,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(4,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(5,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(6,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(7,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        // first element create

        elementIncidence.Info();
        int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
        //myHelpStruc.ElementSetConstitutiveLaw(element,"Material1");
        myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);

        myHelpStruc.NodeInfo(0);

        myHelpStruc.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;
        myHelpStruc.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);

        stiffnessMatrix.RemoveZeroEntries(0,1e-14);
         stiffnessMatrix.Info();
        //NuTo::FullMatrix A(stiffnessMatrix);
        //A.WriteToFile("$HOME/develop/nuto/stiffnessMatrix.txt"," ");

        //grid structure create
        myGrid.CreateNodeGrid("DISPLACEMENTS");

        //myGrid.CreateNodeGrid("COORDINATES");
        std::cout<<"NodeGrid created"<<std::endl;
        int numNodes=myGrid.GetNumNodes();
        std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;
        NuTo::NodeBase* myNode=myGrid.NodeGetNodePtr(numNodes-1);
        std::cout<<"knoten "<<myNode->GetNodeGridNum()<<std::endl; //Gitterknoten


        //set Modul for each color
        NuTo::FullMatrix<double> myMapColorModul(255,1);
        count=0;
        for(count=0;count<131;count++)
//            myMapColorModul(count,0)=0;
            myMapColorModul(count,0)=1000;
        for(count=131;count<161;count++)
            //myMapColorModul(count,0)=8300.;
            myMapColorModul(count,0)=1000.;
        for (count=161;count<255;count++)
            //myMapColorModul(count,0)=11500.;
            myMapColorModul(count,0)=1000.;

        //myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");
        //generiert Knoten, Freiheitsgrade noch nicht gesetzt


         myGrid.CreateElementGrid(stiffnessMatrix,myMapColorModul,"VOXEL8N");
        std::cout<<"ElementGrid created"<<std::endl;

        // boundary conditions
        int NumElementsX = (int) gridDimension[0];
        int NumElementsY = (int) gridDimension[1];
        int NumElementsZ = (int) gridDimension[2];

        NuTo::FullMatrix<double> direction(3,1);
        direction(0,0)= 1;
        direction(1,0)= 0;
        direction(2,0)= 0;
        std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;

        int myNodeNumber=0;
		for (int count = 0;count<myGrid.GetNumNodes();count= count + (NumElementsX + 1))
		{
			std::cout<<__FILE__<<" "<<__LINE__<<" node constraint "<< count <<std::endl;
		   try
			{
				myNodeNumber=myGrid.NodeGetIdFromGridNum(count); //node from type NodeGridCoordinates
				myGrid.ConstraintLinearSetDisplacementNode(myNodeNumber, direction, 0.0);
			}
			catch(NuTo::MechanicsException& e)
			{}
        }

		direction(0,0)= 0;
        direction(1,0)= 0;
        direction(2,0)= 1;
//		for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1);++count)
		for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);++count)
		{
            std::cout<<__FILE__<<" "<<__LINE__<<" node constraint z "<< count <<std::endl;
			try
			{
				myNodeNumber=myGrid.NodeGetIdFromGridNum(count); //node from type NodeGridCoordinates
				myGrid.ConstraintLinearSetDisplacementNode(myNodeNumber, direction, 0.0);
			}
			catch(NuTo::MechanicsException& e)
			{}

		}

		direction(0,0)= 0;
        direction(1,0)= 1;
        direction(2,0)= 0;
        for (int countY = 0; countY < (NumElementsY); ++countY)
        {
//			for (int count = 0;count<(NumElementsX + 1);++count)
			for (int count = 0;count<(NumElementsX + 1)*(NumElementsY+1);++count)
			{
				int node = count+countY*(NumElementsX + 1)*(NumElementsY + 1);
				std::cout<<__FILE__<<" "<<__LINE__<<" node constraint y "<< node <<std::endl;
				try
				{
					myNodeNumber=myGrid.NodeGetIdFromGridNum(node); //node from type NodeGridCoordinates
					myGrid.ConstraintLinearSetDisplacementNode(myNodeNumber, direction, 0.0);
				}
				catch(NuTo::MechanicsException& e)
				{}
			}

        }

        // apply nodes
        if(EnableDisplacementControl)
        {
            std::cout << "Displacement control" << std::endl;
            // boundary displacments
            direction(0,0)= 1;
            direction(1,0)= 0;
            direction(2,0)= 0;
            NuTo::FullMatrix<double> displacements(3,1);
            displacements(0,0)= BoundaryDisplacement;
            displacements(1,0)= 0;
            displacements(2,0)= 0;

             for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
            {
                //std::cout << zCount << std::endl;
                for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
                {
                    int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    //std::cout << "node: " << node << std::endl;
                    int myNodeNumber;
                    int flag=0;
                    try
                    {
                        myNodeNumber=myGrid.NodeGetIdFromGridNum(node); //node from type NodeGridCoordinates
                    }
                    catch(NuTo::MechanicsException& e)
                    {
                        flag=1;
                    }
                   if(flag==0)
                   {
						myGrid.NodeSetDisplacements(myNodeNumber, displacements);
                        myGrid.ConstraintLinearSetDisplacementNode(myNodeNumber, direction, BoundaryDisplacement);
                   }
                }
            }
         }
        else
        {
            std::cout <<__FILE__<<" "<<__LINE__<< "Load control" << "not implemented"<<std::endl;
         /*   //! @TODO: Add special configurations if edge node not exist
            // apply load to nodes
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
        myGrid.SetVerboseLevel(2);
        myGrid.NodeBuildGlobalDofs();

        std::cout<<__FILE__<<" "<<__LINE__<<"  glob dofs "<<myGrid.GetNumDofs()<<std::endl;
        std::cout<<__FILE__<<" "<<__LINE__<<" active dofs "<<myGrid.GetNumActiveDofs()<<std::endl;

        NuTo::FullMatrix<int> voxelLocation(myGrid.GetNumElements(),4);
        int *dofs;
        for (int n=0;n< myGrid.GetNumNodes();++n)
        {
        	NuTo::NodeBase* node=myGrid.NodeGetNodePtr(n);
			dofs=node->GetGlobalDofs();
        }

         myGrid.CalculateVoxelNumAndLocMatrix(voxelLocation);
	     //std::cout<<__FILE__<<" "<<__LINE__<<" VoxelLocationmatrix: \n"<<voxelLocation<< "\n\n"<<std::endl;

       NuTo::CallbackHandlerGrid myCallback;
        std::cout<<__FILE__<<" "<<__LINE__<<"  callback crated"<<std::endl;
        NuTo::ConjugateGradientLinear myOptimizer((unsigned int) myGrid.GetNumActiveDofs());
        std::cout<<__FILE__<<" "<<__LINE__<<"  optimizer created"<<std::endl;
        NuTo::FullMatrix<double> startVector(myGrid.GetNumActiveDofs(),1);
        count=1;
        for(int ii=0;ii<myGrid.GetNumActiveDofs(); ii++)
            startVector(ii,0)=0;

        myOptimizer.SetVerboseLevel(3);
        myOptimizer.SetParameters(startVector);
        std::cout<<__FILE__<<" "<<__LINE__<<"  Parameters set"<<std::endl;
#ifdef ENABLE_MECHANICS
        myOptimizer.SetGridStructure(&myGrid);
//#else
  //      myOptimizer.SetGridStructure();
#endif // ENABLE_MECHANICS

        std::cout<<__FILE__<<" "<<__LINE__<<"  Grid set"<<std::endl;

        NuTo::FullMatrix<double> returnVector(myGrid.GetNumActiveDofs(),1);
        std::cout<<__FILE__<<" "<<__LINE__<<" startVector filled, last value"<<startVector(myGrid.GetNumActiveDofs()-1,0)<<std::endl;
        //set callback routines for the calculation of the objective function, gradient etc
        //this works, because Neural network has been derived from CallbackHandler of the optimization module
        myOptimizer.SetCallback(dynamic_cast<NuTo::CallbackHandler*>(&myCallback));
        myOptimizer.Optimize();
        std::cout<<__FILE__<<" "<<__LINE__<<"  optimiert"<<std::endl;
        // visualize element
		return 0;
        myGrid.AddVisualizationComponentDisplacements();
		myGrid.AddVisualizationComponentEngineeringStrain();
		myGrid.AddVisualizationComponentEngineeringStress();
		myGrid.ExportVtkDataFile("Grid3D.vtk");

        std::cout<<"numpar "<<myOptimizer.GetNumParameters()<<std::endl;
        /*
        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> globStiffnessMatrix;
        myGrid.BuildGlobalCoefficientMatrix0(globStiffnessMatrix, dispForceVector);
        globStiffnessMatrix.RemoveZeroEntries(0,1e-14);
        //NuTo::FullMatrix<double> A(stiffnessMatrix);
        //A.WriteToFile("stiffnessMatrix.txt"," ");
        //stiffnessMatrix.Info();
        dispForceVector.Info();

        // build global external load vector
        NuTo::FullMatrix<double> extForceVector;
        myGrid.BuildGlobalExternalLoadVector(extForceVector);
        //extForceVector.Info();

        // calculate right hand side
        NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector;
        rhsVector.WriteToFile("rhsVector.txt"," ");

        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        NuTo::FullMatrix<double> displacementVector;
        stiffnessMatrix.SetOneBasedIndexing();
        mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
        displacementVector.WriteToFile("displacementVector.txt"," ");

        // write displacements to node
        myGrid.NodeMergeActiveDofValues(displacementVector);

        // calculate residual
        NuTo::FullMatrix<double> intForceVector;
        myGrid.BuildGlobalGradientInternalPotentialVector(intForceVector);
        NuTo::FullMatrix<double> residualVector = extForceVector - intForceVector;
        std::cout << "residual: " << residualVector.Norm() << std::endl;
*/
        // visualize results
//        myGrid.ExportVtkDataFile("StrukturedGrid3D.vtk","DISPLACEMENTS");


    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
