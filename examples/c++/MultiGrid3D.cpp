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
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/grid/MultiGridStructure.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/Jacobi.h"

int main()
{
	bool matrixFreeMethod=0; //0 -EBE, 1- NBN, false=0
//	bool matrixFreeMethod=1; //0 -EBE, 1- NBN

//    double PoissonsRatio = 0.;
    double PoissonsRatio = 0.2;

    // create structure
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    // read entries
	NuTo::StructureGrid myGrid(3); // also creates CallbackHandler
	myGrid.SetVerboseLevel(0);
//	std::string inputFile="InputTest";
	std::string inputFile="../trunk/examples/c++/InputStructureGrid3D";
	myGrid.ImportFromVtkASCIIFileHeader(inputFile);
	//calculate one element stiffness matrix with E=1

	myGrid.SetMatrixFreeMethod(matrixFreeMethod);
	myGrid.SetBasisElementStiffnessMatrix(PoissonsRatio,0);

	//grid structure create
	/**************************************************************************/
	//material values are smaller than threshold value
	//color values form 0 255
	//thresholdMaterialValues: 180, 188 ->last value with material
	int thresholdMaterialValue=188; //last value for "air", voxel which does not get element
	std::vector<double> myMapColorModul(256);

	//set Modul for each color
	double MaterialYoungsModulus=100000.;
	for(int count=0;count<thresholdMaterialValue;count++)
//		myMapColorModul[count]=1.;
		myMapColorModul[count]=MaterialYoungsModulus;
	for(int count=thresholdMaterialValue;count<255;count++)
		myMapColorModul[count]=0.;

	// create grid structure

	myGrid.CreateGrid(thresholdMaterialValue,inputFile,myMapColorModul);
	if (matrixFreeMethod)
	{
		myGrid.SetBasisEdgeStiffnessMatrices(0);
		myGrid.SetNeighborNodesNE();
		myGrid.SetMaterialNumberForEdges();
	}

	//----------------------------------------------------------------------------------------//
	// Boundary condition: all nodes with z=0
	// Boundary condition: set x,y,z direction zero
	//----------------------------------------------------------------------------------------//

////	bool EnableDisplacementControl = false;
//	double BoundaryDisplacement = -1.0;
////	double BoundaryDisplacement = -(rGridDimension[2]*rVoxelSpacing[2]*microTomm)/20.0;
//	std::cout<<"[MultiGrid3D]  Boundary Displacement: "<<BoundaryDisplacement<<std::endl;
////	for z=0,x,y -all ux=0
//	const std::vector<size_t> rGridDimension=myGrid.GetGridDimension();
//	size_t direction=0;
//	// Attention: consider frame nodes
//	size_t rGridLocation[6]={1,rGridDimension[0]-1,1,rGridDimension[1]-1,1,1};
//	double rValue=0;
//	// diplacement vector plus one dof for calculation with non existing neighbor dofs
//	std::vector<double> rDisplVector(3*(3*myGrid.GetNumNodes()+1),0.0);// initialized with zero
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);
//
//	//for z=0,x,y -all uy=0
//	direction=1;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);
//
//	//for z=0,x,y -all uz=0
//	direction=2;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);
//
//	//for z=zmax,x,y -all uz=-1
//	direction=2;
//	rGridLocation[4]=rGridDimension[2]-1;
//	rGridLocation[5]=rGridDimension[2]-1;
//	rValue=BoundaryDisplacement;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//----------------------------------------------------------------------------------------//
	// Boundary condition: x=0, ux=0; y=0, uy=0; z=0,uz=0;z0max,uz=-1
	//----------------------------------------------------------------------------------------//

	double BoundaryDisplacement = -1.0;
//	double BoundaryDisplacement = -(myGrid.GetGridDimension()[2]*myGrid.GetVoxelSpacing()[2])/20.0;
	// diplacement vector plus one dof for calculation with non existing neighbor dofs
	std::vector<double> rDisplVector(3*(3*myGrid.GetNumNodes()+1),0.0);// initialized with zero
	//for z=0,x,y -all ux=0
	const std::vector<size_t> rGridDimension=myGrid.GetGridDimension();
	// Attention: consider frame nodes
	double rValue=0.;
	// rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
	//x=0, ux=0
	size_t direction=0;
	size_t rGridLocation[6]={1,1,1,rGridDimension[1],1,rGridDimension[2]};
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//y=0; uy=0
	direction=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[3]=1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	// z=0; uz=0;
	direction=2;
	rGridLocation[3]=rGridDimension[1];
	rGridLocation[5]=1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//z=zmax, uz=-1
	rValue=BoundaryDisplacement;
	rGridLocation[4]=rGridDimension[2]-1;
	rGridLocation[5]=rGridDimension[2]-1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);


//	myGrid.AnsysInput(rDisplVector);


	size_t numDofs=myGrid.GetNumNodes()*3;
#ifdef SHOW_TIME
end=clock();
#endif
	double mem=sizeof(myGrid);
	myGrid.Info();
	std::cout<<"-   Poissons ratio .............................. "<<PoissonsRatio<<"\n\n";
	std::cout<<"Boundary displacement .......................... "<<BoundaryDisplacement<<"\n";
	std::cout<<"Direction of external force/displacement ....... z\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
#ifdef SHOW_TIME
	std::cout<<"Time for structure initialization............... "<<difftime(end,start)/CLOCKS_PER_SEC<<"(sec)\n";
#endif
	std::cout<<"Allocated memory ............................... "<<mem/1000.<<"(MB)\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";

	std::cout<<"[MultiGrid3D] number of dofs "<<numDofs<<" free: "<<numDofs-myGrid.GetNumConstraints()<<" constraint: "<<myGrid.GetNumConstraints()<<"\n";
	// start analysis
//	std::cout<<__FILE__<<" "<<__LINE__<<"  start analysis"<<std::endl;
	myGrid.SetVerboseLevel(0);
	myGrid.SetMisesWielandt(false);
	size_t numNodes=myGrid.GetNumNodes();

	NuTo::MultiGridStructure myMultiGrid;
	myMultiGrid.SetVerboseLevel(3);
	myMultiGrid.SetStructure(&myGrid);
	myMultiGrid.SetUseMultiGridAsPreconditoner(false);
	myMultiGrid.SetMaxCycle(numNodes);
	myMultiGrid.SetNumPreSmoothingSteps(1);
	myMultiGrid.SetNumPostSmoothingSteps(1);
	// for test purpose of mg
//	std::ifstream test;
//	test.open("result.txt");
//	if(test)	// file is open
//	{
//		rDisplVector.clear();
//		double help=0;
//		size_t count=0;
//		while(!test.eof()) // keep reading untill end-of-file
//		{
//			test>>help;
//			rDisplVector.push_back(help);
//			++count;
//		}
//		test.close();
//		--count; // for last empty line
//
//		if (count!=numDofs)
//			std::cout<<"[MultiGrid3D] Wrong input file.\n";
//	rDisplVector.push_back(0.0);
//	rDisplVector.push_back(0.0);
//	rDisplVector.push_back(0.0);
//	myMultiGrid.SetParameters(rDisplVector);
//	}

	myMultiGrid.Initialize();
	myMultiGrid.Info();
//	std::vector<double> residual(numNodes*3);
//	myGrid.Gradient(rDisplVector,residual);
//	rDisplVector.assign((numNodes+1)*3, 0.0);
//	myGrid.SetParameters(rDisplVector);
//	for(size_t i=0;i<numNodes*3;++i)
//		residual[i]*=-1;
//	myGrid.SetResidual(residual);
////	myMultiGrid.ExportVTKStructuredDataFile(0,"outputGrid0.vtk");
////	myMultiGrid.ExportVTKStructuredDataFile(1,"outputGrid1.vtk");
//
//	std::vector<double> error;
//	std::ifstream input2;
//	double help2=0;
//	// result file only with existing nodes
//	input2.open("result.txt");
//	if(input2)	// file is open
//	{
//		size_t count=0;
//		while(!input2.eof()) // keep reading untill end-of-file
//		{
//			input2>>help2;
//			error.push_back(help2);
//			++count;
//		}
//		input2.close();
//		--count; // for last empty line
//	}
//	myGrid.Gradient(error,residual);
//	std::cout<<"MG: residual of error ";
//	for(size_t i=0;i<numNodes*3;++i)
//		std::cout<<" "<<residual[i];
//	std::cout<<"\n";
//

	myMultiGrid.MultiGridSolve();

// test jacobi
//	NuTo::Jacobi myOptimizer(numNodes*3);
//	myOptimizer.SetVerboseLevel(5);
//	myOptimizer.SetCallback(&myGrid);
//	myOptimizer.Optimize();

//	rDisplVector=myMultiGrid.GetParameters();
	rDisplVector=myGrid.GetParameters();

	std::ofstream file;
    file.open("displacementsMG.txt");
	for(size_t i=0;i<numNodes;++i)
	{
			file<<rDisplVector[3*i]<<"\n";
			file<<rDisplVector[3*i+1]<<"\n";
			file<<rDisplVector[3*i+2]<<"\n";
	}
	file.close();

	file.open("displVTK.txt");
	size_t numGridNodes=(myGrid.GetGridDimension()[0]+1)*(myGrid.GetGridDimension()[1]+1)*(myGrid.GetGridDimension()[2]+1);
	for(size_t i=0;i<numGridNodes;++i)
	{
		size_t nodeId=myGrid.GetNodeId(i);
		if (nodeId==(size_t) myGrid.GetNumNodes())
		{
			file<<0.0<<"\n";
			file<<0.0<<"\n";
			file<<0.0<<"\n";
		}
		else
		{
			file<<rDisplVector[3*nodeId]<<"\n";
			file<<rDisplVector[3*nodeId+1]<<"\n";
			file<<rDisplVector[3*nodeId+2]<<"\n";
		}
	}
	file.close();

	std::vector<double> dispRef;
	std::ifstream input;
	double help=0;
	// result file only with existing nodes
	input.open("result.txt");
	if(input)	// file is open
	{
		size_t count=0;
		while(!input.eof()) // keep reading untill end-of-file
		{
			input>>help;
			dispRef.push_back(help);
			++count;
		}
		input.close();
		--count; // for last empty line

		if (count==numDofs)
		{
			double squareDiffNorm=0;
			double squareRefNorm=0;
			std::ofstream diffFile;
			diffFile.open("displDiff.txt");
			for(size_t i=0;i<numNodes;++i)
			{
					diffFile<<rDisplVector[3*i]-dispRef[3*i]<<"\n";
					diffFile<<rDisplVector[3*i+1]-dispRef[3*i+1]<<"\n";
					diffFile<<rDisplVector[3*i+2]-dispRef[3*i+2]<<"\n";
				squareDiffNorm+=(rDisplVector[3*i]-dispRef[3*i])*(rDisplVector[3*i]-dispRef[3*i]);
				squareDiffNorm+=(rDisplVector[3*i+1]-dispRef[3*i+1])*(rDisplVector[3*i+1]-dispRef[3*i+1]);
				squareDiffNorm+=(rDisplVector[3*i+2]-dispRef[3*i+2])*(rDisplVector[3*i+2]-dispRef[3*i+2]);
				squareRefNorm+=(dispRef[3*i])*(dispRef[3*i]);
				squareRefNorm+=(dispRef[3*i+1])*(dispRef[3*i+1]);
				squareRefNorm+=(dispRef[3*i+2])*(dispRef[3*i+2]);
			}
		std::cout<<"[MultiGrid3D] squared diff norm " <<squareDiffNorm<<std::endl;
		std::cout<<"[MultiGrid3D] error " <<sqrt(squareDiffNorm)/sqrt(squareRefNorm)*100<<" %"<<std::endl;
		}
		else
			std::cout<<"[MultiGrid3D] Comparison with reference results is not possible (wrong size).\n";
	}
	else
		std::cout<<"[MultiGrid3D] Comparison with reference results is not possible (no result file).\n";

	return 0;
}
