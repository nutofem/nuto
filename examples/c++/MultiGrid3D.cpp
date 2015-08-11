// $Id $
#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif

#include <iostream>
#include "stdlib.h"
#include "nuto/base/NuToObject.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/grid/MultiGridStructure.h"
#ifdef ENABLE_OPTIMIZE
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/Jacobi.h"
#include "nuto/optimize/ConjugateGradientGrid.h"
#endif //ENABLE_OPTIMIZE

// arguments for MGCG: nbr of cycles, nbr of pre , nbr of post smooting steps, nbr of grids,
int main(int argc, char *argv[])
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
    std::string inputFile="../../../nuto/examples/c++/InputStructureGrid3D";
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
	// Boundary condition: x=0, ux=0; y=0, uy=0; z=0,uz=0;z0max,uz=-1
	//----------------------------------------------------------------------------------------//

	//	bool EnableDisplacementControl = false;
	double BoundaryDisplacement = 1.0;
//	double BoundaryDisplacement = (double)(myGrid.GetGridDimension()[1])*0.01;
	// diplacement vector plus one dof for calculation with non existing neighbor dofs
	std::vector<double> rDisplVector(3*(3*myGrid.GetNumNodes()+1),0.0);// initialized with zero
	//for z=0,x,y -all ux=0
	const std::vector<size_t> rGridDimension=myGrid.GetGridDimension();
	// Attention: consider frame nodes
	double rValue=0.;
	// rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
	//x=0, ux=0
	size_t direction=0;
	// for symmetric
	size_t rGridLocation[6]={1,1,1,rGridDimension[1],1,rGridDimension[2]};
	// for hole plate with hole
//	size_t rGridLocation[6]={rGridDimension[0]/2,rGridDimension[0]/2,1,rGridDimension[1],1,rGridDimension[2]};
	// xmax
//	size_t rGridLocation[6]={rGridDimension[0]-1,rGridDimension[0],0,rGridDimension[1]+1,0,rGridDimension[2]+1};
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//y=0; uy=0
	direction=1;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=1;
	rGridLocation[3]=1;
	rGridLocation[4]=1;
	rGridLocation[5]=rGridDimension[2];
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,-1,rDisplVector);
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//  uz=0; for all nodes
//	direction=2;
//	rGridLocation[0]=1;
//	rGridLocation[1]=rGridDimension[0];
//	rGridLocation[2]=1;
//	rGridLocation[3]=rGridDimension[1];
//	rGridLocation[4]=1;
//	rGridLocation[5]=rGridDimension[2];
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);
//	std::cout<<"constraints: ";
//	for(size_t i=2;i<numDofs;i+=3)
//		std::cout<<myGrid.GetDisplacementConstaints()[i]<<" ";
//	std::cout<<"\n";


//	//  uz=0; for one nodes
//	direction=2;
//	rGridLocation[0]=rGridDimension[0]-1;
//	rGridLocation[1]=rGridDimension[0];
//	rGridLocation[2]=rGridDimension[1]-1;
//	rGridLocation[3]=rGridDimension[1];
//	rGridLocation[4]=1;
//	rGridLocation[5]=1;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//  uz=0; for one plain
	direction=2;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=1;
	rGridLocation[3]=rGridDimension[1];
	rGridLocation[4]=1;
	rGridLocation[5]=1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);



//	//z=zmax, uz=-1
//	rValue=BoundaryDisplacement;
//	rGridLocation[4]=rGridDimension[2]-1;
//	rGridLocation[5]=rGridDimension[2]-1;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//y=ymax, yz=-1
	direction=1;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=rGridDimension[1]-1;
	rGridLocation[3]=rGridDimension[1];
	rGridLocation[4]=1;
	rGridLocation[5]=rGridDimension[2];
	rValue=BoundaryDisplacement;
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,0,rDisplVector);
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

//	myGrid.AnsysInput(rDisplVector);


//	std::cout<<"constraints: ";
//	for(size_t i=0;i<numDofs;++i)
//		std::cout<<myGrid.GetDisplacementConstaints()[i]<<" ";
//	std::cout<<"\n";

//	myGrid.ExportVTKStructuredDataFile("multigridGeo.vtk");

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

	myGrid.SetVerboseLevel(1);
	myGrid.SetMisesWielandt(false);
	size_t numNodes=myGrid.GetNumNodes();

	NuTo::MultiGridStructure myMultiGrid;
	myMultiGrid.SetVerboseLevel(0);
	myMultiGrid.SetStructure(&myGrid);
	myMultiGrid.SetUseMultiGridAsPreconditoner(false);
	myMultiGrid.SetMaxCycle(numDofs/10);
	//added 01/2015
//	myMultiGrid.SetMaxGrids(3);
	myMultiGrid.SetMaxGrids(atoi(argv[4]));
	myMultiGrid.SetNumPreSmoothingSteps(1);
	myMultiGrid.SetNumPostSmoothingSteps(1);

	myMultiGrid.Initialize();

#ifdef ENABLE_OPTIMIZE
	enum SolMethod // enum for solution method
	{
		MG, 	// multigrid method
		EMG, 	//error equation -  multigrid method  -- does not work so far
		MGCG, 	//multigrid preconditioned conjugate gradient method
	} solMeth=MGCG;

	if(argc==1)//nothing added after programm name
	{
        solMeth=MG;
	}
	else
        solMeth=MGCG;

	if(solMeth==MG)
	{
		std::cout<<"[MultiGrid3D] Solution method is multigrid method. \n";
		std::vector<double> rhs(rDisplVector.size());
		myMultiGrid.SetVerboseLevel(1);
		myMultiGrid.Info();

		myMultiGrid.MultiGridSolve(rDisplVector,rhs);
	}
	else if(solMeth==EMG)
	{
		std::cout<<"[MultiGrid3D] Solution method is multigrid with error equation. \n";
		std::vector<double> residual(rDisplVector.size());
		std::vector<double> error(rDisplVector.size());

//		NuTo::Jacobi myOptimizer(numNodes*3);
//		myOptimizer.SetVerboseLevel(1);
//		myOptimizer.SetCallback( (&myGrid));
//		myOptimizer.Info();
//		myGrid.SetMisesWielandt(false);
//		myGrid.SetRightHandSide(residual);
//		myOptimizer.SetParameters(rDisplVector);
//		myOptimizer.SetMaxIterations(5);
//		myOptimizer.Optimize();
//		rDisplVector=myOptimizer.GetParametersVec();

		myGrid.Gradient(rDisplVector,residual);
		for(size_t i=0;i<numNodes*3;++i)
			residual[i]*=-1;
		std::cout<<"[MultiGrid3D] start residual: ";
		for(size_t i=0;i<numNodes;++i)
			std::cout<<residual[3*i+1]<<" ";
		std::cout<<"\n";
		myMultiGrid.SetVerboseLevel(1);
		myMultiGrid.Info();

		myMultiGrid.MultiGridSolve(error,residual);
		std::cout<<"[MultiGrid3D] solution: ";
		for(size_t i=0;i<numNodes*3;++i)
		{
			rDisplVector[i]+=error[i];
			std::cout<<rDisplVector[i]<<" ";
		}
		std::cout<<"\n";
	}
	else if(solMeth==MGCG)
	{
		std::cout<<"[MultiGrid3D] Solution method is multigrid preconditioned conjugate gradient method. \n";
		NuTo::ConjugateGradientGrid myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myOptimizer.SetMaxIterations((int) numNodes/500);
		myGrid.SetMisesWielandt(false);
		std::vector<double> rhs(rDisplVector.size());
		myGrid.SetParameters(rDisplVector);
		myGrid.SetRightHandSide(rhs);
		myMultiGrid.SetUseMultiGridAsPreconditoner(true);
		//2015-01: levels added
		if(argc==5)
		{
			myMultiGrid.SetMaxCycle(atoi(argv[1]));
			myMultiGrid.SetNumPreSmoothingSteps(atoi(argv[2]));
			myMultiGrid.SetNumPostSmoothingSteps(atoi(argv[3]));
		}
		else
		{
			myMultiGrid.SetMaxCycle(1);
			myMultiGrid.SetNumPreSmoothingSteps(1);
			myMultiGrid.SetNumPostSmoothingSteps(1);
		}
		myOptimizer.SetCallback( (&myMultiGrid));
		myMultiGrid.Info();
		myOptimizer.Info();

		std::ofstream file;
		file.open("sumOutput",std::ofstream::out|std::ofstream::app);
		if(file)
		{
			//output : voxels in one direction - dofs -
			// nbr grids - solMeth - nbr cycles -nbr pre -nbr post - time -its
			file<<rGridDimension[0]-2<<" "<<numDofs<<" "<<
					myMultiGrid.GetNumGrids()<<" ";
			if(solMeth==MG)
				file<<"  MG";
			else if(solMeth==EMG)
				file<<" EMG";
			else if(solMeth==MGCG)
				file<<"MGCG";
			file<<" "<<myMultiGrid.GetMaxCycle()<<" "
							""<<myMultiGrid.GetNumPreSmoothingSteps()<<" "<<myMultiGrid.GetNumPostSmoothingSteps()<<" ";
			file.close();
		}

		myOptimizer.Optimize();

		rDisplVector=myGrid.GetParameters();

		myMultiGrid.ExportVTKStructuredDataFile(0,"outputGrid0.vtk");
	//	myMultiGrid.ExportVTKStructuredDataFile(1,"outputGrid1.vtk");


	}
#else //ENABLE_OPTIMIZE
	std::cout<<"[MultiGrid3D] Solution is not possible. Module optimize is not loaded.\n";
#endif //ENABLE_OPTIMIZE

	std::ofstream file;
    file.open("displacementsMG.txt");
	for(size_t i=0;i<numNodes;++i)
	{
			file<<rDisplVector[3*i]<<"\n";
			file<<rDisplVector[3*i+1]<<"\n";
			file<<rDisplVector[3*i+2]<<"\n";
	}
	file.close();
//
//	file.open("displVTK.txt");
//	size_t numGridNodes=(myGrid.GetGridDimension()[0]+1)*(myGrid.GetGridDimension()[1]+1)*(myGrid.GetGridDimension()[2]+1);
//	for(size_t i=0;i<numGridNodes;++i)
//	{
//		size_t nodeId=myGrid.GetNodeId(i);
//		if (nodeId==(size_t) myGrid.GetNumNodes())
//		{
//			file<<0.0<<"\n";
//			file<<0.0<<"\n";
//			file<<0.0<<"\n";
//		}
//		else
//		{
//			file<<rDisplVector[3*nodeId]<<"\n";
//			file<<rDisplVector[3*nodeId+1]<<"\n";
//			file<<rDisplVector[3*nodeId+2]<<"\n";
//		}
//	}
//	file.close();

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
//			diffFile.open("displDiff.txt");
//			std::cout<<"[MultiGrid3D]  ref Solution ";
			for(size_t i=0;i<numNodes;++i)
			{
//					diffFile<<rDisplVector[3*i]-dispRef[3*i]<<"\n";
//					diffFile<<rDisplVector[3*i+1]-dispRef[3*i+1]<<"\n";
//					diffFile<<rDisplVector[3*i+2]-dispRef[3*i+2]<<"\n";
				squareDiffNorm+=(rDisplVector[3*i]-dispRef[3*i])*(rDisplVector[3*i]-dispRef[3*i]);
				squareDiffNorm+=(rDisplVector[3*i+1]-dispRef[3*i+1])*(rDisplVector[3*i+1]-dispRef[3*i+1]);
				squareDiffNorm+=(rDisplVector[3*i+2]-dispRef[3*i+2])*(rDisplVector[3*i+2]-dispRef[3*i+2]);
				squareRefNorm+=(dispRef[3*i])*(dispRef[3*i]);
				squareRefNorm+=(dispRef[3*i+1])*(dispRef[3*i+1]);
				squareRefNorm+=(dispRef[3*i+2])*(dispRef[3*i+2]);
//				std::cout<<dispRef[3*i]<<" "<<dispRef[3*i+1]<<" "<<dispRef[3*i+2]<<" ";
			}

		std::cout<<"\n[MultiGrid3D] squared diff norm " <<squareDiffNorm<<std::endl;
		std::cout<<"[MultiGrid3D] error " <<sqrt(squareDiffNorm)/sqrt(squareRefNorm)*100<<" %"<<std::endl;
		}
		else
			std::cout<<"[MultiGrid3D] Comparison with reference results is not possible (wrong size).\n";
	}
	else
		std::cout<<"[MultiGrid3D] Comparison with reference results is not possible (no result file).\n";

	return 0;
}
