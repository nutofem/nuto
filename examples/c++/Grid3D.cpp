// $Id$
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
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#ifdef ENABLE_OPTIMIZE
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/Jacobi.h"
#include "nuto/optimize/MisesWielandt.h"
#endif

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
////	double BoundaryDisplacement = -(myGrid.GetGridDimension()[2]*myGrid.GetVoxelSpacing()[2])/20.0;
//	std::cout<<"[NuTo::Grid3D]  Boundary Displacement: "<<BoundaryDisplacement<<std::endl;
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
//
	// end //

	size_t numDofs=myGrid.GetNumNodes()*3;
	size_t numNodes=myGrid.GetNumNodes();
	//----------------------------------------------------------------------------------------//
	// Boundary condition: x=0, ux=0; y=0, uy=0; z=0,uz=0;uymax,uy=-1
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

	//  uz=0; for one plain
	direction=2;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=1;
	rGridLocation[3]=rGridDimension[1];
	rGridLocation[4]=1;
	rGridLocation[5]=1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);


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

/*
 *
	//new boundary conditions min. free  dofs
	//-------------------------------------------------------------------------------------//
	//	bool EnableDisplacementControl = false;
	double BoundaryDisplacement = 1.0;
	// diplacement vector plus one dof for calculation with non existing neighbor dofs
	std::vector<double> rDisplVector(3*(3*numNodes+1),0.0);// initialized with zero
	//for z=0,x,y -all ux=0
	const std::vector<size_t> rGridDimension=myGrid.GetGridDimension();
	// Attention: no considering of frame nodes for octree
	double rValue=0.;
	// rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
	size_t rGridLocation[6]={0};

	//ux=0 all x
	size_t direction=0;
	rGridLocation[0]=0;
	rGridLocation[1]=rGridDimension[0]+1;
	rGridLocation[2]=0;
	rGridLocation[3]=rGridDimension[1]+1;
	rGridLocation[4]=0;
	rGridLocation[5]=rGridDimension[2]+1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//y=0; uy=0
	direction=1;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=1;
	rGridLocation[3]=3;
	rGridLocation[4]=1;
	rGridLocation[5]=rGridDimension[2];

	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//  uz=0 all z
	direction=2;
	rGridLocation[0]=0;
	rGridLocation[1]=rGridDimension[0]+1;
	rGridLocation[2]=0;
	rGridLocation[3]=rGridDimension[1]+1;
	rGridLocation[4]=0;
	rGridLocation[5]=rGridDimension[2]+1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

//y=max, uy=-1
	direction=1;
	rGridLocation[0]=1;
	rGridLocation[1]=rGridDimension[0];
	rGridLocation[2]=rGridDimension[1]-3;
	rGridLocation[3]=rGridDimension[1];
	rGridLocation[4]=1;
	rGridLocation[5]=rGridDimension[2];
	rValue=BoundaryDisplacement;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,1,rDisplVector);
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,0,rDisplVector);
*/
//	myGrid.AnsysInput(rDisplVector);


//	std::cout<<"constraints: ";
//	for(size_t i=0;i<numDofs;++i)
//		std::cout<<myGrid.GetDisplacementConstaints()[i]<<" ";
//	std::cout<<"\n";

#ifdef SHOW_TIME
end=clock();
#endif
	double mem=sizeof(myGrid);
	myGrid.Info();
	std::cout<<"-   Poissons ratio .............................. "<<PoissonsRatio<<"\n\n";
	std::cout<<"Boundary displacement .......................... "<<BoundaryDisplacement<<"\n";
	std::cout<<"Direction of external force/displacement ....... z\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
	std::cout<<"Time for structure initialization............... "<<difftime(end,start)/CLOCKS_PER_SEC<<"(sec)\n";
	std::cout<<"Allocated memory ............................... "<<mem/1000.<<"(MB)\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";

#ifdef ENABLE_OPTIMIZE
	double accuracy=1e-6;
	int precision = 6;

	enum SolMethod // enum for solution method
	{
		JCG, 	// jacobi preconditioned conjugate gradient
		WJCG, 	// weighted jacobi preconditioned conjugate gradient
		EJCG,   // error equation - jacobi conjugate gradient
		J, 		// jacobi method
		EJ, 	//error equation - jacobi method
		M,		//mises method for max. eigenvalue
	} solMeth=JCG;

	std::ofstream file;
	file.open("sumOutput",std::ofstream::out|std::ofstream::app);
	if(file)
	{
		//output : voxels in one direction - dofs -
		// nbr grids - solMeth - nbr cycles -nbr pre -nbr post - time -its
		file<<rGridDimension[0]-2<<" "<<numDofs<<" 1 ";
		if(solMeth==JCG)
			file<<" JCG";
		else
			file<<"  nd";
		file<<" - - - ";
		file.close();
	}

	// start analysis
	if(solMeth==JCG)
	{
		std::cout<<"[Grid3D] Solution method jacobi preconditioned conjugate cradient \n";
		NuTo::ConjugateGradientGrid myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myGrid.SetMisesWielandt(false);
		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.SetAccuracyGradient(accuracy);
		myOptimizer.Info();

		myOptimizer.Optimize();
		rDisplVector=myGrid.GetParameters();
	}
	else if(solMeth==WJCG)
	{
		std::cout<<"[Grid3D] Solution method weighted jacobi preconditioned conjugate cradient \n";
		NuTo::ConjugateGradientGrid myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myGrid.SetMisesWielandt(true);
		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.Info();

		myOptimizer.Optimize();
		rDisplVector=myGrid.GetParameters();
	}
	else if (solMeth==EJCG)
	{
		std::cout<<"[Grid3D] Solution method jacobi conjugate cradient with error equation. \n";
		NuTo::ConjugateGradientGrid myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myGrid.SetMisesWielandt(false);
		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.Info();

		std::vector<double> residual(rDisplVector.size(),0.);
		std::vector<double> error(rDisplVector.size(),0.);
		myGrid.Gradient(rDisplVector,residual);
		for(size_t i=0;i<numNodes*3;++i)
			residual[i]*=-1;

		myGrid.SetParameters(error);
		myGrid.SetRightHandSide(residual);
		myOptimizer.Optimize();
		error=myGrid.GetParameters();
		for(size_t i=0;i<numNodes*3;++i)
			rDisplVector[i]+=error[i];

	}
	else if(solMeth==J)
	{
		std::cout<<"[Grid3D] Solution method jacobi. \n";
		//test jacobi solver
		NuTo::Jacobi myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.Info();
		myGrid.SetMisesWielandt(false);
		std::vector<double> residual(rDisplVector.size(),0.);
		// so
		myGrid.SetRightHandSide(residual);
		myOptimizer.SetParameters(rDisplVector);
		myOptimizer.Optimize();
		rDisplVector=myOptimizer.GetParametersVec();
		// or so
//		myOptimizer.Optimize(rDisplVector,residual);

	}
	else if(solMeth==EJ)
	{
		std::cout<<"[Grid3D] Solution method jacobi with error equation. \n";
		//test jacobi solver
		NuTo::Jacobi myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(1);
		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.Info();
		// res has address of grid residual
		std::vector<double> res((numNodes+1)*3, 0);
		std::vector<double> error((numNodes+1)*3, 0);
		myOptimizer.SetParameters(error);
		myGrid.Gradient(rDisplVector,res);
		for(size_t i=0;i<numNodes*3;++i)
			res[i]*=-1;
		// needed by Jacobi
		myGrid.SetRightHandSide(res);
		myGrid.SetMisesWielandt(false);
		myOptimizer.Optimize();
		error=myOptimizer.GetParametersVec();
		for(size_t i=0;i<numNodes*3;++i)
			rDisplVector[i]=error[i];
	}
	else if(solMeth==M)
	{
		// Convergenc test: lamda_max of M=I-PA <1
		myGrid.SetMisesWielandt(false); // if not, get a infinite loop, for Hessian
		NuTo::MisesWielandt myEigenCalculator(numNodes*3);
		myEigenCalculator.SetVerboseLevel(1);
		myEigenCalculator.SetAccuracyGradient(1e-6);
		myGrid.SetWeightingFactor(1.);
		myEigenCalculator.SetCallback((&myGrid));
		myEigenCalculator.Optimize();
		double lambda_max=myEigenCalculator.GetObjective();
		int precision = 10;
		std::cout.precision(precision);
		std::cout<<"Max. eigenvalue of preconditioned matrix "<<lambda_max<<"\n";
	}
#else //ENABLE_OPTIMIZE
	std::cout<<"[Grid3D] Solution is not possible. Module optimize is not loaded.\n";
#endif //ENABLE_OPTIMIZE

//	std::vector<double> rStrainVector;
//	myGrid.GetEngineeringStrain(rDisplVector, rStrainVector);
//	std::vector<double> rStressVector;
//	myGrid.GetEngineeringStress(rStrainVector,rStressVector);
//	double maxValue=0.;
//	for(size_t i=4;i<numDofs;i+=69)
//	{
//		if(maxValue<rStressVector[i])
//			maxValue=rStressVector[i];
//		else if (maxValue<-1*rStressVector[i])
//			maxValue=-rStressVector[i];
//	std::cout<<"Value of sigma yy ="<<maxValue<<"\n";
//
//	}
//	std::cout<<"maxValue of sigma yy ="<<maxValue<<"\n";
//	size_t dof=0;
//	double stress=0.;
//	int count=0;
//	for(size_t z=1;z<rGridDimension[2]+1;++z)
//	{
//		for(size_t i=rGridDimension[0]*rGridDimension[1]*z+1;i<rGridDimension[0]*rGridDimension[1]*z+rGridDimension[0]+1;i++)
//		{
//			dof=9*myGrid.GetNodeId(i)+4;
//			stress+=rStressVector[dof];
//			++count;
//		}
//	}
//	stress/=count;
//	std::cout<<"moyen value of sigma 0 "<<stress<<"\n";

	myGrid.ExportVTKStructuredDataFile("./outputFile.vtk");

//	std::ofstream file;
	file.precision(precision);
	file.open("displacements.txt");
	for(size_t i=0;i<numNodes;++i)
	{
			file<<rDisplVector[3*i]<<"\n";
			file<<rDisplVector[3*i+1]<<"\n";
			file<<rDisplVector[3*i+2]<<"\n";
	}
	file.close();

//    file.open("strains.txt");
//    size_t it_end=rStrainVector.size();
//	for(size_t i=0;i<it_end;++i)
//			file<<rStrainVector[i]<<"\n";
//	file.close();


	//
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
// output of diff and ref only for VTK
//			std::ofstream diffFile;
//			diffFile.open("displDiffVTK.txt");
//			file.open("displRefVTK.txt");
			for(size_t i=0;i<numNodes;++i)
			{
//					diffFile<<displVector[3*i]-dispRef[3*i]<<"\n";
//					diffFile<<displVector[3*i+1]-dispRef[3*i+1]<<"\n";
//					diffFile<<displVector[3*i+2]-dispRef[3*i+2]<<"\n";
//					file<<dispRef[3*i]<<"\n";
//					file<<dispRef[3*i+1]<<"\n";
//					file<<dispRef[3*i+2]<<"\n";
				squareDiffNorm+=(rDisplVector[3*i]-dispRef[3*i])*(rDisplVector[3*i]-dispRef[3*i]);
				squareDiffNorm+=(rDisplVector[3*i+1]-dispRef[3*i+1])*(rDisplVector[3*i+1]-dispRef[3*i+1]);
				squareDiffNorm+=(rDisplVector[3*i+2]-dispRef[3*i+2])*(rDisplVector[3*i+2]-dispRef[3*i+2]);
				squareRefNorm+=(dispRef[3*i])*(dispRef[3*i]);
				squareRefNorm+=(dispRef[3*i+1])*(dispRef[3*i+1]);
				squareRefNorm+=(dispRef[3*i+2])*(dispRef[3*i+2]);
			}
		std::cout<<"[NuTo::Grid3D] squared diff norm " <<squareDiffNorm<<std::endl;
		std::cout<<"[NuTo::Grid3D] error " <<sqrt(squareDiffNorm)/sqrt(squareRefNorm)*100<<" %"<<std::endl;
		}
		else
			std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (wrong size).\n";
	}
	else
		std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (no result file).\n";
	return 0;
}
