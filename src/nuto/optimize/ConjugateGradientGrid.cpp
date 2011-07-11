// $Id$
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements3D.h"
#include "nuto/mechanics/elements/Voxel8N.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
#define tol 1e-8

int NuTo::ConjugateGradientGrid::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startOpt,endOpt;
    startOpt=clock();
#endif
	double alpha,
		   beta,
		   normGradient,
		   alphaNumerator=0,
		   alphaDenominator,
		   betaNumerator;

    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
		 curCycle(0);           //number of iterations without restart
    bool matrixFreeMethod=false; //false -EBE, true- NBN
    //bool matrixFreeMethod=true; //false -EBE, true- NBN
	if (mVerboseLevel>1)
	{
		std::cout<<__FILE__<<" "<<__LINE__<<" matrixFreeMethod is";
		if (matrixFreeMethod)
			std::cout<<" NBN "<<std::endl;
		else
			std::cout<<" EBE "<<std::endl;
	}


	optimization_return_attributes returnValue;

	FullMatrix<double> gradientOrig(GetNumParameters(),1);
	//FullMatrix<double> hessianOrig(GetNumParameters(),1);
	Eigen::VectorXd prevParameters;
	Eigen::VectorXd gradientScaled;
	Eigen::VectorXd scaleFactorsInv(GetNumParameters());
	Eigen::VectorXd prevGradientOrig(GetNumParameters());
	Eigen::VectorXd gradientNew(GetNumParameters());
	Eigen::VectorXd searchDirectionScaled(GetNumParameters());
	Eigen::VectorXd searchDirectionScaledDiff(GetNumParameters());
	Eigen::VectorXd searchDirectionOrig(GetNumParameters());
	Eigen::VectorXd gradientScaledSave(GetNumParameters());

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< " Para "<< GetNumParameters() << std::endl;
	bool converged(false);
	double mAccuracyGradientScaled = mAccuracyGradient;
	if (mVerboseLevel>3)
		std::cout<<__FILE__<<" "<<__LINE__<<" gradient accuracy "<<mAccuracyGradientScaled <<std::endl;

	int localMaxGradientCalls=2*GetNumParameters();
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);

	//check, if callback handler is set
	if (mpCallbackHandler==0)
		throw OptimizeException("[ConjugateGradientGrid::Optimize] Callback handler not set to determine objective function and derivatives.");

	// calculate objective
	numFunctionCalls++;
	if (numFunctionCalls>mMaxFunctionCalls)
	{
		converged = true;
		returnValue = MAXHESSIANCALLS;
	}

/*
 	if (matrixFreeMethod)
	{
		// save essential data before iterate
		int numDofs=3;
		int numNodes=mpGrid->GetNumNodes();
		double *globArray=new double [9*27*numNodes];
		int *globNodeIds=new int[27*numNodes];
		int *nodeIds=0;
		double *array=new double [9*27];
		bool *globConstraint = new bool[numNodes*numDofs];
		globConstraint = mpGrid->GetConstraintSwitch();
		//loop over all nodes
		for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
		{
			//get pointer to each node
			NodeGrid3D* thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
			//get belonging nodes
			nodeIds=thisNode->GetNodeIds();
			//get pointer to array of part coefficient matrix for all neighbor nodes
			try
			{
				array= thisNode->GetPartCoefficient0(0);
			}
			catch(MechanicsException &)
			{
				std::cout<<"ConjugateGradientGrid "<<__LINE__<<" nodeNum "<<nodeNumber<<std::endl;
				std::cout<<"neighbors: "<< nodeIds[0]<<", "<<nodeIds[1]<<", "<<nodeIds[2]<<", "<<nodeIds[3]<<", "<<nodeIds[4]<<", "<<nodeIds[5]<<", "<<nodeIds[6]<<", "<<nodeIds[7]<<", "<<nodeIds[8]<<", "<<nodeIds[9]<<", "<<nodeIds[10]<<", "<<nodeIds[11]<<", "<<nodeIds[12]<<", "<<nodeIds[13]<<", "<<nodeIds[14]<<", "<<nodeIds[15]<<", "<<nodeIds[16]<<", "<<nodeIds[17]<<", "<<nodeIds[18]<<", "<<nodeIds[19]<<", "<<nodeIds[20]<<", "<<nodeIds[21]<<", "<<nodeIds[22]<<", "<<nodeIds[23]<<", "<<nodeIds[24]<<", "<<nodeIds[25]<<", "<<nodeIds[26]<<std::endl;

			}
			for (int i=0;i<27;++i)
			{
				globNodeIds[27*nodeNumber+i]=nodeIds[i];
				globArray[27*9*nodeNumber+i]= array[i];
			}
		}
	}
*/
	//
	//

	while(!converged)
	{
		numGradientCalls++;
		 if (numGradientCalls>mMaxGradientCalls)
		 {
			 converged = true;
			 returnValue = MAXGRADIENTCALLS;
			 break;
			 //return MAXGRADIENTCALLS;
		 }
		//determine gradient

		// all  GetNumParameters() repeat start equations
		//if (curCycle%GetNumParameters()==0)
		// no update with start equations
		if (curCycle==0)
		{
			  // initialize search direction with steepest descent
			// calculate Hessian for scaling
 			//CalcScalingFactors(numHessianCalls,scaleFactorsInv);
/*

			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}
*/
			//calculate gradient as a start solution
			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" calc start direction"<<std::endl;

			if (!matrixFreeMethod)
				CalculateStartGradient(gradientOrig);
			/*
			FullMatrix<double> gradientOrigEBE(GetNumParameters(),1);
			gradientOrigEBE	=gradientOrig;
			gradientOrig.Resize(GetNumParameters(),1);
			std::cout<<__FILE__<<" "<<__LINE__<<" calc start direction NBN"<<std::endl;

*/
			if (matrixFreeMethod)
			{
				gradientOrig=mvParameters.mEigenMatrix;
				CalculateStartGradientNodeByNodeII(gradientOrig);
			}
				/*
			std::cout << __FILE__<<" "<<__LINE__<<"gradient EBE " <<std::endl ;
			(gradientOrigEBE.Trans()).Info();
			std::cout << __FILE__<<" "<<__LINE__<<"gradient NBN II " <<std::endl ;
			(gradientOrig.Trans()).Info();
			gradientOrigEBE.operator -=(gradientOrig);
			//std::cout << __FILE__<<" "<<__LINE__<<gradientOrigEBE.Norm() <<std::endl ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				//std::cout << std::setw(width)<< gradientScaledSave(count) << "   " ;
				if (abs(gradientOrigEBE(count,0))>tol)
					std::cout << " dof nr. "<< count << "   " << gradientOrigEBE(count,0)<<std::endl;

			}
			if(gradientOrigEBE.Norm() >1e-4)
		    {
		    	std::cout <<  "gradientOrig norm    "<< gradientOrigEBE.Norm() <<std::endl ;
		    	std::cout<<__FILE__<<" "<<__LINE__<<" error different results."<<std::endl;
		    	return 1;
		    }

			std::cout << std::endl;
			//gradientOrigNBN.Info();

*/
			//with scaling
			//gradientScaled = scaleFactorsInv.asDiagonal()*gradientOrig.mEigenMatrix;
			//no scaling
			gradientScaled = gradientOrig.mEigenMatrix;
			gradientNew =  gradientOrig.mEigenMatrix;
			//std::cout<<__FILE__<<" "<<__LINE__<<"scaled grad: \n"<<gradientScaled<<std::endl;

			int precision = 3;
			int width = 10;
			std::cout.precision(precision);
		    if (mVerboseLevel>3)
		    {
				std::cout << std::setw(width)<< "gradient scaled " ;
				for (int count=0; count<GetNumParameters(); count++)
				{
					std::cout << std::setw(width)<< gradientScaled(count) << "   " ;
				}
				std::cout << std::endl;

				std::cout << "gradient new ";
				for (int count=0; count<GetNumParameters(); count++)
				{
					std::cout << std::setw(width)<< gradientNew(count) << "   " ;
				}
				std::cout << std::endl;
		    }

			normGradient = gradientScaled.norm();
			//normGradient = gradientNew.norm();
			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" normGradient "<<normGradient << " accuracy " <<mAccuracyGradientScaled<<std::endl;


			if (normGradient<mAccuracyGradientScaled)
			{
				converged = true;
				returnValue = NORMGRADIENT;
				break;
			}
			alphaNumerator = gradientNew.dot(gradientScaled);
			//update for start solution
			searchDirectionOrig = gradientScaled;
			//initialize searchDirectionScaled with gradientScaled,
			//needed as input for CalculateScaledSearchDirection
			searchDirectionScaled = gradientScaled;
			if (mVerboseLevel>5 && curCycle>0)
				std::cout<< "   Restart after " <<curCycle << " cycles " << std::endl;
			curCycle = 0;
		}

		if (mVerboseLevel>2)
			std::cout<<__FILE__<<" "<<__LINE__<<" calc search direction"<<std::endl;

		int precision = 6;
		int width = 8;
		std::cout.precision(precision);

//		gradientScaledSave=searchDirectionScaled;
		if (!matrixFreeMethod)
			CalculateScaledSearchDirection(searchDirectionScaled);
		//gradientScaledSave=searchDirectionScaled;
		//searchDirectionScaled = gradientScaled;
		if (matrixFreeMethod)
			CalculateScaledSearchDirectionNodeByNodeII(searchDirectionScaled);
			//CalculateScaledSearchDirectionNodeByNodeII(searchDirectionScaled,numNodes,globNodeIds,globArray,globConstraint);
/*
		std::cout << std::setw(width)<< "EBE searchDirectionScaled  " ;
		for (int count=0; count<GetNumParameters(); count++)
		{
			std::cout << std::setw(width)<< searchDirectionScaled(count) << "   " ;
		}
		std::cout << std::endl;
		std::cout << std::setw(width)<< "NBN gradientScaledSave     " ;
		for (int count=0; count<GetNumParameters(); count++)
		{
			std::cout << std::setw(width)<< gradientScaledSave(count) << "   " ;
		}
		std::cout << std::endl;
*/
/*
 	    gradientScaledSave-=searchDirectionScaled;

		for (int count=0; count<GetNumParameters(); count++)
		{
			//std::cout << std::setw(width)<< gradientScaledSave(count) << "   " ;
			if (gradientScaledSave(count)>tol)
			{
				std::cout << std::setw(width)<< "gradientScaledDiff     " ;
				std::cout << " dof nr. "<< count << "   " << gradientScaledSave(count)<<std::endl;
			}
		}
		std::cout << std::endl;



	    if(gradientScaledSave.norm() >1e-4)
	    {
	    	std::cout << std::setw(width)<< "gradientScaledDiff norm    "<< gradientScaledSave.norm() <<std::endl ;
	    	std::cout<<__FILE__<<" "<<__LINE__<<" error different results."<<std::endl;
	    	return 1;
	    }
	    */
	    alphaDenominator = searchDirectionOrig.dot(searchDirectionScaled);
		alpha = alphaNumerator/alphaDenominator;

		// store previous parameter
		prevParameters = mvParameters.mEigenMatrix;
		prevGradientOrig = gradientNew;
		//set new parameter
		mvParameters.mEigenMatrix=prevParameters+alpha*searchDirectionOrig;
		gradientNew = prevGradientOrig - alpha*searchDirectionScaled;
		//gradientScaled = scaleFactorsInv.asDiagonal()*gradientNew;
		//no preconditioning
		gradientScaled = gradientNew;
		betaNumerator = gradientNew.dot(gradientScaled);
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;

		//normGradient = gradientNew.norm()*(double)GetNumParameters();
		normGradient = gradientNew.norm();
		if (normGradient<mAccuracyGradientScaled)
		{
			converged = true;
			returnValue = NORMGRADIENT;
			break;
		}

		if (beta<0)
		{
			std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}

		searchDirectionOrig *=beta;
		searchDirectionOrig +=gradientScaled;

		//update searchDirectionScaledfor input in function
		searchDirectionScaled=searchDirectionOrig;
		if (mVerboseLevel>2)
		{
			int precision = 3;
			int width = 10;
			std::cout.precision(precision);
			std::cout << std::setw(width)<< "gradient scaled " ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				std::cout << std::setw(width)<< gradientScaled(count) << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "gradient orig " ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				std::cout << std::setw(width)<< gradientNew(count,0) << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "searchDirectionScaled " ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				std::cout << std::setw(width)<< searchDirectionScaled(count) << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "searchDirectionOrig " ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				std::cout << std::setw(width)<< searchDirectionOrig(count) << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "displacements " ;
			for (int count=0; count<GetNumParameters(); count++)
			{
				std::cout << std::setw(width)<<mvParameters.mEigenMatrix(count,0) << "   " ;
			}
			std::cout << std::endl;

			std::cout << std::setw(width)<< "alpha "<<alpha<< "beta "<<beta << std::endl;
		}


		if (mVerboseLevel>1 && curIteration%mShowSteps==0)
//			std::cout<< "Iteration " << curIteration <<" with norm grad " << gradientNew.norm()/sqrt(GetNumParameters()) << std::endl;
			std::cout<< "Iteration " << curIteration <<" with norm grad " << gradientNew.norm() << std::endl;

		//increase iteration and curCycle
		curCycle++;
		curIteration++;

		if (curIteration>mMaxIterations)
		{
			converged = true;
			returnValue = MAXITERATIONS;
			break;
		}

		//set new parameters in search direction
		SetParameters(mvParameters);

		numFunctionCalls++;
		if (numFunctionCalls>mMaxFunctionCalls)
		{
			converged = true;
			returnValue = MAXFUNCTIONCALLS;
			break;
		}
	}
	isBuild = true;
	if (mVerboseLevel>0)
	{
		std::cout<< " "  << std::endl;
		std::cout<< "Number of Function Calls......... " << numFunctionCalls << std::endl;
		std::cout<< "Number of Gradient Calls......... " << numGradientCalls << std::endl;
		std::cout<< "Number of Hessian Calls.......... " << numHessianCalls << std::endl;
		std::cout<< "Number of Iterations............. " << curIteration << std::endl;
		std::cout<< "Active convergence criterion..... " ;
		switch (returnValue)
		{
			case MAXFUNCTIONCALLS:
				std::cout<< "Maximum number of function calls reached." << std::endl;
				break;
			case MAXGRADIENTCALLS:
				std::cout<< "Maximum number of gradient calls reached." << std::endl;
				break;
			case MAXHESSIANCALLS:
				std::cout<< "Maximum number of hessian calls reached." << std::endl;
				break;
			case MAXITERATIONS:
				std::cout<< "Maximum number of iterations reached." << std::endl;
				break;
			case NORMGRADIENT:
				std::cout<< "Norm of preconditioned gradient smaller than prescribed value." << std::endl;
				break;
			case MINOBJECTIVE:
				std::cout<< "Objective smaller than prescribed value." << std::endl;
				break;
			case REACHINGMACHINEPRECISION:
				std::cout<< "The machine precision is reached within the initial phase of the linesearch." << std::endl;
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
		}
		std::cout << std::endl;
		int precision = 6;
		int width = 10;
		std::cout.precision(precision);
		std::cout << std::setw(width)<< "displacements " ;
		for (int count=0; count<GetNumParameters(); count++)
		{
			std::cout << std::setw(width)<<mvParameters.mEigenMatrix(count,0) << "   " ;
		}
		std::cout << std::endl;

	}
	//return results to grid structure
#ifdef SHOW_TIME
    endOpt=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] " << difftime(endOpt,startOpt)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
	return returnValue;
}

void NuTo::ConjugateGradientGrid::CalcScalingFactors(int& numHessianCalls,Eigen::VectorXd& scaleFactorsInv)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    //diagonal scaling with scaling factor
	numHessianCalls++;
    double scalefactor=.5;
    for (int count=0; count<GetNumParameters(); ++count)
        scaleFactorsInv(count) =scalefactor;

    /*
    double scalefactor=1;
    for (int count=0; count<GetNumParameters(); ++count)
    {
        if  (hessianOrig(count,count)>1)
        {
            scaleFactorsInv(count) =scalefactor/hessianOrig(count,count);
        	std::cout <<"   "<<scaleFactorsInv(count);
       }
        else
        {
            std::cout<<__FILE__<<" "<<__LINE__<<"hessian value smaller then one."<<std::endl;
        	scaleFactorsInv(count) = 1.;
        }
    }
    */
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalcScalingFactors] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
void NuTo::ConjugateGradientGrid::Hessian(NuTo::FullMatrix<double>& rHessian)const
{
	if (mUseDiagHessian)
		HessianDiag(rHessian);
	else
		throw OptimizeException("[ConjugateGradientGrid::Hessian] Only diagonal hessian does exist.");
}

void NuTo::ConjugateGradientGrid::HessianDiag(NuTo::FullMatrix<double>& rHessian)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
#ifdef ENABLE_MECHANICS
	std::cout<<__FILE__<<" "<<__LINE__<<" in Routine ConjugateGradientGrid::HessianDiag"<<std::endl;
    int numElems=mpGrid->GetNumElements();
    NuTo::FullMatrix<int> *voxelLoc;
    voxelLoc=mpGrid->GetVoxelNumAndLocMatrix();

    int thisvoxelLocation[4]={0};
    int corners[8]={0};
	//array with global dof number of each dof of this element
    int numDofs=3; //for dofs array needed

   	//move all this in element loop
	Voxel8N* thisElement;
    //thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(0));
	int NumParameters =mpGrid->GetNumDofs();

	// initialize Hessian
	rHessian.mEigenMatrix.setZero(NumParameters,NumParameters);
	std::cout<<__FILE__<<" "<<__LINE__<<"numPara "<<NumParameters<<std::endl;
	for (int elementNumber=0;elementNumber<numElems;elementNumber++)
    {
    	//get pointer to each element
		thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementNumber));

		int dofsElem=thisElement->GetNumDofs();
		if (dofsElem!=24)
			throw OptimizeException("[ConjugateGradientGrid::HessianDiag] Number of Dofs is not 24.");
		int dofs[24]={0};
        //get grid location and number of corresponding voxel
        for (int count = 0;count<4;++count)
        	thisvoxelLocation[count]= voxelLoc->GetValue(elementNumber,count);

        //get grid corners of the voxel
        mpGrid->GetCornersOfVoxel(elementNumber, thisvoxelLocation, corners);

        //get the number of the material
        int numStiff = thisElement->GetNumLocalStiffnessMatrix();
        //get the local stiffness matrix
        NuTo::FullMatrix<double> *matrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);

 		//loop over all nodes of one element
        for (int node=0;node<thisElement->GetNumNodes();++node)
        {
			//get pointer to this gridNum node
//			NodeBase* thisNode =mpGrid->NodeGetNodePtrFromGridNum(corners[node]);
        	if (mVerboseLevel>3)
        		std::cout<<__FILE__<<" "<<__LINE__<<" node num "<<node<<" corner "<<corners[node]<<std::endl;
			NodeGrid3D* thisNode =mpGrid->NodeGetNodePtrFromGridNum(corners[node]);
			//which DOFs belonging to this node of this element
			for (int disp = 0;disp<numDofs;++disp)
			{	//save global dof number in local ordering
				dofs[node*numDofs+disp]=(thisNode->GetGlobalDofs())[disp];
				// set diagonal hessian

				//for case when take only active dofs
				//if (dofs[node*numDofs+disp]<mpGrid->NodeGetNumberActiveDofs())
				rHessian.SetValue(dofs[node*numDofs+disp],dofs[node*numDofs+disp],matrix->GetValue(node*numDofs+disp,node*numDofs+disp));
			}
        }
        // Invert Hessian to get Preconditionner
        // but not yet?
        //rHessian.ElementwiseInverse();
    }
//#else
	//throw OptimizeException ( "[ConjugateGradientGrid::HessianDiag] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::HessianDiag] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... calculate start gradient in element-by-element way
void NuTo::ConjugateGradientGrid::CalculateStartGradient(NuTo::FullMatrix<double> &gradientOrig)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif //SHOW_TIME

#ifdef ENABLE_MECHANICS
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradient"<<std::endl;
	int numElems=mpGrid->GetNumElements();

	//array with global dof number of each dof of this element
    int dofs[24]={0};
    //! @TODO replace static_cast
    int numDofs=3; //for dofs array needed
      //! @TODO replace static_cast
    //std::map<int> globtoloc;

   	//move all this in element loop
	Voxel8N* thisElement;
    thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(0));
	int dofsElem=thisElement->GetNumDofs();
	// return vector with all dofs of one element
	FullMatrix<double> locReturn(dofsElem,1);
	// array of all three nodal displacements for all element nodes
	FullMatrix<double> displacements(dofsElem,1);
	// return vector of active dofs
	std::vector<double> activeReturn;
	// global external force vector (active dofs)
	FullMatrix<double>  force(mpGrid->GetNumDofs(),1);

	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;elementNumber++)
    {
    	//get pointer to each element
		thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementNumber));

        //get the number of the material
        int numStiff = thisElement->GetNumLocalStiffnessMatrix();
        //get the local stiffness matrix
        NuTo::FullMatrix<double> *matrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);

 		//loop over all nodes of one element
		FullMatrix<double> locDispValues(3,1);
		int * nodeIds=thisElement->GetNodeIds();
        for (int node=0;node<thisElement->GetNumNodes();++node)
        {
			//get displacements of one node
			mpGrid->NodeGetDisplacements(nodeIds[node],locDispValues);
			/* slower than variant above:
			double locDispValues[3]={0};
			NodeGrid3D* thisNode =mpGrid->NodeGetNodePtrFromGridNum(corners[node]);
			thisNode->GetDisplacements3D(locDispValues);
			 */
			//which DOFs belonging to this node of this element
			for (int disp = 0;disp<numDofs;++disp)
			{
				//save global dof number in local ordering
				//no longer needed with new constraint saving
				//dofs[node*numDofs+disp]=(thisNode->GetGlobalDofs())[disp];
				dofs[node*numDofs+disp]=3*nodeIds[node]+disp;
				//save diplacements of all dofs for one element
				displacements(node*numDofs+disp,0)=locDispValues(disp,0);
			}
        }
        //calculate local return vector with all dofs: r=Ku
        locReturn = matrix->operator *(displacements);

        //update gradient vector for undependant dofs
        for (int count =0; count <dofsElem;++count)
        {
			//when global dof is active
        	//when dof is not constraint, then ...
        	if (mpGrid->NodeGetConstraintSwitch(dofs[count]))
        		//subtract (r=f-Ku) locReturn for active dofs
        		//in dofs[count] is the active dof number of each element dof
        		gradientOrig(dofs[count],0) -= locReturn(count,0);
       }
         //get global external force vector
        //! @TODO add load vector
        // ubpdate this routine
        //mpGrid->BuildGlobalExternalLoadVector(force);
        //add global external force vector
        //gradientOrig+=force;
    }
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradient] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradient] " << difftime(end,start)/CLOCKS_PER_SEC << "sec  " << std::endl;
#endif
}

//! @brief ... calculate start gradient in node-by-node way
void NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNode(NuTo::FullMatrix<double> &gradientOrig)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradientNodeByNode"<<std::endl;
	int numNodes=mpGrid->GetNumNodes();
    int* elementIds;
	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed

	// array of all three nodal displacements
	FullMatrix<double> displacements(24,1);
	// global external force vector (active dofs)
	FullMatrix<double>  force(mpGrid->GetNumDofs(),1);

	FullMatrix<double> locDispValues(3,1);
	//local field with edge location of the node in the eight elements
	int elemEdge[8]={6,7,4,5,2,3,0,1};

	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		//when dof is active
		if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs)||mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+1)||mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+2))
		{
			// return vector with all dofs of one node
			FullMatrix<double> locReturn(numDofs,1);
			//get pointer to each node
			NodeGrid3D* thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
			//mpGrid->NodeGetDisplacements(nodeNumber,locDispValues);
			//get belonging elements
			elementIds=thisNode->GetElementIds();
	//		int numElems=thisNode->GetNumElems();
			//get stiffness
			for (int element=0;element<8;++element)
			{
				NuTo::FullMatrix<double> matrix(3,24);
				Voxel8N* thisElement;
				//element exists
				if (elementIds[element]>=0)
				{
					thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementIds[element]));

					//get the number of the material
					int numStiff = thisElement->GetNumLocalStiffnessMatrix();
					//get the local stiffness matrix

					//Get 3X24 Block with location of node * numDofs (3)
					NuTo::FullMatrix<double> *elemMatrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);
					matrix= (elemMatrix->GetBlock(elemEdge[element]*numDofs,0,3,24));
					//get displacement vector for this element
					//loop over all nodes of one element
					FullMatrix<double> locDispValues(3,1);
					int * nodeIds=thisElement->GetNodeIds();
					for (int node=0;node<thisElement->GetNumNodes();++node)
					{
						mpGrid->NodeGetDisplacements(nodeIds[node],locDispValues);
						displacements(node*3,0)=locDispValues(0,0);
						displacements(node*3+1,0)=locDispValues(1,0);
						displacements(node*3+2,0)=locDispValues(2,0);
					}
					//calculate nodal return vector and sum up
					locReturn.operator += (matrix.operator *(displacements));

				}
			}
			//update gradient vector for undependant dofs
			for (int count =0; count <numDofs;++count)
			{
				//when global dof is active
				//when dof is not constraint, then ...
				if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+count))
					//subtract (r=f-Ku) locReturn for active dofs
					//in dofs[count] is the active dof number of each element dof
					gradientOrig(nodeNumber*numDofs+count,0) -= locReturn(count,0);
			}

		}
		//get global external force vector
		//! @TODO add load vector
		// ubdate this routine
		//mpGrid->BuildGlobalExternalLoadVector(force);
		//add global external force vector
		//gradientOrig+=force;
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNode] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradientNodeByNode] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}

//! @brief ... calculate start gradient in node-by-node way
//! @brief ... variante II: with 3x3 matrix at node
void NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNodeII(NuTo::FullMatrix<double> &gradientOrig)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradientNodeByNodeII"<<std::endl;
	FullMatrix<double> gradientNew(GetNumParameters(),1);
	int numNodes=mpGrid->GetNumNodes();
	//nodeIds here ids of all neighbor nodes
    int* nodeIds;
	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed

	// part of matrix one dimensional 9 fields
	double *array=new double [9*27];

	// global external force vector (active dofs)
//	FullMatrix<double>  force(mpGrid->GetNumDofs(),1);

	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		// return vector with all dofs of one node
		FullMatrix<double> locReturn(numDofs,1);
		//get pointer to each node
		NodeGrid3D* thisNode;
		thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
		//get belonging nodes
		nodeIds=thisNode->GetNodeIds();
//				std::cout<<" node "<<nodeNumber <<""<<std::endl;
		for (int node=0;node<27;++node)
		{
			//NuTo::FullMatrix<double>* matrix;
			//get pointer to array of part coefficient matrix for all neighbor nodes
			array= thisNode->GetPartCoefficient0(0);
//			matrix= thisNode->GetPartCoefficientMatrix0(node);
			//node exists
	/*
			if (nodeNumber==62)
			{
				std::cout<<" node "<<nodeNumber <<"ids :"<<nodeIds[0]<<"  "<<nodeIds[1]<< std::endl;
			}
		*/
			if (nodeIds[node]>=0)
			{
				//FullMatrix<double> locDispValues(3,1);
				//mpGrid->NodeGetDisplacements(nodeIds[node],locDispValues);
				//calculate nodal return vector and sum up
				locReturn(0,0)+=array[node*9+0]*gradientOrig(nodeIds[node]*numDofs,0) + array[node*9+1]*gradientOrig(nodeIds[node]*numDofs+1,0) + array[node*9+2]*gradientOrig(nodeIds[node]*numDofs+2,0);
				locReturn(1,0)+=array[node*9+3]*gradientOrig(nodeIds[node]*numDofs,0) + array[node*9+4]*gradientOrig(nodeIds[node]*numDofs+1,0) + array[node*9+5]*gradientOrig(nodeIds[node]*numDofs+2,0);
				locReturn(2,0)+=array[node*9+6]*gradientOrig(nodeIds[node]*numDofs,0) + array[node*9+7]*gradientOrig(nodeIds[node]*numDofs+1,0) + array[node*9+8]*gradientOrig(nodeIds[node]*numDofs+2,0);
	//				locReturn.operator += (matrix->operator *(locDispValues));
			/*
				if (nodeNumber==62)
				{
					std::cout<<" node "<<nodeNumber <<" dof 0 locReturn :"<<locReturn(0,0)<<" at loc node "<<nodeIds[node]<< std::endl;
					std::cout<<" gradOrig :"<<gradientOrig(nodeIds[node]*numDofs,0)<<" , "<< gradientOrig(nodeIds[node]*numDofs+1,0)<<" , "<<gradientOrig(nodeIds[node]*numDofs+2,0)<<std::endl;
				}
				*/
			}

		}
		//when global dof is active
		//when dof is not constraint, then ...
		//subtract (r=f-Ku) locReturn for active dofs
		if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+0))
			gradientNew(nodeNumber*numDofs+0,0) = -locReturn(0,0);
		if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+1))
			gradientNew(nodeNumber*numDofs+1,0) = -locReturn(1,0);
		if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+2))
			gradientNew(nodeNumber*numDofs+2,0) = -locReturn(2,0);
	}
	gradientOrig=gradientNew;
		//get global external force vector
		//! @TODO add load vector
		// ubdate this routine
		//mpGrid->BuildGlobalExternalLoadVector(force);
		//add global external force vector
		//gradientOrig+=force;

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradientNodeByNodeII] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}

//! @brief ... calculate matix-vector product in element-by-element way
//! @brief ... multiply each element matrix with search direction
void NuTo::ConjugateGradientGrid::CalculateMatrixVectorEBE(bool startSolution, NuTo::FullMatrix<double> &returnVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
#ifdef ENABLE_MECHANICS
	//check if mvParamters is initialized
	if (!&mvParameters)
 		throw OptimizeException("[ConjugateGradientGrid::GetStartGradient] mvParameters not initialized.");

	//set search direction equal mvParameters
	FullMatrix<double>origSearchDirection=mvParameters;

	//new: multiply each time hole element stiffness matrix, reduce instead return vector
	Voxel8N* thisElement=0;
   	thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(0));
	int dofsElem=thisElement->GetNumDofs();

	int numElems=mpGrid->GetNumElements();
    NuTo::FullMatrix<int> *voxelLoc;
    voxelLoc=mpGrid->GetVoxelNumAndLocMatrix();
    int thisvoxelLocation[4]={0};
    int corners[8]={0};
	//array with global dof number of each dof of this element
    int dofs[24]={0};
    int numDofs=3; //for dofs array needed

 	FullMatrix<double> locReturn(dofsElem,1);
	std::vector<double> activeReturn;
	FullMatrix<double> origLocSearch(dofsElem,1);

	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;elementNumber++)
    {
    	//get pointer to each element
		thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementNumber));

        //get grid location and number of corresponding voxel
        for (int count = 0;count<4;++count)
        	thisvoxelLocation[count]= voxelLoc->GetValue(elementNumber,count);

        //get grid corners of the voxel
        mpGrid->GetCornersOfVoxel(elementNumber, thisvoxelLocation, corners);

        //get the number of the material
        int numStiff = thisElement->GetNumLocalStiffnessMatrix();
        //get the local stiffness matrix
        NuTo::FullMatrix<double> *matrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);

 		//loop over all nodes of one element
        for (int node=0;node<thisElement->GetNumNodes();++node)
        {
			//get pointer to this gridNum node
//			NodeBase* thisNode =mpGrid->NodeGetNodePtrFromGridNum(corners[node]);
			NodeGrid3D* thisNode =mpGrid->NodeGetNodePtrFromGridNum(corners[node]);

			//which DOFs belonging to this node of this element
			for (int disp = 0;disp<numDofs;++disp)
				//save global dof number in local ordering
				dofs[node*numDofs+disp]=(thisNode->GetGlobalDofs())[disp];
        }
        //loop over all dofs of one element
        for (int count=0;count<dofsElem;++count)
        {
			//when global dof is active
        	if (dofs[count]<mpGrid->GetNumActiveDofs())
				//then write value of this dof to local vector
        		origLocSearch(count,0)=origSearchDirection(dofs[count],0);
        	//else global dof is dependent
        	else
        		//fill value of local (orig.) search direction vector with zero
        		origLocSearch(count,0)=0;
        }

        //calculate local return vector with all dofs: r=Ku
        locReturn = matrix->operator *(origLocSearch);
        //add result to return vector for active dofs only
        //important: do not take dependent dofs, because they are also calculated (not zero) in order to simplify the routine
        for (int count =0; count <numDofs;++count)
        {
			//when global dof is active
        	if (dofs[count]<mpGrid->GetNumActiveDofs())
        		//add locReturn for active dofs
        		returnVector(dofs[count],0) += locReturn(count,0);
        }
    }
#else
		throw OptimizeException ( "[ConjugateGradientGrid::CalculateMatrixVectorEBE] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateMatrixVectorEBE] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... calculate scaled search direction multiplied with stiffness matrix in element-by-element way for each step
void NuTo::ConjugateGradientGrid::CalculateScaledSearchDirection(Eigen::VectorXd& searchDirectionScaled)
{
/*
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
*/
#ifdef ENABLE_MECHANICS
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" "<<"in CalculateScaledSearchDirection"<<std::endl;
    int numElems=mpGrid->GetNumElements();
	//array with global dof number of each dof of this element
	int dofs[24]={0};
	int numDofs=3; //for dofs array neede

	//move all this in element loop
	Voxel8N* thisElement;
	thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(0));
	int dofsElem=thisElement->GetNumDofs();
	// return vector of all dofs
	FullMatrix<double> allReturn(mpGrid->GetNumDofs(),1);
	// return vector with all dofs of one element
	FullMatrix<double> locReturn(dofsElem,1);

	//local search direction
	//replaced std::vector by FullMatrix for multiplication, but tests shows same results also for vector
	//std::vector<double> searchDirectionLocal(24);
	FullMatrix<double> searchDirectionLocal(dofsElem,1);

	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;elementNumber++)
	{
		//get pointer to each element
		thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementNumber));

		//get the number of the material
		int numStiff = thisElement->GetNumLocalStiffnessMatrix();
		//get the local stiffness matrix
		NuTo::FullMatrix<double> *matrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);
		int * nodeIds=thisElement->GetNodeIds();

		//loop over all nodes of one element
		for (int node=0;node<8;++node)
		{

			//save global dof number in local ordering
			dofs[node*numDofs]=3*nodeIds[node];
			dofs[node*numDofs+1]=3*nodeIds[node]+1;
			dofs[node*numDofs+2]=3*nodeIds[node]+2;
			//save local search direction
			searchDirectionLocal(node*numDofs,0) = searchDirectionScaled(dofs[node*numDofs]);
			searchDirectionLocal(node*numDofs+1,0) = searchDirectionScaled(dofs[node*numDofs+1]);
			searchDirectionLocal(node*numDofs+2,0) = searchDirectionScaled(dofs[node*numDofs+2]);
		}
		//calculate local return vector
		locReturn = matrix->operator *(searchDirectionLocal);
 		//reduce return vector for element with dependant dofs

		for (int count =0; count <dofsElem;++count)
		{
			//when global dof is active
			if (mpGrid->NodeGetConstraintSwitch(dofs[count]))
				//subtract (r=f-Ku) locReturn for active dofs
				allReturn(dofs[count],0) += locReturn(count,0);
	   }
 	}
	searchDirectionScaled = allReturn.mEigenMatrix;
#else
 	throw OptimizeException ( "[ConjugateGradientGrid::CalculateScaledSearchDirection] Modul Mechanics is not loaded." );

#endif // ENABLE_MECHANICS
/*
#ifdef SHOW_TIME
    end=clock();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
     diffn=diff(startn,endn);
     if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirection            ] "<< diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec"<<std::endl;
        //std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirection] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
*/
}

//! @brief ... calculate search direction in node-by-node way
void NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNode(Eigen::VectorXd& searchDirectionScaled)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateScaledSearchDirectionNodeByNode"<<std::endl;
	int numNodes=mpGrid->GetNumNodes();
    int* elementIds;
	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed
	// return vector of all dofs
	FullMatrix<double> allReturn(mpGrid->GetNumDofs(),1);
	//local field with edge location of the node in the eight elements
	int elemEdge[8]={6,7,4,5,2,3,0,1};

	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		// array of all three nodal values of searchDirection
		FullMatrix<double> searchDirectionLocal(24,1);
		// return vector with all dofs of one node
		FullMatrix<double> locReturn(numDofs,1);
		//get pointer to each node
		NodeGrid3D* thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
		//get belonging elements
		elementIds=thisNode->GetElementIds();
		//get stiffness
		for (int element=0;element<8;++element)
		{
			//element exists
			if (elementIds[element]>=0)
			{
				Voxel8N* thisElement;
				thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(elementIds[element]));

				//get the number of the material
				int numStiff = thisElement->GetNumLocalStiffnessMatrix();
				//get the local stiffness matrix

				//Get 3X24 Block with location of node * numDofs (3)
				NuTo::FullMatrix<double> *elemMatrix = mpGrid->GetLocalCoefficientMatrix0(numStiff);
				NuTo::FullMatrix<double> matrix(3,24);
				matrix= (elemMatrix->GetBlock(elemEdge[element]*numDofs,0,3,24));
				//loop over all nodes of one element
				int * nodeIds=thisElement->GetNodeIds();
				for (int node=0;node<thisElement->GetNumNodes();++node)
				{
					//which DOFs belonging to this node of this element
					for (int disp = 0;disp<numDofs;++disp)
					{
						// if this dof is not constraint
							if (mpGrid->NodeGetConstraintSwitch(3*nodeIds[node]+disp))
						{
//			std::cout<<__LINE__<<"| node*numDofs+disp "<< node*numDofs+disp<<" 3*nodeIds[node]+disp "<<3*nodeIds[node]+disp<<std::endl;
							searchDirectionLocal(node*numDofs+disp,0)= searchDirectionScaled(3*nodeIds[node]+disp);
						}
						else
							searchDirectionLocal(node*numDofs+disp,0) = 0;
					}

				}
				//calculate nodal return vector and sum up
				locReturn.operator += (matrix.operator *(searchDirectionLocal));
			}
		}
		//update gradient vector for undependant dofs
		for (int count =0; count <numDofs;++count)
		{
			//when global dof is active
			//when dof is not constraint, then ...
			if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+count))
			{	//subtract (r=f-Ku) locReturn for active dofs
				//in dofs[count] is the active dof number of each element dof
//				allReturn(nodeNumber*numDofs+count,0) += locReturn(count,0);
				allReturn(nodeNumber*numDofs+count,0) = locReturn(count,0);
			}
		}

	}
	searchDirectionScaled = allReturn.mEigenMatrix;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNode] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradientNodeByNode] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}
/*
//! @brief ... calculate search direction in node-by-node way
void NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII(Eigen::VectorXd& searchDirectionScaled)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateScaledSearchDirectionNodeByNodeII"<<std::endl;
	int numNodes=mpGrid->GetNumNodes();
	//nodeIds here ids of all neighbor nodes
    int* nodeIds;
 	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed
	// return vector of all dofs
	FullMatrix<double> allReturn(mpGrid->GetNumDofs(),1);

	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		// return vector with all dofs of one node
		FullMatrix<double> locReturn(numDofs,1);
		//get pointer to each node
		NodeGrid3D* thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
		//get belonging nodes
		nodeIds=thisNode->GetNodeIds();
		for (int node=0;node<27;++node)
		{
			//node exists
			if (nodeIds[node]>=0)
			{
				// array of all three nodal values of searchDirection
				FullMatrix<double> searchDirectionLocal(numDofs,1);
				NuTo::FullMatrix<double>* matrix;
				matrix= thisNode->GetPartCoefficientMatrix0(node);
				//which DOFs belonging to this node of this element
				for (int count = 0;count<numDofs;++count)
				{
					searchDirectionLocal(count,0)= searchDirectionScaled(numDofs*nodeIds[node]+count);
				}
				//calculate nodal return vector and sum up
				locReturn.operator += (matrix->operator *(searchDirectionLocal));
			}
		}
		//update gradient vector for undependant dofs
		for (int count =0; count <numDofs;++count)
		{
			//when global dof is active
			//when dof is not constraint, then ...
			if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+count))
			{
				allReturn(nodeNumber*numDofs+count,0) = locReturn(count,0);
			}
		}


	}
	searchDirectionScaled = allReturn.mEigenMatrix;
#ifdef SHOW_TIME
    end=clock();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
    diffn=diff(startn,endn);
    if (mShowTime)
       std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] "<< diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec"<<std::endl;
 //      std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}
*/

//! @brief ... calculate search direction in node-by-node way
void NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII(Eigen::VectorXd& searchDirectionScaled)
{
#ifdef ENABLE_MECHANICS
/*
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
*/
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateScaledSearchDirectionNodeByNodeII"<<std::endl;
	int numNodes=mpGrid->GetNumNodes();
	//nodeIds here ids of all neighbor nodes
    int* nodeIds;
 	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed
	// return vector of all dofs
	FullMatrix<double> allReturn(mpGrid->GetNumDofs(),1);
	// part of matrix one dimensional 9 fields
	double *array=new double [9*27];

	bool *constraint = new bool[numNodes*numDofs];
	constraint = mpGrid->GetConstraintSwitch();
	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		//get pointer to each node
		NodeGrid3D* thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
		//get belonging nodes
		nodeIds=thisNode->GetNodeIds();
		//get pointer to array of part coefficient matrix for all neighbor nodes
		array= thisNode->GetPartCoefficient0(0);
		for (int node=0;node<27;++node)
		{
			//node exists
			if (nodeIds[node]>=0)
			{
				if (constraint[nodeNumber*numDofs])
					allReturn(nodeNumber*numDofs,0) += array[node*9+0]*searchDirectionScaled(numDofs*nodeIds[node]+0) + array[node*9+1]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+2]*searchDirectionScaled(numDofs*nodeIds[node]+2);
				if (constraint[nodeNumber*numDofs+1])
					allReturn(nodeNumber*numDofs+1,0) += array[node*9+3]*searchDirectionScaled(numDofs*nodeIds[node]+0) + array[node*9+4]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+5]*searchDirectionScaled(numDofs*nodeIds[node]+2);
				if (constraint[nodeNumber*numDofs+2])
					allReturn(nodeNumber*numDofs+2,0) += array[node*9+6]*searchDirectionScaled(numDofs*nodeIds[node]+0) +  array[node*9+7]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+8]*searchDirectionScaled(numDofs*nodeIds[node]+2);
			}
		}
	}
	searchDirectionScaled = allReturn.mEigenMatrix;
/*
#ifdef SHOW_TIME
    end=clock();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
    diffn=diff(startn,endn);
    if (mShowTime)
       std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] "<< diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec"<<std::endl;
 //      std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
*/
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}
void NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII(Eigen::VectorXd& searchDirectionScaled,int numNodes, int* globNodeIds, double* globArray, bool* constraint)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateScaledSearchDirectionNodeByNodeII"<<std::endl;
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateScaledSearchDirectionNodeByNodeII"<<std::endl;
	//nodeIds here ids of all neighbor nodes
    int* nodeIds;
 	//array with global dof number of each dof of this node
    int numDofs=3; //for dofs array needed
	// return vector of all dofs
	FullMatrix<double> allReturn(numDofs*numNodes,1);
	// part of matrix one dimensional 9 fields
	double *array=new double [9*27];

	//loop over all nodes
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		//get belonging nodes
		nodeIds=&(globNodeIds[nodeNumber*27]);
		std::cout<<" node ids of node "<<nodeNumber <<": "<<std::endl;
		//get pointer to array of part coefficient matrix for all neighbor nodes
		array= &(array[nodeNumber*27*9]);
		for (int node=0;node<27;++node)
		{
			std::cout<<nodeIds[node]<<" ," ;
			//node exists
			if (nodeIds[node]>=0)
			{
				if (constraint[nodeNumber*numDofs])
					allReturn(nodeNumber*numDofs,0) += array[node*9+0]*searchDirectionScaled(numDofs*nodeIds[node]+0) + array[node*9+1]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+2]*searchDirectionScaled(numDofs*nodeIds[node]+2);
				if (constraint[nodeNumber*numDofs+1])
					allReturn(nodeNumber*numDofs+1,0) += array[node*9+3]*searchDirectionScaled(numDofs*nodeIds[node]+0) + array[node*9+4]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+5]*searchDirectionScaled(numDofs*nodeIds[node]+2);
				if (constraint[nodeNumber*numDofs+2])
					allReturn(nodeNumber*numDofs+2,0) += array[node*9+6]*searchDirectionScaled(numDofs*nodeIds[node]+0) +  array[node*9+7]*searchDirectionScaled(numDofs*nodeIds[node]+1) + array[node*9+8]*searchDirectionScaled(numDofs*nodeIds[node]+2);
			}
		}
		std::cout<<std::endl;
	}
	searchDirectionScaled = allReturn.mEigenMatrix;
#ifdef SHOW_TIME
    end=clock();
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
    diffn=diff(startn,endn);
    if (mShowTime)
       std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] "<< diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec"<<std::endl;
 //      std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateScaledSearchDirectionNodeByNodeII] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
}



#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Save ( const std::string &filename, std::string rType)const
{
	try
	{
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
		std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
		std::string tmpStr ( GetTypeId() );
		std::string baseClassStr = tmpStr.substr ( 4,100 );
		if (rType=="BINARY")
		{
			boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
			oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_oarchive ota ( ofs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else
		{
			throw MathException ( "[FullMatrix::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception e )
	{
		std::string s ( std::string ( "[FullMatrix::Save]File save exception in boost - " ) +std::string ( e.what() ) );
		std::cout << s << "\n";
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
	catch ( ... )
	{
		throw MathException ( "[Matrix::Save]Unhandled exception." );
	}
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Restore ( const std::string &filename,  std::string rType)
{
    try
    {
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        std::string tmpString;
		if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else
		{
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MathException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MathException ( e.what() );
    }
    catch ( ... )
    {
        throw MathException ( "[Matrix::Restore]Unhandled exception." );
    }
}
#endif // ENABLE_SERIALIZATION

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name ConjugateGradientGrid
std::string NuTo::ConjugateGradientGrid::GetTypeId()const
{
    return std::string("ConjugateGradientGrid");
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ConjugateGradientGrid::Info () const
{
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}

