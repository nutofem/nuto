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
    std::clock_t start,end;
    start=clock();
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

	optimization_return_attributes returnValue;

	FullMatrix<double> gradientOrig(GetNumParameters(),1);
	FullMatrix<double> hessianOrig(GetNumParameters(),1);
	Eigen::VectorXd prevParameters;
	Eigen::VectorXd gradientScaled;
	Eigen::VectorXd scaleFactorsInv(GetNumParameters());
	Eigen::VectorXd prevGradientOrig(GetNumParameters());
	Eigen::VectorXd gradientNew(GetNumParameters());
	Eigen::VectorXd searchDirectionScaled(GetNumParameters());
	Eigen::VectorXd searchDirectionOrig(GetNumParameters());

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< " Para "<< GetNumParameters() << std::endl;
	bool converged(false);
	double mAccuracyGradientScaled = mAccuracyGradient;
	SetMaxGradientCalls(GetNumParameters()*GetNumParameters());
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
		if (curCycle%GetNumParameters()==0)
		{
			  // initialize search direction with steepest descent
			// calculate Hessian for scaling
			CalcScalingFactors(numHessianCalls,scaleFactorsInv);
			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}
			//calculate gradient as a start solution
			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" calc start direction"<<std::endl;
			CalculateStartGradient(gradientOrig);
			//gradientOrig.Info();
			//with scaling
			gradientScaled = scaleFactorsInv.asDiagonal()*gradientOrig.mEigenMatrix;
			//no scaling
			//gradientScaled = gradientOrig.mEigenMatrix;
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
		curCycle = 0;
		}
		searchDirectionOrig = gradientScaled;
		//initialize searchDirectionScaled with gradientScaled,
		//needed as input for CalculateScaledSearchDirection
		searchDirectionScaled = gradientScaled;

		if (mVerboseLevel>2)
			std::cout<<__FILE__<<" "<<__LINE__<<" calc search direction"<<std::endl;
		CalculateScaledSearchDirection(searchDirectionScaled);
		alphaDenominator = searchDirectionOrig.dot(searchDirectionScaled);
		alpha = alphaNumerator/alphaDenominator;
		if (mVerboseLevel>5 && curCycle>0)
			std::cout<< "   Restart after " <<curCycle << " cycles, " << std::endl;

		// store previous parameter
		prevParameters = mvParameters.mEigenMatrix;
		prevGradientOrig = gradientNew;
		//set new parameter
		mvParameters.mEigenMatrix=prevParameters+alpha*searchDirectionOrig;
		gradientNew = prevGradientOrig - alpha*searchDirectionScaled;
		gradientScaled = scaleFactorsInv.asDiagonal()*gradientNew;
		//no preconditioning
		//gradientScaled = gradientNew;
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
		curCycle = 0;

		if (beta<0)
		{
			std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}

		searchDirectionOrig *=beta;
		searchDirectionOrig +=gradientScaled;
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
		std::cout<< "Norm of preconditioned gradient.. " << gradientScaled.norm() << std::endl;
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
		int precision = 3;
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
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
    double scalefactor=1;
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
#endif
#ifdef ENABLE_MECHANICS
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradient"<<std::endl;
    int numElems=mpGrid->GetNumElements();
    NuTo::FullMatrix<int> *voxelLoc;
    voxelLoc=mpGrid->GetVoxelNumAndLocMatrix();

    int thisvoxelLocation[4]={0};
    int corners[8]={0};
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
	FullMatrix<double>  force(mpGrid->GetNumActiveDofs(),1);

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

			double locDispValues[3]={0};
			//get displacements of one node
			thisNode->GetDisplacements3D(locDispValues);
			//which DOFs belonging to this node of this element
			for (int disp = 0;disp<numDofs;++disp)
			{
				//save global dof number in local ordering
				dofs[node*numDofs+disp]=(thisNode->GetGlobalDofs())[disp];
				//save diplacements of all dofs
				displacements(node*numDofs+disp,0)=locDispValues[disp];
			}
        }
        //calculate local return vector with all dofs: r=Ku
        locReturn = matrix->operator *(displacements);

        //update gradient vector for undependant dofs
        for (int count =0; count <dofsElem;++count)
        {
			//when global dof is active
        	if (dofs[count]<mpGrid->GetNumActiveDofs())
        		//subtract (r=f-Ku) locReturn for active dofs
        		gradientOrig(dofs[count],0) -= locReturn(count,0);
       }

        //get global external force vector
        mpGrid->BuildGlobalExternalLoadVector(force);
        //add global external force vector
        gradientOrig+=force;

    }
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradient] Modul Mechanics is not loaded." );
#endif // ENABLE_MECHANICS
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradient] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
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
				//! TODO save global dofs number of each element?!
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
#ifdef SHOW_TIME
    std::clock_t start,end,between;
    start=clock();
#endif
#ifdef ENABLE_MECHANICS
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" "<<"in CalculateScaledSearchDirection"<<std::endl;
   int numElems=mpGrid->GetNumElements();
	NuTo::FullMatrix<int> *voxelLoc;
	voxelLoc=mpGrid->GetVoxelNumAndLocMatrix();

	int thisvoxelLocation[4]={0};
	int corners[8]={0};
	//array with global dof number of each dof of this element
	int dofs[24]={0};
	int numDofs=3; //for dofs array needed
	//std::map<int> globtoloc;

	//move all this in element loop
	Voxel8N* thisElement;
	thisElement= static_cast<Voxel8N*>( mpGrid->ElementGetElementPtr(0));
	int dofsElem=thisElement->GetNumDofs();
	// return vector with all dofs of one element
	FullMatrix<double> locReturn(dofsElem,1);
	// return vector of all dofs
	FullMatrix<double> activeReturn(mpGrid->GetNumActiveDofs(),1);
	//local search direction
	std::vector<double> searchDirectionLocal(24);
#ifdef SHOW_TIME
	between=clock();
//	if (mShowTime)
//	    std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirection] " << difftime(between,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

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
			{
				//save global dof number in local ordering
				dofs[node*numDofs+disp]=(thisNode->GetGlobalDofs())[disp];
				//save local search direction
				if (dofs[node*numDofs+disp]<mpGrid->GetNumActiveDofs())
					  searchDirectionLocal[node*numDofs+disp] = searchDirectionScaled(dofs[node*numDofs+disp]);
				else
					searchDirectionLocal[node*numDofs+disp] = 0;
			}
		}
		//calculate local return vector
		locReturn = matrix->operator *(searchDirectionLocal);
		//reduce return vector for element with dependant dofs
		for (int count =0; count <dofsElem;++count)
		{
			//when global dof is active
			if (dofs[count]<mpGrid->GetNumActiveDofs())
				//subtract (r=f-Ku) locReturn for active dofs
				activeReturn(dofs[count],0) += locReturn(count,0);
	   }
	}
	searchDirectionScaled = activeReturn.mEigenMatrix;
#else
 	throw OptimizeException ( "[ConjugateGradientGrid::CalculateScaledSearchDirection] Modul Mechanics is not loaded." );

#endif // ENABLE_MECHANICS
#ifdef SHOW_TIME
    end=clock();
//    if (mShowTime)
  //      std::cout<<"[NuTo::ConjugateGradientGrid::CalculateScaledSearchDirection] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
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

