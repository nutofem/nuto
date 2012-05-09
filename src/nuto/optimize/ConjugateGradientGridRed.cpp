// $Id $

//! @author Andrea Keßler, ISM
//! @date December 2011
//! @brief ... conjugate gradient grid with reduced vectors only existing dofs


#include "nuto/optimize/ConjugateGradientGridRed.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
#define tol 1e-8

int NuTo::ConjugateGradientGridRed::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startOpt,endOpt;
    startOpt=clock();
#endif
	double alpha=0,
		   beta=0,
		   alphaNumerator=0,
		   alphaDenominator=0,
		   betaNumerator=0,
		   squaredNorm=0;

    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
		 curCycle(0);           //number of iterations without restart

	if (mVerboseLevel>5)
	{
		std::cout<<__FILE__<<" "<<__LINE__<<" mMatrixFreeMethod is";
		if (mMatrixFreeMethod)
		{
			std::cout<<" NBN \n";
		}
		else
			std::cout<<" EBE \n";
	}

	optimization_return_attributes returnValue;

	// u  = parameters
	// r  = gradient
	// p  = preconditioner
	// pr = preconditioned gradient = p r
	// d  = search direction
	// h  = scaled search direction = K d

#ifdef ENABLE_EIGEN
	std::cout<<"[NuTo::ConjugateGradientGridRed::Optimize] Eigen vectors are used.\n";
	myType u;
	u.MapAligned(mParameters.data(),mNumParameters);
	u.push_back(0.);
	myType f(mNumParameters);
	if(mForces.size()!=0)
		f.MapAligned(mForces.data(),mNumParameters);

	myType r(mNumParameters+1);
	myType pr;
	myType p(mNumParameters+1);
	myType h(mNumParameters+1);
	myType d(mNumParameters+1);
#else
	//std::vector
	std::cout<<"[NuTo::ConjugateGradientGridRed::Optimize] std::vector is used.\n";
	myType &u=mParameters;
	myType &f=mForces;
	myType r(mNumParameters);
	myType pr(mNumParameters);
	myType p(mNumParameters,1);
	myType h(mNumParameters);
	myType d(mNumParameters);
#endif
	int precision = 6;
	int width = 10;

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< " Para "<< mNumParameters << std::endl;
	bool converged(false);
	double mAccuracyGradientScaled = mAccuracyGradient;
	std::cout<<"[NuTo::ConjugateGradientGridRed::Optimize] gradient accuracy "<<mAccuracyGradientScaled <<std::endl;

	int localMaxGradientCalls=mNumParameters;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);

	// calculate objective
	numFunctionCalls++;
	if (numFunctionCalls>mMaxFunctionCalls)
	{
		converged = true;
		returnValue = MAXHESSIANCALLS;
	}

	//0 set if no precondition
//	std::fstream outputTime;
//	std::string filename = "timeOutput";
//  outputTime.open(filename,std::fstream::out|std::fstream::app);
//  outputTime<<" 0 ";
//	outputTime.close();
	// calculate Diag preonditioner;
	HessianDiag(p);

	//1 set if no scaling
//	std::fstream outputTime;
//	std::string filename = "timeOutput";
//    outputTime.open(filename,std::fstream::out|std::fstream::app);
//    outputTime<<" 1 ";
//	outputTime.close();
	CalcScalingFactors(numHessianCalls,p);

//	 print p
//	std::cout<<" p ";
//	for(int i=0;i<mNumParameters;i+=3)
//		std::cout<<p[i]<<" ";
//	std::cout<<"\n";


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

		// all  mNumParameters repeat start equations
		//if (curCycle%mNumParameters==0)
		// no upate with start equations
		if (curCycle==0)
		{
//			CalcScalingFactors(numHessianCalls,p);
			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}
			//calculate gradient as a start solution
			if (mVerboseLevel>3)
				std::cout<<__FILE__<<" "<<__LINE__<<" calc start direction"<<std::endl;
			if(!mMatrixFreeMethod)
				CalculateMatrixVectorProductEBE(u,r);
			else
				CalculateMatrixVectorProductNBN(u,r);

			//get global external force vector
			if (mVerboseLevel>3)
			{
				if(mForces.empty())
				std::cout<<" no external force.\n";
			else
				std::cout<<" with external force.\n";
			}

#ifdef ENABLE_EIGEN
			r*=-1;
			if(!mForces.empty())
				r+=f;
			//or r=KU-f test
//			r-=f;
#else
			for(size_t i=0;i<mNumParameters;++i)
			{
				if(r[i]!=0)
					r[i]*=-1;
				if(!mForces.empty())
					r[i]+=f[i];
			}
#endif
			// print r
//			std::cout<<" r ";
//			for(size_t i=0;i<mNumParameters;++i)
//				std::cout<<r[i]<<" ";
//			std::cout<<"\n";

#ifdef ENABLE_EIGEN
			pr = p*r;
			//upate for start solution
			d = pr;
#else
//			pr=r;
			for(size_t i=0;i<mNumParameters;++i)
			{
//				pr[i]=r[i];
 				pr[i]=p[i]*r[i];
				//upate for start solution
 				d[i]=pr[i];
			}
#endif

#ifdef ENABLE_EIGEN
			alphaNumerator = r.dot(pr);
#else
			alphaNumerator=0;
			for(size_t i=0;i<mNumParameters;++i)
 				alphaNumerator+=r[i]*pr[i];
#endif
			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" normGradient "<<alphaNumerator << " accuracy " <<mAccuracyGradientScaled<<std::endl;

			if (mVerboseLevel>5 && curCycle>0)
				std::cout<< "   Restart after " <<curCycle << " cycles " << std::endl;
			curCycle = 0;
		std::cout.precision(precision);
//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" a1 "<<alphaNumerator<<"\n";
		}


		//begin with cycles bigger 0

		if (mVerboseLevel>2)
			std::cout<<__FILE__<<" "<<__LINE__<<" calc search direction"<<std::endl;


		for(size_t i=0;i<mNumParameters;++i)
			h[i]=0;

		if(!mMatrixFreeMethod)
			CalculateMatrixVectorProductEBE(d,h);
		else
			CalculateMatrixVectorProductNBN(d,h);
#ifdef ENABLE_EIGEN
	    alphaDenominator = d.dot(h);
		alpha = alphaNumerator/alphaDenominator;
		//set new parameter
		u+= alpha*d;
		r-=alpha*h;
		pr = p*r;
		//no preconditioning
//		pr = r;
		betaNumerator = r.dot(pr);
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;
		d *=beta;
		d +=pr;
#else
	    alphaDenominator=0;
	    squaredNorm=0;
		for(size_t i=0;i<mNumParameters;++i)
			alphaDenominator+=d[i]*h[i];
		alpha = alphaNumerator/alphaDenominator;
		betaNumerator =0;
		for(size_t i=0;i<mNumParameters;++i)
		{
			u[i]+= alpha*d[i];
			r[i]-=alpha*h[i];
			pr[i]=p[i]*r[i];
			squaredNorm+=r[i]*r[i];
			betaNumerator+=r[i]*pr[i];
		}
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;
		for(size_t i=0;i<mNumParameters;++i)
		{
			d[i] *=beta;
			d[i] +=pr[i];
		}
#endif
//		std::cout.precision(precision);
//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> h "<<h[3]<<" alpha "<<alpha<<" r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" b1 "<<betaNumerator<<" d "<<d[3]<<" u "<<u[3]<<"\n";


		if (mVerboseLevel>2)
		{
			std::cout.precision(precision);
//			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm gradient squared "<<squaredNorm<<std::endl;
			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm gradient squared "<<betaNumerator<<std::endl;
			std::cout << std::setw(width)<< "pr " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< pr[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "r " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< r[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "h " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< h[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "d " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< d[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "u " ;
			for (size_t count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<<u[count] << "   " ;
			}
			std::cout << std::endl;

			std::cout << std::setw(width)<< "alpha "<<alpha<< "beta "<<beta << std::endl;
		}

//		if (squaredNorm<mAccuracyGradientScaled*mAccuracyGradientScaled)
		if (betaNumerator<mAccuracyGradientScaled*mAccuracyGradientScaled)
		{
			if (mVerboseLevel>2)
			{
				std::cout<< "CONVERGED " << std::endl;
			}
			converged = true;
			returnValue = NORMGRADIENT;
			break;
		}

		if (beta<0)
		{
			std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}


		if (mVerboseLevel>1 && curIteration%mShowSteps==0)
			std::cout<< "Iteration " << curIteration <<" with norm grad squared " << betaNumerator << std::endl;

		//increase iteration and curCycle
		curCycle++;
		curIteration++;

		if (curIteration>mMaxIterations)
		{
			converged = true;
			returnValue = MAXITERATIONS;
			break;
		}

		numFunctionCalls++;
		if (numFunctionCalls>mMaxFunctionCalls)
		{
			converged = true;
			returnValue = MAXFUNCTIONCALLS;
			break;
		}
	}
	isBuild = true;
	// Get force vector
//	r.clear();
//	CalculateReactionForcesEBE(u,r);
//	std::cout<<" Reaction forces \n";
//	for (size_t count=0; count<mNumParameters; count++)
//	{
//		std::cout << std::setw(width)<< r[count] << "\n   " ;
//	}
//	std::cout<<"\n";

#ifdef SHOW_TIME
    endOpt=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGridRed::Optimize] " << difftime(endOpt,startOpt)/CLOCKS_PER_SEC << "sec" << std::endl;
	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
 	outputTime<<(difftime(endOpt,startOpt)/CLOCKS_PER_SEC)<<"   "<<curIteration<<"\n";
	outputTime.close();

#endif
	std::cout<< "Number of Iterations............. " << curIteration << std::endl;
	if (mVerboseLevel>0)
	{
		std::cout<< " "  << std::endl;
		std::cout<< "Number of Function Calls......... " << numFunctionCalls << std::endl;
		std::cout<< "Number of Gradient Calls......... " << numGradientCalls << std::endl;
		std::cout<< "Number of Hessian Calls.......... " << numHessianCalls << std::endl;
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
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;
		std::cout.precision(precision);
		std::cout << std::setw(width)<< "displacements " ;
		for (size_t count=0; count<mNumParameters; count++)
		{
			std::cout << std::setw(width)<<u[count] << "   " ;

		}
		std::cout << std::endl;
	}
//	mParameters=u;

	return returnValue;
}

void NuTo::ConjugateGradientGridRed::CalcScalingFactors(int& numHessianCalls,myType &p)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
    //diagonal scaling with scaling factor
	++numHessianCalls;
    double scalefactor=0.0000000001;
	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
    outputTime<<scalefactor<<"  ";
    outputTime.close();
//    if(mVerboseLevel>0)
    std::cout<<"[ConjugateGradientGridRed::CalcScalingFactors] scale factor "<<scalefactor<<"\n";
    for (size_t count=0; count<mNumParameters; ++count)
        p[count] *=scalefactor;

//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGridRed::CalcScalingFactors] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//#endif
}

void NuTo::ConjugateGradientGridRed::HessianDiag(myType& rHessianDiag)const
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
	if (!mMatrixFreeMethod)
	{
		size_t numElems=mVoxelId.size();
		std::cout<<"[ConjugateGradientGridRed::HessianDiag] EBE"<<std::endl;
		std::vector<size_t> nodeNumbers(8);

		for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
		{
			double elemStiff=mYoungsModulus[mMaterialOfElem[elementNumber]]*mBaseStiffness[0];
#ifdef NODESATELEM
			nodeNumbers[node]=mNodesAtElem[8*elementNumber+node];
#else
			CalculateNodesAtElement(elementNumber,nodeNumbers);
#endif
			for(int node=0;node<8;++node)
			{
				rHessianDiag[3*nodeNumbers[node]]+=elemStiff;
				rHessianDiag[3*nodeNumbers[node]+1]+=elemStiff;
				rHessianDiag[3*nodeNumbers[node]+2]+=elemStiff;
			}
		}
	}
	else
	{
		std::cout<<"[ConjugateGradientGridRed::HessianDiag] NBN"<<std::endl;
		size_t numNodes=mEdgeId.size();
		for (size_t nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
		{
			double elemStiff[3]={0,0,0};
			for (size_t edge=28;edge<36;++edge)
			{
					elemStiff[0]+=mYoungsModulus[mMaterialOfElem[64*nodeNumber+edge]]*mEdgeStiffness[9*edge];
					elemStiff[1]+=mYoungsModulus[mMaterialOfElem[64*nodeNumber+edge]]*mEdgeStiffness[9*edge+4];
					elemStiff[2]+=mYoungsModulus[mMaterialOfElem[64*nodeNumber+edge]]*mEdgeStiffness[9*edge+8];
//					std::cout<<"elem Stiff of "<<nodeNumber <<" : "<<elemStiff[0]<<" "<<elemStiff[1]<<" "<<elemStiff[2]<<" \n";
			}
			rHessianDiag[3*nodeNumber]=elemStiff[0];
			rHessianDiag[3*nodeNumber+1]=elemStiff[1];
			rHessianDiag[3*nodeNumber+2]=elemStiff[2];
		}

	}

	for (size_t i=0;i<mNumParameters;++i)
	{
		if (rHessianDiag[i]!=0)
		{
					rHessianDiag[i]=1./rHessianDiag[i];
		}
	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGridRed::HessianDiag] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief ... calculate matrix-vector-product in element-by-element way
//! @brief ... with local vectors
//! @param ... u - prarmeters input, r - gradient output
void NuTo::ConjugateGradientGridRed::CalculateMatrixVectorProductEBE(myType &u,myType &r)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME

	size_t numElems=mVoxelId.size();
	int dofsElem=24;
	// global external force vector (active dofs)
    //	myType  force(mNumParameters,1);

	std::vector<double> residual (24);
	std::vector<double> displacement(24);
	std::vector<size_t> nodeNumbers(8);

	//loop over all elements
	for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
	{
		//calculate local return vector with all dofs: r=Ku
#ifndef NODESATELEM
		CalculateNodesAtElement(elementNumber,nodeNumbers);
#endif
		double youngsModulus=mYoungsModulus[mMaterialOfElem[elementNumber]];
		for(int node=0;node<8;++node)
		{
#ifdef NODESATELEM
			nodeNumbers[node]=mNodesAtElem[8*elementNumber+node];
#endif
			displacement[3*node]=u[3*nodeNumbers[node]];
			displacement[3*node+1]=u[3*nodeNumbers[node]+1];
			displacement[3*node+2]=u[3*nodeNumbers[node]+2];
		}
		for(int i=0;i<dofsElem;++i)
		{
			residual[i]=0;
			for(int j=0;j<dofsElem;++j)
				residual[i]+=youngsModulus*mBaseStiffness[dofsElem*i+j]*displacement[j];
		}
		for(int node=0;node<8;++node)
		{
			if (!mDofIsConstraint[3*nodeNumbers[node]])
				r[3*nodeNumbers[node]]+=residual[3*node];
			if (!mDofIsConstraint[3*nodeNumbers[node]+1])
				r[3*nodeNumbers[node]+1]+=residual[3*node+1];
			if (!mDofIsConstraint[3*nodeNumbers[node]+2])
				r[3*nodeNumbers[node]+2]+=residual[3*node+2];
//			std::cout<<"elem "<<elementNumber<<"  r["<<node<<"] = "<<r[3*nodeId] <<" "<<r[3*nodeId+1] <<" "<<r[3*nodeId+2] <<"\n";
		}
	}

//#ifdef NEIGHBORS
//	//loop over all elements
//	for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
//    {
//		//calculate local return vector with all dofs: r=Ku
//		for(int node=0;node<8;++node)
//		{
//			nodeId=mNodeId[mEdgesAtVoxel[8*mVoxelId[elementNumber]+node]];
//			displacement[3*node]=u[3*nodeId];
//			displacement[3*node+1]=u[3*nodeId+1];
//			displacement[3*node+2]=u[3*nodeId+2];
//		}
////		std::cout<<" mat "<<mMaterialOfElem[elementNumber];
////		std::cout<<" E "<<mYoungsModulus[mMaterialOfElem[elementNumber]];
////		std::cout<<" mBaseStiffness[0] "<<mBaseStiffness[0];
////		std::cout<<" mBaseStiffness[575] "<<mBaseStiffness[575]<<" \n";
//		for(int i=0;i<dofsElem;++i)
//		{
//			residual[i]=0;
//			for(int j=0;j<dofsElem;++j)
//			{
//				residual[i]+=mYoungsModulus[mMaterialOfElem[elementNumber]]*mBaseStiffness[dofsElem*i+j]*displacement[j];
//			}
//		}
//		for(int node=0;node<8;++node)
//		{
//			nodeId=mNodeId[mEdgesAtVoxel[8*mVoxelId[elementNumber]+node]];
//			if (!mDofIsConstraint[3*nodeId])
//				r[3*nodeId]+=residual[3*node];
//			if (!mDofIsConstraint[3*nodeId+1])
//				r[3*nodeId+1]+=residual[3*node+1];
//			if (!mDofIsConstraint[3*nodeId+2])
//				r[3*nodeId+2]+=residual[3*node+2];
////			std::cout<<"elem "<<elementNumber<<"  r["<<node<<"] = "<<r[3*nodeId] <<" "<<r[3*nodeId+1] <<" "<<r[3*nodeId+2] <<"\n";
//		}
//    }
//#else
//#endif
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGridRed::CalculateMatrixVectorProductEBE] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

void NuTo::ConjugateGradientGridRed::CalculateReactionForcesEBE(myType &u,myType &f)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME

	size_t numElems=mElemExist.size();
	int dofsElem=24;
	// global external force vector (active dofs)
    //	myType  force(mNumParameters,1);
	int nodeId=0;
	int numStiff=0;

	std::vector<double> residual (24);
	std::vector<double> displacement(24);


	//loop over all elements
	for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
    {
		if (mElemExist[elementNumber])
		{
			//get the number of the material
			numStiff = mMaterialOfElem[elementNumber];
			//calculate local return vector with all dofs: f=Ku
			for(int node=0;node<8;++node)
			{
				nodeId=mNodesAtElem[8*elementNumber+node];
				displacement[3*node]=u[3*nodeId];
				displacement[3*node+1]=u[3*nodeId+1];
				displacement[3*node+2]=u[3*nodeId+2];
			}
			for(int i=0;i<dofsElem;++i)
			{
				residual[i]=0;
				for(int j=0;j<dofsElem;++j)
					residual[i]+=mYoungsModulus[numStiff]*mBaseStiffness[dofsElem*i+j]*displacement[j];
			}
			for(int node=0;node<8;++node)
			{
				nodeId=mNodesAtElem[8*elementNumber+node];
				f[3*nodeId]+=residual[3*node];
				f[3*nodeId+1]+=residual[3*node+1];
				f[3*nodeId+2]+=residual[3*node+2];
//			std::cout<<"elem "<<elementNumber<<"  f["<<node<<"] = "<<f[3*nodeId] <<" "<<f[3*nodeId+1] <<" "<<f[3*nodeId+2] <<"\n";
			}
			// Variants-Test:
			// if-clauses after calculation
			// local u and f faster, because of cache
		}
	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGridRed::CalculateMatrixVectorProductEBE] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief ... calculate matrix-vector-product in node-by-node way
//! @brief ... with local vectors
//! @param ... u - prarmeters input, r - gradient output
void NuTo::ConjugateGradientGridRed::CalculateMatrixVectorProductNBN(myType &u,myType &r)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME
//	std::cout<<"[ConjugatGradientGrid::CalculateMatrixVectorProductNBN]\n";
	size_t numNodes=mEdgeId.size();
//	int  nodeNumbers[64];
	const int nodeNbrOfEdge[64]={0,1,1,2,3,3,4,4,4,4,5,5,6,7,7,8,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,15,15,16,16,16,16,17,17,18,19,19,20,21,21,22,22,22,22,23,23,24,25,25,26};
//	const int nbrOfEdges[27]={1,2,1,2,4,2,1,2,1,2,4,2,4,8,4,2,4,2,1,2,1,2,4,2,1,2,1};
	//loop over all elements
//	std::cout<<" NBN : numNodes "<<numNodes<<"\n";
	for (size_t node=0;node<numNodes;++node)
    {
//		double* youngsModulus=&mYoungsModulus[mMaterialOfElem[64*node]];

		size_t neighbors[64] ;
		double youngsModulus[64];
		for (size_t i=0;i<64;++i)
		{
			neighbors[i]=mNodeId[mEdgeId[node]+mNeighborNodes[nodeNbrOfEdge[i]]];
			youngsModulus[i] = mYoungsModulus[mMaterialOfElem[64*node+i]];
		}


//		std::cout<<" NBN : grid "<<nodeNumber<<" nodeid "<<node<<" grid neigh: ";
		r[3*node+0]=0;
		r[3*node+1]=0;
		r[3*node+2]=0;
		//27 neighbor nodes
		//calculate local return vector: r=Ku
		for(int i=0;i<64;++i)
		{
//					r[3*node+0]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]];
//					r[3*node+0]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+1]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+1];
//					r[3*node+0]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+2]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+2];
//
//					r[3*node+1]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+3]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]];
//					r[3*node+1]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+4]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+1];
//					r[3*node+1]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+5]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+2];
//
//					r[3*node+2]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+6]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]];
//					r[3*node+2]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+7]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+1];
//					r[3*node+2]+=mYoungsModulus[mMaterialOfElem[64*node+i]]*mEdgeStiffness[9*i+8]*u[3*mNodeId[nodeNumber+mNeighborNodes[nodeNbrOfEdge[i]]]+2];

//					r[3*node+0]+=*youngsModulus++*mEdgeStiffness[9*i]*u[3*neighbors[i]];
					r[3*node+0]+=youngsModulus[i]*mEdgeStiffness[9*i]*u[3*neighbors[i]];
					r[3*node+0]+=youngsModulus[i]*mEdgeStiffness[9*i+1]*u[3*neighbors[i]+1];
					r[3*node+0]+=youngsModulus[i]*mEdgeStiffness[9*i+2]*u[3*neighbors[i]+2];

					r[3*node+1]+=youngsModulus[i]*mEdgeStiffness[9*i+3]*u[3*neighbors[i]];
					r[3*node+1]+=youngsModulus[i]*mEdgeStiffness[9*i+4]*u[3*neighbors[i]+1];
					r[3*node+1]+=youngsModulus[i]*mEdgeStiffness[9*i+5]*u[3*neighbors[i]+2];

					r[3*node+2]+=youngsModulus[i]*mEdgeStiffness[9*i+6]*u[3*neighbors[i]];
					r[3*node+2]+=youngsModulus[i]*mEdgeStiffness[9*i+7]*u[3*neighbors[i]+1];
					r[3*node+2]+=youngsModulus[i]*mEdgeStiffness[9*i+8]*u[3*neighbors[i]+2];
		}
		r[3*node+0]*=!mDofIsConstraint[3*node];
		r[3*node+1]*=!mDofIsConstraint[3*node+1];
		r[3*node+2]*=!mDofIsConstraint[3*node+2];

//		std::cout<<"\n";
	}


//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGridRed::CalculateMatrixVectorProductNBN] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief ... calculate start gradient in node-by-node way
//! @brief ... variante II: with 3x3 matrix at node
void NuTo::ConjugateGradientGridRed::CalculateStartGradientNodeByNode(myType &r)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradientNodeByNodeII"<<std::endl;
//	myType gradientNew(mNumParameters,1);
//	int numNodes=mpGrid->GetNumNodes();
//	//nodeIds here ids of all neighbor nodes
//    int* nodeIds;
//	//array with global dof number of each dof of this node
//    int numDofs=3; //for dofs array needed
//
//	// part of matrix one dimensional 9 fields
//	double *array=new double [9*27];
//
//	NodeGrid3D* thisNode=0;
//	// global external force vector (active dofs)
////	myType  force(mpGrid->GetNumDofs(),1);
//	myType locReturn(numDofs,1);
//
//	//loop over all nodes
//	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
//    {
//		// return vector with all dofs of one node
//		//get pointer to each node
//		thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
//		if(thisNode)
//		{
//			locReturn.Resize(numDofs,1);
//			//get belonging nodes
//			nodeIds=thisNode->GetNodeIds();
//	//				std::cout<<" node "<<nodeNumber <<""<<std::endl;
//			for (int node=0;node<27;++node)
//			{
//				//myType* matrix;
//				//get pointer to array of part coefficient matrix for all neighbor nodes
//				array= thisNode->GetPartCoefficient0(0);
//	//			matrix= thisNode->GetPartCoefficientMatrix0(node);
//				//node exists
//		/*
//				if (nodeNumber==62)
//				{
//					std::cout<<" node "<<nodeNumber <<"ids :"<<nodeIds[0]<<"  "<<nodeIds[1]<< std::endl;
//				}
//			*/
//				if (nodeIds[node]>=0)
//				{
//					//myType locDispValues(3,1);
//					//mpGrid->NodeGetDisplacements(nodeIds[node],locDispValues);
//					//calculate nodal return vector and sum up
//					//when global dof is active
//					//when dof is not constraint, then ...
//					//subtract (r=f-Ku) locReturn for active dofs
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+0))
//						gradientNew(nodeNumber*numDofs+0,0)+=array[node*9+0]*r(nodeIds[node]*numDofs,0) + array[node*9+1]*r(nodeIds[node]*numDofs+1,0) + array[node*9+2]*r(nodeIds[node]*numDofs+2,0);
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+1))
//						gradientNew(nodeNumber*numDofs+1,0)+=array[node*9+3]*r(nodeIds[node]*numDofs,0) + array[node*9+4]*r(nodeIds[node]*numDofs+1,0) + array[node*9+5]*r(nodeIds[node]*numDofs+2,0);
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+2))
//						gradientNew(nodeNumber*numDofs+2,0)+=array[node*9+6]*r(nodeIds[node]*numDofs,0) + array[node*9+7]*r(nodeIds[node]*numDofs+1,0) + array[node*9+8]*r(nodeIds[node]*numDofs+2,0);
//		//				locReturn.operator += (matrix->operator *(locDispValues));
//				/*
//					if (nodeNumber==62)
//					{
//						std::cout<<" node "<<nodeNumber <<" dof 0 locReturn :"<<locReturn(0,0)<<" at loc node "<<nodeIds[node]<< std::endl;
//						std::cout<<" gradOrig :"<<r(nodeIds[node]*numDofs,0)<<" , "<< r(nodeIds[node]*numDofs+1,0)<<" , "<<r(nodeIds[node]*numDofs+2,0)<<std::endl;
//					}
//					*/
//				}
//
//			}
//		}
//    }
//	r=gradientNew;
//		//get global external force vector
//		//! @TODO add load vector
//		// ubdate this routine
//		//mpGrid->BuildGlobalExternalLoadVector(force);
//		//add global external force vector
//		//r+=force;

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGridRed::CalculateStartGradientNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGridRed::CalculateStartGradientNodeByNodeII] Modul Optimize is not loaded." );
#endif // ENABLE_MECHANICS
}



#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGridRed::Save ( const std::string &filename, std::string rType)const
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
			throw MathException ( "[ConjugateGradientGridRed::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception &e )
	{
		std::string s ( std::string ( "[ConjugateGradientGridRed::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
void NuTo::ConjugateGradientGridRed::Restore ( const std::string &filename,  std::string rType)
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
                throw OptimizeException ( "[NuTo::ConjugateGradientGridRed::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
                throw OptimizeException ( "[NuTo::ConjugateGradientGridRed::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
                throw OptimizeException ( "[NuTo::ConjugateGradientGridRed::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

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
//! @return    class name ConjugateGradientGridRed
std::string NuTo::ConjugateGradientGridRed::GetTypeId()const
{
    return std::string("ConjugateGradientGridRed");
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ConjugateGradientGridRed::Info () const
{
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}
void  NuTo::ConjugateGradientGridRed::AnsysInput(size_t numNodes,std::vector<size_t>& nodeId,boost::dynamic_bitset<> &rDofIsConstraint,std::vector<double>& youngsModulus,size_t* rGridDimension,double* rVoxelSpacing,std::vector<int>& materialOfElem,std::vector<int>& allNodesAtElem,std::vector<double>& parameters)
{
	// open file
	std::ofstream file;
    file.open("ansysInput");
    file<<"!Ansys Input File: \n DIM: "<<rGridDimension[0]<<" DOFS: "<<numNodes*3<<"\n";
    file<<"/prep7 \net,1,solid185 \nkeyopt,1,2,3 \n";

   int countMat=youngsModulus.size();
    for(int i=1;i<countMat;++i)
    {
        file<<"mp,ex,"<<i<<","<<youngsModulus[i]<<" \n";
        file<<"mp,prxy,"<<i<<",.2 \n";

    }
    int node=0;
    for (size_t z=0;z<rGridDimension[2]+1;++z)
    {
    	for (size_t y=0;y<rGridDimension[1]+1;++y)
    	{
    		for (size_t x=0;x<rGridDimension[0]+1;++x)
			{
    			if(nodeId[node]<numNodes)
    			{
    				file<<"n,"<<nodeId[node]+1<<","<<x*rVoxelSpacing[0]<<","<<y*rVoxelSpacing[1]<<","<<z*rVoxelSpacing[2]<<"\n";
   				if (rDofIsConstraint[3*nodeId[node]])
    				{
    					file<<"nsel,s,node,,"<<nodeId[node]+1<<"\n";
    				   	file<<"d,all,ux,"<<parameters[3*nodeId[node]]<<"\n";
     				}
   				if (rDofIsConstraint[3*nodeId[node]+1])
    				{
    					file<<"nsel,s,node,,"<<nodeId[node]+1<<"\n";
    				   	file<<"d,all,uy,"<<parameters[3*nodeId[node]+1]<<"\n";
     				}
   				if (rDofIsConstraint[3*nodeId[node]+2])
    				{
    					file<<"nsel,s,node,,"<<nodeId[node]+1<<"\n";
    				   	file<<"d,all,uz,"<<parameters[3*nodeId[node]+2]<<"\n";
     				}
   			}
    			node++;
			}
    	}
    }
    file<<"alls\n";
    int mat=1;
    file<<"type,1 \n mat,"<<mat<<" \n";
    int numElems=(int) materialOfElem.size();
    // attention for more materials
    for (int i=0;i<numElems;++i)
    {
		if(materialOfElem[i]!=mat)
		{
			mat=materialOfElem[i]+1;
			file<<"mat,"<<mat<<"\n";
		}
		file << "en," << i + 1 << "," << allNodesAtElem[8 * i] + 1
				<< "," << nodeId[allNodesAtElem[8 * i + 1]] + 1 << ","
				<< allNodesAtElem[8 * i + 2] + 1 << ","
				<< allNodesAtElem[8 * i + 3] + 1 << ","
				<< allNodesAtElem[8 * i + 4] + 1 << ","
				<< allNodesAtElem[8 * i + 5] + 1 << ","
				<< allNodesAtElem[8 * i + 6] + 1 << ","
				<< allNodesAtElem[8 * i + 7] + 1 << "\n";

    }
    file.close();
}

//! @brief ..calculate node numbers at one element
void NuTo::ConjugateGradientGridRed::CalculateNodesAtElement(size_t elementNumber,std::vector<size_t>& nodeNumbers)const
{
	size_t rVoxelLocation[3];
	size_t residual1=mVoxelId[elementNumber]%((mGridDimension[0])*(mGridDimension[1]));
	rVoxelLocation[0]=residual1%(mGridDimension[0]);
	rVoxelLocation[1]=residual1/(mGridDimension[0]);
	rVoxelLocation[2]=mVoxelId[elementNumber]/((mGridDimension[0])*(mGridDimension[1]));

	nodeNumbers[0] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[1]+1) + rVoxelLocation[0]];
	nodeNumbers[1] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[1]+1) + rVoxelLocation[0]+1];
	nodeNumbers[2] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0]+1];
	nodeNumbers[3] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0]];
	nodeNumbers[4] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0]];
	nodeNumbers[5] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0]+1];
	nodeNumbers[6] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0]+1];
	nodeNumbers[7] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0]];

//	std::cout<<" nodesAtElem "<<elementNumber<<": "<<nodeNumbers[0]<<" "<<nodeNumbers[1]<<" "<<nodeNumbers[2]<<" "<<nodeNumbers[3]<<" "<<nodeNumbers[4]<<" "<<nodeNumbers[5]<<" "<<nodeNumbers[6]<<" "<<nodeNumbers[7]<<"\n";


}
