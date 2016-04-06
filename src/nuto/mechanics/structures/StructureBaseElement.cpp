// $Id$

#ifdef SHOW_TIME
    #include <ctime>
#endif

# ifdef _OPENMP
    #include <omp.h>
# endif

#include <eigen3/Eigen/Eigenvalues>

#include <assert.h>
#include <boost/tokenizer.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputFullVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementOutputDummy.h"

//! @brief calls ElementCoefficientMatrix_0,
//! renaming only for clarification in mechanical problems for the end user
void NuTo::StructureBase::ElementStiffness(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult ,
		NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsRow,
		NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsColumn)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	try
	{
		ElementCoefficientMatrix(rElementId, 0, rResult, rGlobalDofsRow, rGlobalDofsColumn);
	}
	catch(NuTo::MechanicsException& e)
	{
		std::stringstream ss;
		ss << rElementId;
		std::string s = ss.str(); //Gets you a C++ STL string
		e.AddMessage("[NuTo::StructureBase::ElementStiffness] Error calculating stiffness of element "
				+ s + ".");
		throw e;
	}
	catch (...)
	{
		std::stringstream ss;
		ss << rElementId;
		std::string s = ss.str(); //Gets you a C++ STL string
		throw MechanicsException("[NuTo::StructureBase::ElementStiffness] Error calculating stiffness of element "
				+ s + ".");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementStiffness] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! and compares it to the matrix using central differences
//! for a mechanical problem, this corresponds to the stiffness matrix
//! @param rDelta  delta step for finite differences
//! @return element with maximum error
int NuTo::StructureBase::ElementTotalCoefficientMatrix_0_Check(double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDifference)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> tmpDifference;
    double maxError(0);
    int maxElement(-1);
	std::vector<int> globalDofsRow,globalDofsColumn;
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessAnalytic;
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessCDF;
//	bool symmetryFlag;
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	std::cout << "Element Id " << this->ElementGetId(elementVector[countElement]) << "\n" << "\n";

//        	elementVector[countElement]->CalculateCoefficientMatrix_0(stiffnessAnalytic, globalDofsRow, globalDofsColumn, symmetryFlag);
        	std::cout << "stiffnessAnalytic " << "\n" << stiffnessAnalytic << "\n" << "\n";

        	stiffnessCDF.Resize(stiffnessAnalytic.GetNumRows(),stiffnessAnalytic.GetNumColumns());
        	this->ElementCoefficientMatrix_0_Resforce(elementVector[countElement],rDelta,stiffnessCDF);
        	std::cout << "stiffnessCDF " << "\n" << stiffnessCDF << "\n" << "\n";

        	//check the maximum error
        	tmpDifference = (stiffnessCDF-stiffnessAnalytic)*(1./stiffnessAnalytic.Abs().Max());
        	std::cout << "difference "<< "\n" << tmpDifference <<"\n" << "\n";
        	double curError =  tmpDifference.Abs().Max();
        	std::cout << "error " << curError << "\n";
        	if (curError>maxError)
        	{
        		maxError = curError;
        		maxElement = countElement;
        		rDifference = tmpDifference;
        	}
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalCoefficientMatrix_0_Check] Error checking stiffness for "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalCoefficientMatrix_0_Check] Error checking stiffness for "
        			   + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalCoefficientMatrix_0_Check] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return this->ElementGetId(elementVector[maxElement]);
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! and compares it to the matrix using central differences
//! for a mechanical problem, this corresponds to the stiffness matrix
//! @param rElementId element
//! @param rDelta  delta step for finite differences
//! @return maximum difference between analytical and central difference method
double NuTo::StructureBase::ElementCoefficientMatrix_0_Check(int rElementId, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDifference)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0_Check] First update of tmp static data required.");
    }

//    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    double maxError;


    try
    {
    	NuTo::FullVector<int, Eigen::Dynamic> globalDofsRow,globalDofsColumn;
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessAnalytic;
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessCDF;

    	this->ElementCoefficientMatrix(rElementId, 0, stiffnessAnalytic, globalDofsRow, globalDofsColumn);
    	globalDofsRow.Info();
    	globalDofsColumn.Info();
    	std::cout << "stiffnessAnalytic " << "\n" << stiffnessAnalytic << "\n" << "\n";

    	stiffnessCDF.Resize(stiffnessAnalytic.GetNumRows(),stiffnessAnalytic.GetNumColumns());


    	this->ElementCoefficientMatrix_0_Resforce(rElementId,rDelta,stiffnessCDF);
    	std::cout << "stiffnessCDF " << "\n" << stiffnessCDF << "\n" << "\n";

    	//check the maximum error
    	rDifference = (stiffnessCDF-stiffnessAnalytic)*(1./stiffnessAnalytic.Abs().Max());
    	std::cout << "rDifference "<< "\n" << rDifference <<"\n" << "\n";
        maxError =  rDifference.Abs().Max();
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
    	e.AddMessage("[NuTo::StructureBase::ElementCoefficientMatrix_0_Check] Error checking element matrix for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementCoefficientMatrix_0_Check] Error checking element matrix for element " + ss.str() + ".");
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementCoefficientMatrix_0_Check] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return maxError;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! @param rElementId elementId
//! @param rDelta  delta step for finite differences
//! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
//! @return maximum difference between analytical and central difference method
void NuTo::StructureBase::ElementCoefficientMatrix_0_Resforce(int rElementId, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& stiffnessCDF)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    ElementCoefficientMatrix_0_Resforce(elementPtr,rDelta,stiffnessCDF);

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementCoefficientMatrix_0_Resforce] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! @param rElementPtr element
//! @param rDelta  delta step for finite differences
//! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
//! @return maximum difference between analytical and central difference method
void NuTo::StructureBase::ElementCoefficientMatrix_0_Resforce(ElementBase* rElementPtr, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& stiffnessCDF)
{
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> resforce1;
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> resforce2;

	// if other DOFs than displacements are implemented, a routine should be implemented for all elements, that
	// specifies wich DOFs are used

	int numNonlocalElements = rElementPtr->GetNumNonlocalElements();

	std::vector<const NuTo::ElementBase*>  nonlocalElements;
	if (numNonlocalElements==0 || stiffnessCDF.GetNumRows()==stiffnessCDF.GetNumColumns())
	{
		//local formulation
		nonlocalElements.push_back(rElementPtr);
	}
	else
	{
		nonlocalElements = rElementPtr->GetNonlocalElements();
	}

	//determine the initial trial resforce vector
	// if stiffnessAnalytic is square with nonlocal!=0, then all ips are elastic
	std::vector<int> globalDofsRow;
	//update tmpstatic data of all nonlocal elements
	if (mHaveTmpStaticData)
	{
		for (unsigned int countElement2=0; countElement2<nonlocalElements.size(); countElement2++)
		{
//            const_cast<NuTo::ElementBase*>(nonlocalElements[countElement2])->UpdateStaticData(NuTo::Element::TMPSTATICDATA);
		}
	}
//	rElementPtr->CalculateGradientInternalPotential(resforce1,globalDofsRow);
	//std::cout << "resforce 1" << std::endl;
	//resforce1.Trans().Info(12,10);

	int curCol(0);

	for (unsigned int countElement=0; countElement<nonlocalElements.size(); countElement++)
	{
		for (int countNode=0; countNode<nonlocalElements[countElement]->GetNumNodes(); countNode++)
		{
			NodeBase *theNode=const_cast<NuTo::ElementBase*>(nonlocalElements[countElement])->GetNode(countNode);
			Eigen::Matrix<double, Eigen::Dynamic, 1> disp;
			int numDisp = theNode->GetNumDisplacements();
			assert(numDisp<=3);

			for (int countDisp=0; countDisp<numDisp; countDisp++)
			{
			    if (theNode->GetNumDisplacements() != 1 and theNode->GetNumDisplacements() != 2 and theNode->GetNumDisplacements() != 3)
                    throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0_Resforce] Only nodes with 1,2 or 3 displacement components considered.");

			    disp = theNode->GetDisplacements();
			    disp[countDisp]+=rDelta;
				theNode->SetDisplacements(disp);

				//update tmpstatic data of nonlocal elements
				if (mHaveTmpStaticData)
				{
//                     const_cast<NuTo::ElementBase*>(nonlocalElements[countElement])->UpdateStaticData(NuTo::Element::TMPSTATICDATA);
				}

				if (nonlocalElements[countElement]==rElementPtr)
				{
					//calculate new residual vector and afterwards reset the displacements
//					rElementPtr->CalculateGradientInternalPotential(resforce2, globalDofsRow);

					disp[countDisp]-=rDelta;
					theNode->SetDisplacements(disp);

					//update tmpstatic data of nonlocal elements
					if (mHaveTmpStaticData)
					{
//						const_cast<NuTo::ElementBase*>(nonlocalElements[countElement])->UpdateStaticData(NuTo::Element::TMPSTATICDATA);
					}
				}
				else
				{
					//reset the displacements and then calculate the resforce and the reset the tmp static data
					disp[countDisp]-=rDelta;
                    theNode->SetDisplacements(disp);

//					rElementPtr->CalculateGradientInternalPotential(resforce2, globalDofsRow);

					//update tmpstatic data of nonlocal elements
					if (mHaveTmpStaticData)
					{
//						const_cast<NuTo::ElementBase*>(nonlocalElements[countElement])->UpdateStaticData(NuTo::Element::TMPSTATICDATA);
					}

				}

				//std::cout << "resforce 2" << std::endl;
				//resforce2.Trans().Info(12,10);
				//assemble the matrix
				if (curCol>=stiffnessCDF.GetNumColumns())
				{
					std::cout << "num nonlocal elements " << nonlocalElements.size() << "\n";
					std::cout << "cur col " << curCol << " size of matrix " << stiffnessCDF.GetNumColumns() << "\n";
					throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0_Resforce] Allocated matrix for calculation of stiffness with finite differences has illegal size.");
				}
				stiffnessCDF.col(curCol) = (resforce2-resforce1)*(1./rDelta);
				curCol++;
			}
		}
	}
}

//! @brief Checks the global CoefficientMatrix_0 and the internal forces vector by the internal energy
//! @return false, if stiffness is not correct
bool NuTo::StructureBase::CheckStiffness()
{

    throw NuTo::MechanicsException("[NuTo::StructureBase::CheckStiffness] Check of internal forces vector currently. "
            " not correctly implemented. Use NuTo::StructureBase::CoefficientMatrix_0() to check the stiffness only.");

    //be carefull this routine performs a node merge, which modifies the displacements
    //as a result, the stiffness calculated here might be different from the one if you just call the stiffness routine
    //this is especially true for the first step of the Newton iteration in displacement control situation
    //where the stiffness of the old state is calulated on purpose and the multiplied by the difference between prescribed dependent dofs and actual dependent dofs

    mLogger << "test of stiffness still included, node merge is called!!! " << "\n";
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    NuTo::FullVector<double,Eigen::Dynamic> displacementsActiveDOFsCheck;
    NuTo::FullVector<double,Eigen::Dynamic> displacementsDependentDOFsCheck;

#ifdef SHOW_TIME
    bool oldShowtime = mShowTime;
    mShowTime = false;
#endif

    //recalculate stiffness
    this->NodeExtractDofValues(0,displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    //mLogger << "active dof values " << "\n";
    //displacementsActiveDOFsCheck.Trans().Info(12,4);
    //mLogger << "dependent dof values " << "\n";
    //displacementsDependentDOFsCheck.Trans().Info(12,4);
    this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();
    this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    //this->ConstraintInfo(10);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
    //std::cout<<"stiffness matrix" << "\n";
    //stiffnessMatrixCSRVector2Full.Info(10,3);
    double interval(-1e-10);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrixCSRVector2_CDF(stiffnessMatrixCSRVector2.GetNumRows(), stiffnessMatrixCSRVector2.GetNumColumns());
    NuTo::FullVector<double,Eigen::Dynamic> intForceVector1, intForceVector2, intForceVectorCDF(stiffnessMatrixCSRVector2.GetNumRows());
    double energy1,energy2;
    this->NodeExtractDofValues(0,displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();
    this->BuildGlobalGradientInternalPotentialVector(intForceVector1);
    //this->NodeInfo(10);
    energy1 = this->ElementTotalGetInternalEnergy();
    energy1 += this->ConstraintTotalGetTotalEnergy();
    //std::cout << "check stiffness:: energy1 "<<  energy1 << "\n";

    for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
    {
    	displacementsActiveDOFsCheck(count,0)+=interval;
        this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
        this->ElementTotalUpdateTmpStaticData();
        this->BuildGlobalGradientInternalPotentialVector(intForceVector2);
        //std::cout << "check stiffness:: intForceVector2"<< "\n";
        //intForceVector2.Trans().Info(10,6);
        //this->ConstraintInfo(10);
        energy2 = this->ElementTotalGetInternalEnergy();
        energy2 += this->ConstraintTotalGetTotalEnergy();
        stiffnessMatrixCSRVector2_CDF.SetColumn(count,(intForceVector2-intForceVector1)*(1./interval));
        intForceVectorCDF(count,0) = (energy2-energy1)/interval;
        displacementsActiveDOFsCheck(count,0)-=interval;
    }
    this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();

#ifdef SHOW_TIME
    mShowTime=oldShowtime;

#endif
    if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).cwiseAbs().maxCoeff()>1e-1)
    {
        if (stiffnessMatrixCSRVector2Full.GetNumRows()<100)
        {
			mLogger << "globalStiffnessMatrix algo" << "\n";
			mLogger << stiffnessMatrixCSRVector2Full;
			mLogger << "\n" << "globalStiffnessMatrix cdf" << "\n";
			mLogger << stiffnessMatrixCSRVector2_CDF;
			mLogger << "\n" << "error" << "\n";
			mLogger << (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full);
        }
        //extract the first 5x5 block
        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> blockAlgo(stiffnessMatrixCSRVector2Full.GetBlock(0,0,5,5));
        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> blockCDF(stiffnessMatrixCSRVector2_CDF.GetBlock(0,0,5,5));
        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> blockDelta(blockAlgo-blockCDF);

        double maxError;
        int row,col;
        maxError =(stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).cwiseAbs().maxCoeff(&row,&col);
        mLogger << "maximum error stiffness is " << maxError << " at (" << row << "," << col << ") with abs value in correct matrix " << stiffnessMatrixCSRVector2Full(row,col) << "\n";

        if (stiffnessMatrixCSRVector2Full.GetNumRows()<100)
        {
			mLogger << "\n" << "intForceVector algo" << "\n";
			mLogger << intForceVector1.Trans();
			mLogger << "\n" << "intForceVector cdf" << "\n";
			mLogger << intForceVectorCDF.Trans();
			mLogger << "\n" << "error" << "\n";
			mLogger << NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(intForceVector1-intForceVectorCDF).cwiseAbs().transpose();
        }
        maxError = (intForceVector1-intForceVectorCDF).cwiseAbs().maxCoeff(&row,&col);
        mLogger << "maximum error resforce is " << maxError << " at (" << row << "," << col << ") " << "\n";

        //throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).cwiseAbs().maxCoeff()>1e1)
        {
        	NodeInfo(10);
            mLogger << "stiffness ist wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<< "\n";
            exit(0);
        }
        return false;
    }
    else
    {
        mLogger << "stiffness is OK "<< "\n";
        return true;
    }

}

//! @brief Compares the global CoefficientMatrix_0 to the matrix using central differences
//! @param rDelta delta step for finite differences
//! @param rPrintResult 'true' prints the result
bool NuTo::StructureBase::CheckCoefficientMatrix_0(double rDelta, bool rPrintResult)
{
    //be carefull this routine performs a node merge, which modifies the displacements
    //as a result, the stiffness calculated here might be different from the one if you just call the stiffness routine
    //this is especially true for the first step of the Newton iteration in displacement control situation
    //where the stiffness of the old state is calulated on purpose and the multiplied by the difference between prescribed dependent dofs and actual dependent dofs

    // ttitsche: deleted all the internal energy stuff, since it is not implemented.

#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

//    mLogger << "test of stiffness still included, node merge is called!!! " << "\n";
    NuTo::SparseMatrixCSRVector2General<double> coeffMatrix_CSR;
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    NuTo::FullVector<double,Eigen::Dynamic> displacementsActiveDOFsCheck;
    NuTo::FullVector<double,Eigen::Dynamic> displacementsDependentDOFsCheck;

#ifdef SHOW_TIME
    bool oldShowtime = mShowTime;
    mShowTime = false;
#endif

    // calculate stiffness
    this->NodeExtractDofValues(0,displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();
    this->BuildGlobalCoefficientMatrix0(coeffMatrix_CSR, dispForceVector);


    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> coeffMatrix(coeffMatrix_CSR);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> coeffMatrix_CDF(mNumActiveDofs, mNumActiveDofs);
    NuTo::FullVector<double,Eigen::Dynamic> gradIntPotential1, gradIntPotential2, gradIntPotentialMod;

    this->ElementTotalUpdateTmpStaticData();
    this->BuildGlobalGradientInternalPotentialVector(gradIntPotential1);

    // if this routine is called before the boundary conditions are applied,
    // the dispForceVector will be not equal to zero.
    // during the equilibrium iterations, it will be zero.
    gradIntPotentialMod = gradIntPotential1 - dispForceVector;

    for (int iCol=0; iCol<displacementsActiveDOFsCheck.GetNumRows(); iCol++)
    {
        displacementsActiveDOFsCheck(iCol)+=rDelta;
        this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);

        this->ElementTotalUpdateTmpStaticData();

        this->BuildGlobalGradientInternalPotentialVector(gradIntPotential2);

        coeffMatrix_CDF.SetColumn(iCol,(gradIntPotential2-gradIntPotentialMod)*(1./rDelta));




        displacementsActiveDOFsCheck(iCol)-=rDelta;
        // the NodeMerge is called in the next loop-step and
        // (for the last step) right after this loop
    }
    this->NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();

#ifdef SHOW_TIME
    mShowTime=oldShowtime;
#endif

    int row,col;
    double maxError =(coeffMatrix_CDF-coeffMatrix).cwiseAbs().maxCoeff(&row,&col);


    bool stiffnesIsCorrect = maxError < 1.e-1;

    if (not stiffnesIsCorrect)
    {
        mLogger << "[NuTo::StructureBase::CheckCoefficientMatrix_0] maximum error stiffness is " << maxError << " at (" << row << "," << col << ") with abs value in correct matrix " << coeffMatrix(row,col) << "\n";

        if (mNumActiveDofs < 100 and rPrintResult)
        {
            mLogger << "global coefficient matrix algo" << "\n";
            mLogger << coeffMatrix;
            mLogger << "\n";
            mLogger << "global coefficient matrix cdf" << "\n";
            mLogger << coeffMatrix_CDF;
            mLogger << "\n";
            mLogger << "error" << "\n";
            mLogger << (coeffMatrix_CDF-coeffMatrix);
            NodeInfo(10);
        }
        mLogger << "[NuTo::StructureBase::CheckCoefficientMatrix_0]               !!! WRONG !!! "<< "\n";
    }
    else
    {
        mLogger << "[NuTo::StructureBase::CheckCoefficientMatrix_0]               !!!  OK  !!! "<< "\n";
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::CheckCoefficientMatrix_0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif

    return stiffnesIsCorrect;
}



//! @brief Compares the element CoefficientMatrix_0 to the matrix using central differences
//! for all elements
//! @param rDelta delta step for finite differences
bool NuTo::StructureBase::ElementCheckCoefficientMatrix_0(double rDelta)
{

#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
    bool oldShowtime = mShowTime;
    mShowTime = false;
#endif

    // get all element pointers
    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    int numElements = elementVector.size();

    // store elementIds and their max difference in vectors for output
    NuTo::FullVector<double, Eigen::Dynamic> elementDiffs(numElements);
    NuTo::FullVector<int, Eigen::Dynamic> elementIds(numElements);

    bool stiffnessIsCorrectTotal = true;

    // check every element
    for (int iElement = 0; iElement < numElements; ++iElement) {

        int elementId = this->ElementGetId(elementVector[iElement]);
        elementIds(iElement) = elementId;

        // check specific element with elementId
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> difference;
        bool stiffnessIsCorrect = this->ElementCheckCoefficientMatrix_0(rDelta, elementId, difference, true);

        // a single incorrect stiffness causes the total variable to become and stay false
        stiffnessIsCorrectTotal = stiffnessIsCorrectTotal and stiffnessIsCorrect;

        // calculate maximal difference
        double maxDiff = difference.cwiseAbs().maxCoeff();
        elementDiffs(iElement) = maxDiff;
    }

#ifdef SHOW_TIME
    mShowTime=oldShowtime;
#endif

    int maxElementIndex, maxCol;
    double maxCoeff = elementDiffs.Max(maxElementIndex, maxCol);

    if (stiffnessIsCorrectTotal)
    {
        mLogger << "[NuTo::StructureBase::ElementCheckCoefficientMatrix_0]        !!!  OK  !!! "<< "\n";
    } else {
        mLogger << "[NuTo::StructureBase::ElementCheckCoefficientMatrix_0] max diff = "
                << maxCoeff << " in elementId "<< elementIds(maxElementIndex) << "\n";
        mLogger << "[NuTo::StructureBase::ElementCheckCoefficientMatrix_0]        !!! WRONG !!! "<< "\n";
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementCheckCoefficientMatrix_0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif

    return stiffnessIsCorrectTotal;
}

//! @brief Compares the element CoefficientMatrix_0 to the matrix using central differences
//! for a specific element
//! @param rDelta delta step for finite differences
//! @param rElementId element id
//! @param rDifference difference to CDF solution
bool NuTo::StructureBase::ElementCheckCoefficientMatrix_0(
        double rDelta, int rElementId,
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rDifference,
        bool rPrintResult)
{

    NuTo::FullVector<double, Eigen::Dynamic> gradientInternal, gradientInternal2;
    NuTo::FullVector<int, Eigen::Dynamic> globalRowDofs;
    NuTo::FullVector<int, Eigen::Dynamic> globalColDofs;

    // calculate coefficient matrix
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffMatrix;
    this->ElementCoefficientMatrix(rElementId, 0, coeffMatrix, globalRowDofs, globalColDofs);

    // get reference
    this->ElementGradientInternalPotential(rElementId, gradientInternal, globalRowDofs);

    // get active dof values
    NuTo::FullVector<double, Eigen::Dynamic> actDofValues, depDofValues;
    this->NodeExtractDofValues(actDofValues, depDofValues);

    int numRows = coeffMatrix.rows();
    int numCols = coeffMatrix.cols();

    // calculate coefficient matrix via CDF
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffMatrix_CDF(numRows, numCols);
    for (int iCol = 0; iCol < numRows; ++iCol) {
        // apply rDelta to the corresponding Dof
        int index = globalRowDofs(iCol);
        if (index < mNumActiveDofs)
            actDofValues(index) += rDelta;
        else
            depDofValues(index-mNumActiveDofs) += rDelta;

        this->NodeMergeDofValues(actDofValues, depDofValues);

        // this->ElementTotalUpdateTmpStaticData(); ??

        this->ElementGradientInternalPotential(rElementId, gradientInternal2, globalRowDofs);

        coeffMatrix_CDF.SetColumn(iCol, (gradientInternal2-gradientInternal)/rDelta);

        // restore the orignial dof values
        if (index < mNumActiveDofs)
            actDofValues(index) -= rDelta;
        else
            depDofValues(index-mNumActiveDofs) -= rDelta;

    }
    this->NodeMergeDofValues(actDofValues, depDofValues);

    rDifference = coeffMatrix - coeffMatrix_CDF;
    double maxDiff = rDifference.cwiseAbs().maxCoeff();

    bool stiffnessIsCorrect = maxDiff < 1.e-1;

    if (rPrintResult)
    {
        if (stiffnessIsCorrect)
        {
            mLogger << "[NuTo::StructureBase::ElementCheckCoefficientMatrix_0] Element ID: " << rElementId << " OK \n";
        }
        else
        {
            mLogger << "[NuTo::StructureBase::ElementCheckCoefficientMatrix_0] Element ID: " << rElementId << " WRONG \n";
            mLogger << "element coefficient matrix algo" << "\n";
            mLogger << coeffMatrix;
            mLogger << "\n";
            mLogger << "element coefficient matrix CDF" << "\n";
            mLogger << coeffMatrix_CDF;
            mLogger << "\n";
            mLogger << "error" << "\n";
            mLogger << rDifference;
        }
    }

    return stiffnessIsCorrect;
}

//! @brief calculates the coefficient matrix for the rTimeDerivative derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness, damping or mass matrix (rTimeDerivative=0,1,2)
void NuTo::StructureBase::ElementCoefficientMatrix(int rElementId,
		                        int rTimeDerivative,
                                NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
                                NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsRow,
                                NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsColumn)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix_0] First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    std::vector<int> globalDofsRow, globalDofsColumn;

    try
    {
    	//define output
    	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	switch(rTimeDerivative)
    	{
    	case 0:
        	boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>( elementOutput )( Element::HESSIAN_0_TIME_DERIVATIVE );
    		break;
    	case 1:
        	boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>( elementOutput )( Element::HESSIAN_1_TIME_DERIVATIVE );
    		break;
    	case 2:
        	boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>( elementOutput )( Element::HESSIAN_2_TIME_DERIVATIVE );
    		break;
    	default:
    		throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix] only time derivatives 0,1 and 2 supported.");
    	}
    	boost::assign::ptr_map_insert<ElementOutputVectorInt>( elementOutput )( Element::GLOBAL_ROW_DOF );
    	boost::assign::ptr_map_insert<ElementOutputVectorInt>( elementOutput )( Element::GLOBAL_COLUMN_DOF );


		//evaluate output
    	elementPtr->Evaluate(elementOutput);

    	//assign output
    	switch(rTimeDerivative)
    	{
    	case 0:
    		rResult = elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second->GetFullMatrixDouble();
    		break;
    	case 1:
    		rResult = elementOutput.find(Element::HESSIAN_1_TIME_DERIVATIVE)->second->GetFullMatrixDouble();
    		break;
    	case 2:
    		rResult = elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second->GetFullMatrixDouble();
    		break;
    	default:
    		throw MechanicsException("[NuTo::StructureBase::ElementCoefficientMatrix] only time derivatives 0,1 and 2 supported.");
    	}

		globalDofsRow = elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt();
		globalDofsColumn = elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt();
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
    	e.AddMessage("[NuTo::StructureBase::ElementCoefficientMatrix] Error building element matrix for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementCoefficientMatrix] Error building element matrix for element " + ss.str() + ".");
    }

    //cast to FullMatrixInt
    rGlobalDofsRow.Resize(globalDofsRow.size());
    memcpy(rGlobalDofsRow.data(),&globalDofsRow[0],globalDofsRow.size()*sizeof(int));

    rGlobalDofsColumn.Resize(globalDofsColumn.size());
    memcpy(rGlobalDofsColumn.data(),&globalDofsColumn[0],globalDofsColumn.size()*sizeof(int));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementCoefficientMatrix_0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::StructureBase::ElementGradientInternalPotential(int rElementId,
		NuTo::FullVector<double, Eigen::Dynamic>& rResult,
		NuTo::FullVector<int, Eigen::Dynamic>& rGlobalDofsRow)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGradientInternalPotential] First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
    	//define output
    	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>( elementOutput )( Element::INTERNAL_GRADIENT );
    	boost::assign::ptr_map_insert<ElementOutputVectorInt>( elementOutput )( Element::GLOBAL_ROW_DOF );

		//evaluate output
    	elementPtr->Evaluate(elementOutput);

    	//assign output
		rResult = elementOutput.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble();
		std::vector<int>& globalDofsRow = elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt();
		//cast to FullMatrixInt
	    rGlobalDofsRow.Resize(globalDofsRow.size());
	    memcpy(rGlobalDofsRow.data(),&globalDofsRow[0],globalDofsRow.size()*sizeof(int));
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
    	e.AddMessage("[NuTo::StructureBase::ElementGradientInternalPotential] Error building element vector for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGradientInternalPotential] Error buildung element vector for element " + ss.str() + ".");
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGradientInternalPotential] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the constitutive law of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    try
    {
    	ElementSetConstitutiveLaw(elementPtr,itConstitutive->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the constitutive law of a single element
//! @param rElementIdent identifier for the element
//! @param rIp  id of integration point
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId,int rIp, int rConstitutiveLawIdent)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    try
    {
    	ElementSetConstitutiveLaw(elementPtr,rIp,itConstitutive->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}
//! @brief sets the constitutive law of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetConstitutiveLaw(int rGroupIdent, int rConstitutiveLawIdent)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(itElement->second,itConstitutive->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law  for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law for element " + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the constitutive law of a all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(elementVector[countElement],itConstitutive->second);
        	if (elementVector[countElement]->GetNumNonlocalElements()>0)
        	    elementVector[countElement]->DeleteNonlocalElements();
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law for element "
        			   + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the constitutive law of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive)
{
	rElement->SetConstitutiveLaw(rConstitutive);
}

//! @brief sets the constitutive law of a single ip at an element
//! @param rElement element pointer
//! @param rIp number of integration point
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement,int rIp, ConstitutiveBase* rConstitutive)
{
    //std::cout<< "[NuTo::StructureBase::ElementSetConstitutiveLaw]" << "\n";
	rElement->SetConstitutiveLaw(rIp,rConstitutive);
}


//! @brief sets the section of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the section
void NuTo::StructureBase::ElementSetSection(int rElementId, int rSectionId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetSection] Section with the given identifier does not exist.");

    try
    {
    	ElementSetSection(elementPtr,itSection->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetSection] Error setting section for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetSection] Error setting section for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementSetSection] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the section of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetSection(int rGroupIdent, int rSectionId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Section with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetSection(itElement->second,itSection->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element " + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief sets the section for all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetSection(int rSectionId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetSection] Section with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetSection(elementVector[countElement],itSection->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetSection] Error setting section  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetSection] Error setting section for element "
        			   + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalSetSection] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief modifies the material of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetSection(ElementBase* rElement, SectionBase* rSection)
{
    rElement->SetSection(rSection);
}

//! @brief modifies the interpolation type of a single element
//! @param rElementId ... element number
//! @param rInterpolationTypeId ... interpolation type id
void NuTo::StructureBase::ElementSetInterpolationType(int rElementId, int rInterpolationTypeId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetInterpolationType] Interpolation type with the given identifier does not exist.");

    try
    {
        ElementSetInterpolationType(elementPtr, itInterpolationType->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetInterpolationType] Error setting interpolation type for element "
            + ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
        throw NuTo::MechanicsException
           ("[NuTo::StructureBase::ElementSetInterpolationType] Error setting interpolation type for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementSetInterpolationType] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif

}

//! @brief modifies the interpolation type of a group of elements
//! @param rGroupId ... identifier for the group of elements
//! @param rInterpolationTypeId ... interpolation type id
void NuTo::StructureBase::ElementGroupSetInterpolationType(int rGroupId, int rInterpolationTypeId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Interpolation type with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
            ElementSetInterpolationType(itElement->second, itInterpolationType->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetInterpolationType] Error setting interpolation type for element "
                + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupSetInterpolationType] Error setting interpolation type for element " + ss.str() + ".");
        }
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif

}

//! @brief modifies the interpolation type of a single element
//! @param rElement element pointer
//! @param rInterpolationType interpolation type
void NuTo::StructureBase::ElementSetInterpolationType(ElementBase* rElement, InterpolationType* rInterpolationType)
{
    rElement->SetInterpolationType(rInterpolationType);
}


//! @brief returns the enum of string identifier for an integration type
//! @param rIpDataTypeStr string
//! @return enum
NuTo::IpData::eIpDataType NuTo::StructureBase::ElementGetEnumIntegrationType(const std::string& rIpDataTypeStr)
{
    // get ip data type
    std::string upperCaseIpDataTypeStr;
    std::transform(rIpDataTypeStr.begin(), rIpDataTypeStr.end(), std::back_inserter(upperCaseIpDataTypeStr), (int(*)(int)) toupper);

    NuTo::IpData::eIpDataType ipDataType;
    if (upperCaseIpDataTypeStr=="NOIPDATA")
    {
    	ipDataType = NuTo::IpData::NOIPDATA;
    }
    else if (upperCaseIpDataTypeStr=="STATICDATA")
	{
    	ipDataType = NuTo::IpData::STATICDATA;
	}
    else if (upperCaseIpDataTypeStr=="STATICDATANONLOCAL")
    {
    	ipDataType = NuTo::IpData::STATICDATANONLOCAL;
    }
    else
    {
    	throw MechanicsException("[NuTo::Structure::ElementGetEnumIntegrationType] Ip data type "+upperCaseIpDataTypeStr +" does not exist.");
    }
    return ipDataType;
}


//! @brief calculates the engineering strain
//! @param rElemIdent  identifier for the element
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringStrain(int rElementId, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetEngineeringStrain] First update of tmp static data required.");
    }
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
		// define variables storing the element contribution
		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRAIN);

		//evaluate the strain
    	elementPtr->Evaluate(elementOutput);

    	rEngineeringStrain = elementOutput.find(Element::IP_DATA)->second->GetFullMatrixDouble();
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringStrain] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetEngineeringStrain] Error getting engineering strain for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGetEngineeringStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the engineering plastic strain
//! @param rElemIdent  identifier for the element
//! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringPlasticStrain(int rElementId, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringPlasticStrain)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] First update of tmp static data required.");
    }
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    try
    {
		// define variables storing the element contribution
		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_PLASTIC_STRAIN);

		//evaluate the strain
    	elementPtr->Evaluate(elementOutput);

    	rEngineeringPlasticStrain = elementOutput.find(Element::IP_DATA)->second->GetFullMatrixDouble();
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] Error getting engineering strain for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGetEngineeringPlasticStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the engineering stress
//! @param rElemIdent  identifier for the element
//! @param rEngineeringStress engineering stress (return value, always 6xnumIp matrix)
void NuTo::StructureBase::ElementGetEngineeringStress(int rElementId, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetEngineeringStress] First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
		// define variables storing the element contribution
		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRESS);

		//evaluate the strain
    	elementPtr->Evaluate(elementOutput);

    	rEngineeringStress = elementOutput.find(Element::IP_DATA)->second->GetFullMatrixDouble();
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetEngineeringStress] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetEngineeringStress] Error getting engineering strain for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGetEngineeringStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


//! @brief calculates the damage
//! @param rElemIdent  identifier for the element
//! @param rEngineeringStress damage (return value, always 1xnumIp matrix)
void NuTo::StructureBase::ElementGetDamage(int rElementId, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDamage)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetDamage] First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
		// define variables storing the element contribution
		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::DAMAGE);

		//evaluate the strain
    	elementPtr->Evaluate(elementOutput);

    	rDamage = elementOutput.find(Element::IP_DATA)->second->GetFullMatrixDouble();
     }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetDamage] Error getting engineering strain for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetDamage] Error getting engineering strain for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGetDamage] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the global integration point coordinates
//! @param rElemIdent  identifier for the element
//! @param rCoordinates integration point coordinates (return value, always 3xnumIp matrix)
void NuTo::StructureBase::ElementGetIntegrationPointCoordinates(int rElementId, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
		//evaluate the coordinates
    	rCoordinates.Resize(3,elementPtr->GetNumIntegrationPoints());
    	for (int count=0; count<elementPtr->GetNumIntegrationPoints(); count++)
    	{
    	    Eigen::Vector3d coords = elementPtr->GetGlobalIntegrationPointCoordinates(count);
    	    rCoordinates.SetBlock(0, count, coords);
    	}
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] Error getting integration point coordinates for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] Error getting integration point coordinates for element " + ss.str() + ".");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGetEngineeringStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the maximum damage in all elements
//! @param rElemIdent  identifier for the element
//! @return max damage value
double NuTo::StructureBase::ElementTotalGetMaxDamage()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementTotalGetMaxDamage] First update of tmp static data required.");
    }
	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);

	double maxDamage(0);
	// define variables storing the element contribution
	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::DAMAGE);

	for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
			elementVector[countElement]->Evaluate(elementOutput);
			FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& damage = elementOutput.find(Element::IP_DATA)->second->GetFullMatrixDouble();
			if (damage.Max()>maxDamage)
				maxDamage = damage.Max();
		}
		catch(NuTo::MechanicsException &e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			e.AddMessage("[NuTo::StructureBase::ElementTotalGetMaxDamage] Error getting damage for element "
				+ ss.str() + ".");
			throw e;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			throw NuTo::MechanicsException
			   ("[NuTo::StructureBase::ElementTotalGetMaxDamage] Error getting engineering strain for element " + ss.str() + ".");
		}
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalGetMaxDamage] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
	return maxDamage;
}

//! @brief updates the history data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementTotalUpdateStaticData()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        try
        {
            Error::eError error = this->ElementTotalUpdateTmpStaticData();
            if (error!=Error::SUCCESSFUL)
            	return error;
        }
        catch (NuTo::Exception& e)
        {
        	e.AddMessage("[NuTo::StructureBase::ElementTotalUpdateStaticData] error building tmp static data.");
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateStaticData] error building tmp static data.");
        }
    }

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
	Error::eError errorGlobal (Error::SUCCESSFUL);
    int exception(0);
    std::string exceptionStringTotal;

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputDummy>( elementOutput )( Element::UPDATE_STATIC_DATA );

#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
            Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
            if (error!=Error::SUCCESSFUL)
            {
            	if (errorGlobal==Error::SUCCESSFUL)
            		errorGlobal = error;
            	else if (errorGlobal!=error)
            		throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateStaticData] elements have returned multiple different error codes, can't handle that.");
            }
		}
		catch(NuTo::Exception& e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal(e.ErrorMessage()
					+"[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
	}
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalUpdateStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalUpdateStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return errorGlobal;
}

//! @brief updates the history data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementTotalUpdateTmpStaticData()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

    Error::eError errorGlobal (Error::SUCCESSFUL);

    //std::cout << "do we really have tmp static data " << mHaveTmpStaticData << std::endl;
    if (mHaveTmpStaticData)
	{
		std::vector<ElementBase*> elementVector;
		GetElementsTotal(elementVector);
	    int exception(0);
	    std::string exceptionStringTotal;

		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
		boost::assign::ptr_map_insert<ElementOutputDummy>( elementOutput )( Element::UPDATE_TMP_STATIC_DATA );

#ifdef _OPENMP
		if (mNumProcessors!=0)
			omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
	    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
		{
			try
			{
				Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
				if (error!=Error::SUCCESSFUL)
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] elements have returned multiple different error codes, can't handle that.");
				}
			}
			catch(NuTo::Exception& e)
			{
				std::stringstream ss;
				ss << ElementGetId(elementVector[countElement]);
				std::string exceptionStringLocal(e.ErrorMessage()
						+"[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
			catch(...)
			{
				std::stringstream ss;
				ss << ElementGetId(elementVector[countElement]);
				std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
		}
		if(exception>0)
		{
			throw MechanicsException(exceptionStringTotal);
		}
	}
	//std::cout << "NuTo::StructureBase::ElementTotalUpdateTmpStaticData " << mUpdateTmpStaticDataRequired << "\n";
	mUpdateTmpStaticDataRequired = false;
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return errorGlobal;
}

//! @brief saves static data of all elements
NuTo::Error::eError NuTo::StructureBase::ElementFatigueSaveStaticData()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
	Error::eError errorGlobal (Error::SUCCESSFUL);
    int exception(0);
    std::string exceptionStringTotal;

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputDummy>( elementOutput )( Element::FATIGUE_SAVE_STATIC_DATA );

#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
            Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
            if (error!=Error::SUCCESSFUL)
            {
            	if (errorGlobal==Error::SUCCESSFUL)
            		errorGlobal = error;
            	else if (errorGlobal!=error)
            		throw MechanicsException("[NuTo::StructureBase::ElementFatigueSaveStaticData] elements have returned multiple different error codes, can't handle that.");
            }
		}
		catch(NuTo::Exception& e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal(e.ErrorMessage()
					+"[NuTo::StructureBase::ElementFatigueSaveStaticData] Error saving static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal("[NuTo::StructureBase::ElementFatigueSaveStaticData] Error saving static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
	}
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueSaveStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueSaveStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return errorGlobal;
}

//! @brief restores static data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementFatigueRestoreStaticData()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
	Error::eError errorGlobal (Error::SUCCESSFUL);
    int exception(0);
    std::string exceptionStringTotal;

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputDummy>( elementOutput )( Element::FATIGUE_RESTORE_STATIC_DATA );

#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
            Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
            if (error!=Error::SUCCESSFUL)
            {
            	if (errorGlobal==Error::SUCCESSFUL)
            		errorGlobal = error;
            	else if (errorGlobal!=error)
            		throw MechanicsException("[NuTo::StructureBase::ElementFatigueRestoreStaticData] elements have returned multiple different error codes, can't handle that.");
            }
		}
		catch(NuTo::Exception& e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal(e.ErrorMessage()
					+"[NuTo::StructureBase::ElementFatigueRestoreStaticData] Error restoring static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal("[NuTo::StructureBase::ElementFatigueRestoreStaticData] Error saving static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
	}
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueRestoreStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueRestoreStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return errorGlobal;
}

//! @brief saves static data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementFatigueExtrapolateStaticData()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
	Error::eError errorGlobal (Error::SUCCESSFUL);
    int exception(0);
    std::string exceptionStringTotal;

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputDummy>( elementOutput )( Element::FATIGUE_EXTRAPOLATE_STATIC_DATA );

#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
            Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
            if (error!=Error::SUCCESSFUL)
            {
            	if (errorGlobal==Error::SUCCESSFUL)
            		errorGlobal = error;
            	else if (errorGlobal!=error)
            		throw MechanicsException("[NuTo::StructureBase::ElementFatigueExtrapolateStaticData] elements have returned multiple different error codes, can't handle that.");
            }
		}
		catch(NuTo::Exception& e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal(e.ErrorMessage()
					+"[NuTo::StructureBase::ElementFatigueExtrapolateStaticData] Error extrapolating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal("[NuTo::StructureBase::ElementFatigueExtrapolateStaticData] Error extrapolating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
	}
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueExtrapolateStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementFatigueExtrapolateStaticData] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return errorGlobal;
}

//! @brief calculates the average stress
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
//! @param rEngineeringStress  average stress (return value)
void NuTo::StructureBase::ElementTotalGetAverageStress(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStress;
    rEngineeringStress.Resize(6,1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
            elementVector[elementCount]->GetIntegratedStress(elementEngineeringStress);
            rEngineeringStress+=elementEngineeringStress;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStress] Error calculating integrated stress  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStress] Error calculating integrated stress  for element " + ss.str() + ".");
        }
    }
    rEngineeringStress*=1./rVolume;

#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalGetAverageStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalGetAverageStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
}

//! @brief calculates the average stress
//! @param rGroupId  group number
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
//! @param rEngineeringStress  average stress (return value)
void NuTo::StructureBase::ElementGroupGetAverageStress(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStress] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStress] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStress;
    rEngineeringStress.Resize(6,1);

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	itElement->second->GetIntegratedStress(elementEngineeringStress);
            rEngineeringStress+=elementEngineeringStress;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetAverageStress] Error calculating integrated stress  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetAverageStress] Error calculating integrated stress  for element " + ss.str() + ".");
        }
    }
    rEngineeringStress*=(1./rVolume);

#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementGroupGetAverageStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementGroupGetAverageStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
}


//! @brief calculates the average strain
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero strain, but still nonzero area)
//! @param rEngineeringStraiu  average strain (return value)
void NuTo::StructureBase::ElementTotalGetAverageStrain(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStrain;
    rEngineeringStrain.Resize(6,1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
            elementVector[elementCount]->GetIntegratedStrain(elementEngineeringStrain);
            rEngineeringStrain+=elementEngineeringStrain;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=1./rVolume;

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalGetAverageStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief calculates the average strain
//! @param rGroupId  group number
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero strain, but still nonzero area)
//! @param rEngineeringStrain  average strain (return value)
void NuTo::StructureBase::ElementGroupGetAverageStrain(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStrain] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStrain] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStrain;
    rEngineeringStrain.Resize(6,1);

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	itElement->second->GetIntegratedStrain(elementEngineeringStrain);
            rEngineeringStrain+=elementEngineeringStrain;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=(1./rVolume);


#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupGetAverageStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


void NuTo::StructureBase::ElementGroupGetMembers(int rGroupId, NuTo::FullVector<int,Eigen::Dynamic>& rMembers)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetMembers] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetMembers] Group is not an element group.");
    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);

    rMembers.Resize(elementGroup->GetNumMembers());
    int countElement(0);
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++,countElement++)
    {
       	rMembers[countElement] = itElement->first;
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupGetMembers] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


//! @brief calculates the total energy of the system
//! @return total energy
double NuTo::StructureBase::ElementTotalGetInternalEnergy()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    double totalEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    std::vector<const ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalGetInternalEnerg] not yet implemented on ip level.");
//            elementVector[elementCount]->GetIpData(NuTo::IpData::INTERNAL_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                totalEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetInternalEnergy] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetInternalEnergy] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger <<"[NuTo::StructureBase::ElementTotalGetTotalEnergy] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return totalEnergy;
}

//! @brief calculates the total energy of a group of elements
//! @return total energy
double NuTo::StructureBase::ElementGroupGetTotalEnergy(int rGroupId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Group is not an element group.");
    const Group<ElementBase> *elementGroup = dynamic_cast<const Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    double totalEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] not yet implemented on ip level.");
        	//itElement->second->GetIpData(NuTo::IpData::INTERNAL_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                totalEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupGetTotalEnergy] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return totalEnergy;
}

//! @brief calculates the elastic energy of the system
//! @return elastic energy
double NuTo::StructureBase::ElementTotalGetElasticEnergy()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    double elasticEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    std::vector<const ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalGetElasticEnergy] not yet implemented on ip level.");
        	//elementVector[elementCount]->GetIpData(NuTo::IpData::ELASTIC_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                elasticEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementTotalGetElasticEnergy] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return elasticEnergy;
}

//! @brief calculate the largest element eigenvalue for a group of elements solving the generalized eigenvalue problem Ku=lambda Mu
//! this is used for the estimation of the critical time step
double NuTo::StructureBase::ElementGroupCalculateLargestElementEigenvalue(int rGroupId)
{
    std::vector< NuTo::ElementBase*> elementVector;
    NuTo::GroupBase* grp_PtrBase = this->GroupGetGroupPtr(rGroupId);
    Group<NuTo::ElementBase> *grp_Ptr = grp_PtrBase->AsGroupElement();
    assert(grp_Ptr!=0);
	this->GetElementsByGroup(grp_Ptr,elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}

//! @brief calculate the critical time step for all elements solving the generalized eigenvalue problem Ku=lambda Mu
double NuTo::StructureBase::ElementTotalCalculateLargestElementEigenvalue()
{
    std::vector< ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}


//! @brief calculate the critical time step for a vector of elements solving the generalized eigenvalue problem Ku=lambda Mu
double NuTo::StructureBase::ElementCalculateLargestElementEigenvalue(const std::vector< ElementBase*>& rElementVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif

    Error::eError errorGlobal (Error::SUCCESSFUL);

    int exception(0);
    std::string exceptionStringTotal;

	double maxGlobalEigenValue(0);


#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
#endif //_OPENMP
	{

		boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

		boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>( elementOutput )( Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE );
		boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>( elementOutput )( Element::HESSIAN_0_TIME_DERIVATIVE );

		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;


#ifdef _OPENMP
    	#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
		for (unsigned int countElement=0;  countElement<rElementVector.size();countElement++)
		{
			try
			{
				Error::eError error = rElementVector[countElement]->Evaluate(elementOutput);
				if (error!=Error::SUCCESSFUL)
#ifdef _OPENMP
#pragma omp critical
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] multiple elements have returned different error codes, can't handle that.");
				}
#else //_OPENMP
				errorGlobal = error;
#endif //_OPENMP

				NuTo::FullVector<double,Eigen::Dynamic>&  lumpedMass(elementOutput.find(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->second->GetFullVectorDouble());
				NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  stiffness(elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second->GetFullMatrixDouble());

				//assuming the stiffness matrix is symmetric
				//std::cout << "lumped mass in element routine\n" <<  lumpedMass << std::endl;
				//std::cout << "element stiffness\n" << stiffness << std::endl;

				//eigenSolver.compute(stiffness);
				//std::cout << "eigenvalues element stiffness\n" << eigenSolver.eigenvalues() << std::endl;

				//std::cout << "eigenmatrix element\n" << lumpedMass.cwiseInverse().asDiagonal()*stiffness << std::endl;

				//invert the lumped mass matrix
				eigenSolver.compute(stiffness,lumpedMass.asDiagonal());
				//std::cout << "eigenvalues in element routine\n" <<  eigenSolver.eigenvalues() << std::endl;

				double maxElementEigenValue = eigenSolver.eigenvalues().maxCoeff();
				if (maxElementEigenValue>maxGlobalEigenValue)
				{
#ifdef _OPENMP
#pragma omp critical
#endif //_OPENMP
					maxGlobalEigenValue = maxElementEigenValue;
				}
			}
			catch(NuTo::Exception& e)
			{
				std::stringstream ss;
				ss << ElementGetId(rElementVector[countElement]);
				std::string exceptionStringLocal(e.ErrorMessage()
						+"[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] Error calculating critical time step for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
			catch(...)
			{
				std::stringstream ss;
				ss << ElementGetId(rElementVector[countElement]);
				std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] Error calculating critical time step for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
		}
	} //end of parallel region
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
    if (errorGlobal!=Error::SUCCESSFUL)
    {
	    throw MechanicsException("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] error calculating critical time step.");
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else //_OPENMP
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif //_OPENMP
#endif // SHOW_TIME

    return maxGlobalEigenValue;

}

//! @brief calculates the volume of the elements
//! @param rGroupId  group number
//! @return volume of the structure in 3D /area in 2D/ length in 1D
double NuTo::StructureBase::ElementGroupGetVolume(int rGroupId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetVolume] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetVolume] Group is not an element group.");
    const Group<ElementBase> *elementGroup = dynamic_cast<const Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    double totalVolume(0);

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	totalVolume += itElement->second->GetIntegrationPointVolume().sum();
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetVolume] Error calculating volume for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetVolume] Error calculating volume for element " + ss.str() + ".");
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupGetVolume] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
    return totalVolume;
}


#ifdef ENABLE_VISUALIZE
//! @brief ... adds all the elements in the vector to the data structure that is finally visualized
void NuTo::StructureBase::ElementTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    std::vector< ElementBase*> elementVec;
    this->GetElementsTotal(elementVec);
    ElementVectorAddToVisualize(rVisualize,rVisualizationList,elementVec);
}

//! @brief ... adds all the elements in a group to the data structure that is finally visualized
void NuTo::StructureBase::ElementGroupAddToVisualize(int rGroupId, VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
	// find group by name
	Group<ElementBase>* elementGroup = this->GroupGetGroupPtr(rGroupId)->AsGroupElement();
	std::vector< ElementBase*> elementVec;
	this->GetElementsByGroup(elementGroup,elementVec);
    ElementVectorAddToVisualize(rVisualize,rVisualizationList,elementVec,mGroupVisualizationType.at(rGroupId));

}

//! @brief ... adds all the elements in the vector to the data structure that is finally visualized
void NuTo::StructureBase::ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements, const VisualizeBase::eVisualizationType rVisualizationType)
{
    // build global tmp static data
    if (mHaveTmpStaticData and mUpdateTmpStaticDataRequired)
    	throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Update of tmpStaticData required first.");

    switch (rVisualizationType)
    {
        case VisualizeBase::VORONOI_CELL:
            for (auto const & iElePtr : rElements)
                iElePtr->Visualize(rVisualize, rVisualizationList);
            break;
        case VisualizeBase::EXTRAPOLATION_TO_NODES:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeExtrapolateToNodes(rVisualize, rVisualizationList);
            break;
        case VisualizeBase::POINTS:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeIntegrationPointData(rVisualize, rVisualizationList);
            break;
        default:
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization type not implemented.");
    }
}





#endif //VISUALIZE
