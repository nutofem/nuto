// $Id: ConstraintNodeGroupDisplacements2D.cpp 265 2010-06-08 08:47:00Z arnold2 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <math.h>
#include <iostream>


#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsMultiscale2D.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintLinearFineScaleDisplacementsPeriodic2D.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::GetNumLinearConstraints()const
{
    return 2*(mSlaveNodesRightBoundary.size()+mSlaveNodesTopBoundary.size())+2;
}

 //! @brief constructor
NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::ConstraintLinearFineScaleDisplacementsPeriodic2D(const Group<NodeBase>* rGroupBoundary,const EngineeringStrain2D& rStrain) :  ConstraintLinear()
{
	mGroupBoundaryNodes = rGroupBoundary;
	mStrain = rStrain;
    SetBoundaryVectors();
}

//!@brief set the strain of the periodic boundary conditions
//!@param rStrain strain (e_xx,e_yy,gamma_xy)
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetStrain(const EngineeringStrain2D& rStrain)
{
    mStrain = rStrain;
}

#define tolerance 1e-8
//!@brief calculate the border vectors in counterclockwise direction
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors()
{
    //determine min and max coordinates
	double minX, maxX, minY, maxY;
	double coordinates[2];
	if (mGroupBoundaryNodes->GetNumMembers()<4)
	{
		throw MechanicsException("[NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors] number of boundary nodes is less than 4, check your code.");
	}

	Group<NodeBase>::const_iterator itNode=mGroupBoundaryNodes->begin();
	(itNode->second)->GetCoordinates2D(coordinates);
	minX =  coordinates[0];
	maxX =  coordinates[0];
	minY =  coordinates[1];
	maxY =  coordinates[1];
	for (itNode=mGroupBoundaryNodes->begin()++; itNode!=mGroupBoundaryNodes->end();itNode++)
    {
    	(itNode->second)->GetCoordinates2D(coordinates);
        if (minX>coordinates[0])
        	minX=coordinates[0];
        if (maxX<coordinates[0])
        	maxX=coordinates[0];
        if (minY>coordinates[1])
        	minY=coordinates[1];
        if (maxY<coordinates[1])
        	maxY=coordinates[1];
    }
    //calculate master nodes left boundary (green)
    mMasterNodesLeftBoundary.resize(0);
    mSlaveNodesRightBoundary.resize(0);
    mMasterNodesBottomBoundary.resize(0);
    mSlaveNodesTopBoundary.resize(0);
    for (Group<NodeBase>::const_iterator itNode=mGroupBoundaryNodes->begin(); itNode!=mGroupBoundaryNodes->end();itNode++)
    {
    	(itNode->second)->GetCoordinates2D(coordinates);
    	if (fabs(coordinates[0]-minX)<tolerance)
            mMasterNodesLeftBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacementsMultiscale2D*>(itNode->second));
    	if (fabs(coordinates[0]-maxX)<tolerance)
    		mSlaveNodesRightBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacementsMultiscale2D*>(itNode->second));
    	if (fabs(coordinates[1]-minY)<tolerance)
            mMasterNodesBottomBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacementsMultiscale2D*>(itNode->second));
    	if (fabs(coordinates[1]-maxY)<tolerance && fabs(coordinates[0]-minX)>tolerance)
    		mSlaveNodesTopBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacementsMultiscale2D*>(itNode->second));
    }

    if (mMasterNodesLeftBoundary.size()<2)
    	throw MechanicsException("[NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors] left boundary has less than 2 nodes.");
    if (mSlaveNodesRightBoundary.size()<1)
    	throw MechanicsException("[NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors] right boundary has less than 1 nodes.");
    if (mMasterNodesBottomBoundary.size()<2)
    	throw MechanicsException("[NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors] bottom boundary has less than 2 nodes.");
    if (mSlaveNodesTopBoundary.size()<1)
    	throw MechanicsException("[NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::SetBoundaryVectors] top boundary has less than 1 nodes.");

    sort(mMasterNodesLeftBoundary.begin(), mMasterNodesLeftBoundary.end(), less_YCoordinate2D());
    sort(mSlaveNodesRightBoundary.begin(), mSlaveNodesRightBoundary.end(), less_YCoordinate2D());
    sort(mMasterNodesBottomBoundary.begin(), mMasterNodesBottomBoundary.end(), less_XCoordinate2D());
    sort(mSlaveNodesTopBoundary.begin(), mSlaveNodesTopBoundary.end(), less_XCoordinate2D());

    /*
     //Info about the nodes
	for (int countBoundary=0; countBoundary<4; countBoundary++)
	{
		std::vector<NodeCoordinatesDisplacementsMultiscale2D*>* nodeVectorPtr;
		switch(countBoundary)
		{
		case 0:
			nodeVectorPtr = &mMasterNodesLeftBoundary;
			std::cout << "mMasterNodesLeftBoundary " << std::endl;
			break;
		case 1:
			nodeVectorPtr = &mSlaveNodesRightBoundary;
			std::cout << "mSlaveNodesRightBoundary " << std::endl;
			break;
		case 2:
			nodeVectorPtr = &mMasterNodesBottomBoundary;
			std::cout << "mMasterNodesBottomBoundary " << std::endl;
			break;
		default:
			nodeVectorPtr = &mSlaveNodesTopBoundary;
			std::cout << "mSlaveNodesTopBoundary " << std::endl;
			break;
		}
		for (unsigned int countNodes=0; countNodes<(*nodeVectorPtr).size(); countNodes++)
		{
			double coordinates[2];
			(*nodeVectorPtr)[countNodes]->GetCoordinates2D(coordinates);
			std::cout << "  " << coordinates[0] << " " << coordinates[1] << std::endl;
		}
    }
*/
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
#define MIN_CONSTRAINT 1e-6
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const
{
    //only du linear interpolation between nodes on the boundary, this is not exact for quadrativ elements, but the effort is much reduced
    //**************************************************************
	//add  constraints for all the slave nodes of the right boundary
	//**************************************************************
	assert(mMasterNodesLeftBoundary.size()>1);
	unsigned int nextMasterNodecount(1);

	NodeCoordinatesDisplacementsMultiscale2D* curMasterNodePtr(mMasterNodesLeftBoundary[0]);
	double coordinatesCurMaster[2];
	curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

	NodeCoordinatesDisplacementsMultiscale2D* nextMasterNodePtr(mMasterNodesLeftBoundary[nextMasterNodecount]);;
	double coordinatesNextMaster[2];
	nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

	//double deltaDisp[2];
	for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
	{
		NodeCoordinatesDisplacementsMultiscale2D* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
		double coordinatesSlave[2];
		curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

		while (coordinatesNextMaster[1]<coordinatesSlave[1] && nextMasterNodecount+1<mMasterNodesLeftBoundary.size())
		{
			 curMasterNodePtr = nextMasterNodePtr;
			 coordinatesCurMaster[0] = coordinatesNextMaster[0];
			 coordinatesCurMaster[1] = coordinatesNextMaster[1];
			 nextMasterNodecount++;

			 assert (nextMasterNodecount<mMasterNodesLeftBoundary.size());
			 nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
			 nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
		}
		//slave is between two master nodes

		//calculate weighting function for each master node
		double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlave[1]);

		//deltaDisp[0] = (coordinatesCurMaster[0]-coordinatesSlave[0])*mStrain.mEngineeringStrain[0];
		//deltaDisp[1] = (coordinatesCurMaster[0]-coordinatesSlave[0])*0.5*mStrain.mEngineeringStrain[2];

		//constrain x direction
		rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(0),1);
		if (fabs(w)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(0),-w);
		if (fabs(w-1.)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(0),w-1.);
		//rRHS(curConstraintEquation,0) = deltaDisp[0];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(0) << "=" << deltaDisp[0] << "\n";
		curConstraintEquation++;

		//constrain y direction
		rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(1),1);
		if (fabs(w)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(1),-w);
		if (fabs(w-1.)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(1),w-1.);
		//rRHS(curConstraintEquation,0) = deltaDisp[1];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(1) << "=" << deltaDisp[1] << "\n";
		curConstraintEquation++;
	}

	//**************************************************************
	// add  constraints for all the slave nodes of the top boundary
	//**************************************************************
	assert(mMasterNodesBottomBoundary.size()>1);
	nextMasterNodecount = 1;

	curMasterNodePtr = mMasterNodesBottomBoundary[0];
	curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

	nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
	nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

	for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
	{
		NodeCoordinatesDisplacementsMultiscale2D* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
		double coordinatesSlave[2];
		curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

		while (coordinatesNextMaster[0]<coordinatesSlave[0] && nextMasterNodecount+1<mMasterNodesBottomBoundary.size())
		{
			 curMasterNodePtr = nextMasterNodePtr;
			 coordinatesCurMaster[0] = coordinatesNextMaster[0];
			 coordinatesCurMaster[1] = coordinatesNextMaster[1];
			 nextMasterNodecount++;

			 assert (nextMasterNodecount<mMasterNodesBottomBoundary.size());
			 nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
			 nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
		}
		//slave is between two master nodes

		//calculate weighting function for each master node
		double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlave[0]);

		//deltaDisp[0] = (coordinatesCurMaster[1]-coordinatesSlave[1])*0.5*mStrain.mEngineeringStrain[2] ;
		//deltaDisp[1] = (coordinatesCurMaster[1]-coordinatesSlave[1])*mStrain.mEngineeringStrain[1] ;

		//constrain x direction
		rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(0),1);
		if (fabs(w)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(0),-w);
		if (fabs(w-1.)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(0),w-1.);
		//rRHS(curConstraintEquation,0) = deltaDisp[0];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(0) << "=" << deltaDisp[0] << "\n";
		curConstraintEquation++;

		//constrain y direction
		rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(1),1);
		if (fabs(w)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(1),-w);
		if (fabs(w-1.)>MIN_CONSTRAINT)
			rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(1),w-1.);
		//rRHS(curConstraintEquation,0) = deltaDisp[1];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(1) << "=" << deltaDisp[1] << "\n";
		curConstraintEquation++;
	}

	//add constraint in order to avoid rigid body translations
	//instead of fixing a single node, I just fix the sum of x(lower left)+x(lower right)=0 and y(lower left)+x(upper left)=0
	//constrain x direction
	rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesBottomBoundary[0]->GetDofFineScaleDisplacement(0),1);
	rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesBottomBoundary[mMasterNodesBottomBoundary.size()-1]->GetDofFineScaleDisplacement(0),1);
	//std::cout << "constraint " << curConstraintEquation << ":" << mMasterNodesBottomBoundary[0]->GetDofFineScaleDisplacement(0) << "+" << mMasterNodesBottomBoundary[mMasterNodesBottomBoundary.size()-1]->GetDofFineScaleDisplacement(0) << "=" << 0 << "\n";
	//rRHS(curConstraintEquation,0) = 0;
	curConstraintEquation++;

	//constrain y direction
	rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesLeftBoundary[0]->GetDofFineScaleDisplacement(1),1);
	rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesLeftBoundary[mMasterNodesLeftBoundary.size()-1]->GetDofFineScaleDisplacement(1),1);
	//rRHS(curConstraintEquation,0) = 0;
	//std::cout << "constraint " << curConstraintEquation << ":" << mMasterNodesLeftBoundary[0]->GetDofFineScaleDisplacement(1) << "+" << mMasterNodesLeftBoundary[mMasterNodesLeftBoundary.size()-1]->GetDofFineScaleDisplacement(1) << "=" << 0 << "\n";
	curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::GetRHS(int& curConstraintEquation,NuTo::FullMatrix<double>& rRHS)const
{
    //only du linear interpolation between nodes on the boundary, this is not exact for quadrativ elements, but the effort is much reduced
    //**************************************************************
	//add  constraints for all the slave nodes of the right boundary
	//**************************************************************
	assert(mMasterNodesLeftBoundary.size()>1);
	unsigned int nextMasterNodecount(1);

	NodeCoordinatesDisplacementsMultiscale2D* curMasterNodePtr(mMasterNodesLeftBoundary[0]);
	double coordinatesCurMaster[2];
	curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

	NodeCoordinatesDisplacementsMultiscale2D* nextMasterNodePtr(mMasterNodesLeftBoundary[nextMasterNodecount]);;
	double coordinatesNextMaster[2];
	nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

	double deltaDisp[2];
	for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
	{
		NodeCoordinatesDisplacementsMultiscale2D* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
		double coordinatesSlave[2];
		curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

		while (coordinatesNextMaster[1]<coordinatesSlave[1] && nextMasterNodecount+1<mMasterNodesLeftBoundary.size())
		{
			 curMasterNodePtr = nextMasterNodePtr;
			 coordinatesCurMaster[0] = coordinatesNextMaster[0];
			 coordinatesCurMaster[1] = coordinatesNextMaster[1];
			 nextMasterNodecount++;

			 assert (nextMasterNodecount<mMasterNodesLeftBoundary.size());
			 nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
			 nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
		}
		//slave is between two master nodes

		//calculate weighting function for each master node
		//double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlave[1]);

		deltaDisp[0] = (coordinatesCurMaster[0]-coordinatesSlave[0])*mStrain.mEngineeringStrain[0];
		deltaDisp[1] = (coordinatesCurMaster[0]-coordinatesSlave[0])*0.5*mStrain.mEngineeringStrain[2];

		//constrain x direction
		//rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(0),1);
		//if (fabs(w)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(0),-w);
		//if (fabs(w-1.)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(0),w-1.);
		//rRHS(curConstraintEquation,0) = deltaDisp[0];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(0) << "=" << deltaDisp[0] << "\n";
		curConstraintEquation++;

		//constrain y direction
		//rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(1),1);
		//if (fabs(w)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(1),-w);
		//if (fabs(w-1.)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(1),w-1.);
		rRHS(curConstraintEquation,0) = deltaDisp[1];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(1) << "=" << deltaDisp[1] << "\n";
		curConstraintEquation++;
	}

	//**************************************************************
	// add  constraints for all the slave nodes of the top boundary
	//**************************************************************
	assert(mMasterNodesBottomBoundary.size()>1);
	nextMasterNodecount = 1;

	curMasterNodePtr = mMasterNodesBottomBoundary[0];
	curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

	nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
	nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

	for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
	{
		NodeCoordinatesDisplacementsMultiscale2D* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
		double coordinatesSlave[2];
		curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

		while (coordinatesNextMaster[0]<coordinatesSlave[0] && nextMasterNodecount+1<mMasterNodesBottomBoundary.size())
		{
			 curMasterNodePtr = nextMasterNodePtr;
			 coordinatesCurMaster[0] = coordinatesNextMaster[0];
			 coordinatesCurMaster[1] = coordinatesNextMaster[1];
			 nextMasterNodecount++;

			 assert (nextMasterNodecount<mMasterNodesBottomBoundary.size());
			 nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
			 nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
		}
		//slave is between two master nodes

		//calculate weighting function for each master node
		//double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlave[0]);

		deltaDisp[0] = (coordinatesCurMaster[1]-coordinatesSlave[1])*0.5*mStrain.mEngineeringStrain[2] ;
		deltaDisp[1] = (coordinatesCurMaster[1]-coordinatesSlave[1])*mStrain.mEngineeringStrain[1] ;

		//constrain x direction
		//rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(0),1);
		//if (fabs(w)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(0),-w);
		//if (fabs(w-1.)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(0),w-1.);
		rRHS(curConstraintEquation,0) = deltaDisp[0];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(0) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(0) << "=" << deltaDisp[0] << "\n";
		curConstraintEquation++;

		//constrain y direction
		//rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofFineScaleDisplacement(1),1);
		//if (fabs(w)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofFineScaleDisplacement(1),-w);
		//if (fabs(w-1.)>MIN_CONSTRAINT)
		//	rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofFineScaleDisplacement(1),w-1.);
		rRHS(curConstraintEquation,0) = deltaDisp[1];
		//std::cout << "constraint " << curConstraintEquation << ":" << -w << "*" << curMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << w-1. << "*" << nextMasterNodePtr->GetDofFineScaleDisplacement(1) << "+" << curSlaveNodePtr->GetDofFineScaleDisplacement(1) << "=" << deltaDisp[1] << "\n";
		curConstraintEquation++;
	}

	//add constraint in order to avoid rigid body translations
	//instead of fixing a single node, I just fix the sum of x(lower left)+x(lower right)=0 and y(lower left)+x(upper left)=0
	//constrain x direction
	//rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesBottomBoundary[0]->GetDofFineScaleDisplacement(0),1);
	//rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesBottomBoundary[mMasterNodesBottomBoundary.size()-1]->GetDofFineScaleDisplacement(0),1);
	//std::cout << "constraint " << curConstraintEquation << ":" << mMasterNodesBottomBoundary[0]->GetDofFineScaleDisplacement(0) << "+" << mMasterNodesBottomBoundary[mMasterNodesBottomBoundary.size()-1]->GetDofFineScaleDisplacement(0) << "=" << 0 << "\n";
	rRHS(curConstraintEquation,0) = 0;
	curConstraintEquation++;

	//constrain y direction
	//rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesLeftBoundary[0]->GetDofFineScaleDisplacement(1),1);
	//rConstraintMatrix.AddEntry(curConstraintEquation,mMasterNodesLeftBoundary[mMasterNodesLeftBoundary.size()-1]->GetDofFineScaleDisplacement(1),1);
	rRHS(curConstraintEquation,0) = 0;
	//std::cout << "constraint " << curConstraintEquation << ":" << mMasterNodesLeftBoundary[0]->GetDofFineScaleDisplacement(1) << "+" << mMasterNodesLeftBoundary[mMasterNodesLeftBoundary.size()-1]->GetDofFineScaleDisplacement(1) << "=" << 0 << "\n";
	curConstraintEquation++;

}

//calculate weighting function for each master node
double NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::CalculateWeightFunction(double rCoordinateCurMaster, double rCoordinateNextMaster, double rCoordinateSlave)const
{
    assert(rCoordinateNextMaster!=rCoordinateCurMaster);
    return 1.-(rCoordinateSlave-rCoordinateCurMaster)/(rCoordinateNextMaster-rCoordinateCurMaster);
}

//calculate delta displacement in x and y direction from the applied strain and the nodal position
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::CalculateDeltaDisp(double rCoordinates[2], double rDeltaDisp[2])const
{
    rDeltaDisp[0] = mStrain.mEngineeringStrain[0]*rCoordinates[0] + mStrain.mEngineeringStrain[2]*rCoordinates[1];
    rDeltaDisp[1] = mStrain.mEngineeringStrain[1]*rCoordinates[1] + mStrain.mEngineeringStrain[2]*rCoordinates[0];
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::Info(unsigned short rVerboseLevel) const
{
	    std::cout << "NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D, strain. " <<  mStrain.mEngineeringStrain[0] << " " <<  mStrain.mEngineeringStrain[1] << " "<<  mStrain.mEngineeringStrain[2] <<std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearFineScaleDisplacementsPeriodic2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mSlaveNodesRightBoundary)
       & BOOST_SERIALIZATION_NVP(mSlaveNodesTopBoundary)
       & BOOST_SERIALIZATION_NVP(mMasterNodesLeftBoundary)
       & BOOST_SERIALIZATION_NVP(mMasterNodesBottomBoundary)
       & BOOST_SERIALIZATION_NVP(mStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearFineScaleDisplacementsPeriodic2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D)
#endif // ENABLE_SERIALIZATION
