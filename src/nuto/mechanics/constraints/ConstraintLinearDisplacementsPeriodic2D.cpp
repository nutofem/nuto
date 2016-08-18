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
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearDisplacementsPeriodic2D::GetNumLinearConstraints()const
{
    return 2*(mSlaveNodesRightBoundary.size()+mSlaveNodesTopBoundary.size());
}

 //! @brief constructor
NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D(const StructureBase* rStructure, double rAngle,
        const EngineeringStrain<2>& rStrain,NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rCrackOpening, double rRadiusToCrackWithoutConstraints,
        const Group<NodeBase>* rGroupTop,const Group<NodeBase>* rGroupBottom,
        const Group<NodeBase>* rGroupLeft, const Group<NodeBase>* rGroupRight) :  ConstraintLinear()
{
    mStructure = rStructure,

    mRadiusToCrackWithoutConstraints = rRadiusToCrackWithoutConstraints;

    SetCrackOpening(rCrackOpening);

    mStrain = rStrain;

    mGroupTop = rGroupTop;
    mGroupBottom = rGroupBottom;
    mGroupLeft = rGroupLeft;
    mGroupRight = rGroupRight;

    //determine the corner nodes
    Group<NodeBase>* newGroup = dynamic_cast<Group<NodeBase>*>(rGroupTop->Intersection (rGroupLeft));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Top left corner node can not be determined.");
    mLeftUpperCorner = newGroup->begin()->second;
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupTop->Intersection (rGroupRight));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Top right corner node can not be determined.");
    mRightUpperCorner = newGroup->begin()->second;
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupBottom->Intersection (rGroupLeft));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Bottom left corner node can not be determined.");
    mLeftLowerCorner = newGroup->begin()->second;
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupBottom->Intersection (rGroupRight));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Bottom right corner node can not be determined.");
    mRightLowerCorner = newGroup->begin()->second;
    delete newGroup;

    //get box coordinates
    if (mLeftUpperCorner->GetNum(Node::COORDINATES)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Upper left node does not have 2 coordinates.");
    if (mRightUpperCorner->GetNum(Node::COORDINATES)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Upper right node does not have 2 coordinates.");
    if (mLeftLowerCorner->GetNum(Node::COORDINATES)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Lower left node does not have 2 coordinates.");
    if (mRightLowerCorner->GetNum(Node::COORDINATES)!=2)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Lower right node does not have 2 coordinates.");


    Eigen::Matrix<double, 2, 1> LeftUpperCoordinates  = mLeftUpperCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> RightUpperCoordinates = mRightUpperCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> LeftLowerCoordinates  = mLeftLowerCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> RightLowerCoordinates = mRightLowerCorner->Get(Node::COORDINATES);

    if (mStructure->GetVerboseLevel()>0)
    {
        std::cout << "Left upper corner is node " << mStructure->NodeGetId(mLeftUpperCorner) << std::endl;
        std::cout << "Left lower corner is node " << mStructure->NodeGetId(mLeftLowerCorner) << std::endl;
        std::cout << "Right upper corner is node " << mStructure->NodeGetId(mRightUpperCorner) << std::endl;
        std::cout << "Right lower corner is node " << mStructure->NodeGetId(mRightLowerCorner) << std::endl;
    }

    //check box
    if (LeftUpperCoordinates[0]!=LeftLowerCoordinates[0])
    {
        std::cout << "LeftUpperCoordinates " << LeftUpperCoordinates[0] << " " << LeftUpperCoordinates[1] << std::endl;
        std::cout << "LeftLowerCoordinates " << LeftLowerCoordinates[0] << " " << LeftLowerCoordinates[1] << std::endl;
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Left boundary is not correct - check your node groups.");
    }
    if (RightUpperCoordinates[0]!=RightLowerCoordinates[0])
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] Right boundary is not correct - check your node groups.");

    if (RightUpperCoordinates[0]<=LeftUpperCoordinates[0])
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] left boundary coordinate is larger than right boundary coordinate.");

    if (LeftUpperCoordinates[1]!=RightUpperCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] upper boundary is not correct - check your node groups.");

    if (LeftLowerCoordinates[0]!=RightLowerCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] lower boundary is not correct - check your node groups.");

    if (LeftUpperCoordinates[1]<=LeftLowerCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::ConstraintLinearDisplacementsPeriodic2D] upper boundary coordinate is larger than lower boundary coordinate.");


    SetAngle(rAngle);
}

//!@brief set the angle of the periodic boundary conditions
//!@param rRHS new right hand side
void NuTo::ConstraintLinearDisplacementsPeriodic2D::SetAngle(double rAngle)
{
    mAngle = rAngle;
    while (mAngle>225 || mAngle<45)
    {
        if (mAngle>225)
        {
            mAngle-=180;
            mCrackOpening[0] *= -1.;
            mCrackOpening[1] *= -1.;
        }
        if (mAngle<45)
        {
            mAngle+=180;
            mCrackOpening[0] *= -1.;
            mCrackOpening[1] *= -1.;
        }
    }
   SetBoundaryVectors();
}

//!@brief sets/modifies the average strain applied to the boundary
//!@param rAngle angle in deg
void NuTo::ConstraintLinearDisplacementsPeriodic2D::SetCrackOpening(const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCrackOpening)
{
    if (rCrackOpening.GetNumRows()!=2 || rCrackOpening.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D::SetCrackOpening] crack opening should be a (2,1) matrix.");
    mCrackOpening[0] = rCrackOpening(0,0);
    mCrackOpening[1] = rCrackOpening(1,0);
}

//!@brief set the strain of the periodic boundary conditions
//!@param rStrain strain (e_xx,e_yy,gamma_xy)
void NuTo::ConstraintLinearDisplacementsPeriodic2D::SetStrain(const EngineeringStrain<2>& rStrain)
{
    mStrain = rStrain;
}

#define PI 3.14159265359

//!@brief calculate the border vectors in counterclockwise direction
void NuTo::ConstraintLinearDisplacementsPeriodic2D::SetBoundaryVectors()
{
    Eigen::Matrix<double, 2, 1> LeftUpperCoordinates = mLeftUpperCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> LeftLowerCoordinates = mLeftLowerCorner->Get(Node::COORDINATES);

    //calculate length of specimen
    double length = LeftUpperCoordinates[1]-LeftLowerCoordinates[1];
    //calculate master nodes left boundary (green)
    mMasterNodesLeftBoundary.resize(0);
    mMasterNodesLeftBoundary.reserve(mGroupLeft->GetNumMembers());
    for (Group<NodeBase>::const_iterator itNode=mGroupLeft->begin(); itNode!=mGroupLeft->end();itNode++)
    {
        mMasterNodesLeftBoundary.push_back(itNode->second);
    }
    sort(mMasterNodesLeftBoundary.begin(), mMasterNodesLeftBoundary.end(), less_YCoordinate2D());

    //calculate master nodes bottom boundary (yellow-orange-blue)
    mMasterNodesBottomBoundary.resize(0);
    mMasterNodesBottomBoundary.reserve(mGroupBottom->GetNumMembers());
    for (Group<NodeBase>::const_iterator itNode=mGroupBottom->begin(); itNode!=mGroupBottom->end();itNode++)
    {
        mMasterNodesBottomBoundary.push_back(itNode->second);
    }
    sort(mMasterNodesBottomBoundary.begin(), mMasterNodesBottomBoundary.end(), less_XCoordinate2D());

    //calculate slave nodes right boundary (green)
    mSlaveNodesRightBoundary.resize(0);
    if (mAngle>=45 && mAngle<135)
    {
        double crackShift;
        crackShift = length*tan((90-mAngle)*PI/180.);
        mSlaveNodesRightBoundary.reserve(mGroupRight->GetNumMembers());
        for (Group<NodeBase>::const_iterator itNode=mGroupRight->begin(); itNode!=mGroupRight->end();itNode++)
        {
            double coordinate = itNode->second->Get(Node::COORDINATES)[1];
            double DeltaX((length-crackShift)*0.5);
            double DeltaY(length-coordinate);

            if (DeltaX*DeltaX + DeltaY*DeltaY >=mRadiusToCrackWithoutConstraints * mRadiusToCrackWithoutConstraints)
                mSlaveNodesRightBoundary.push_back(itNode->second);
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        mSlaveNodesRightBoundary.reserve(mGroupRight->GetNumMembers()-1);

        for (Group<NodeBase>::const_iterator itNode=mGroupRight->begin(); itNode!=mGroupRight->end();itNode++)
        {
            if (itNode->second!=mRightLowerCorner)
            {
                double coordinate = itNode->second->Get(Node::COORDINATES)[1];
                if (std::abs(coordinate-(length+crackShift)*0.5)>=mRadiusToCrackWithoutConstraints)
                    mSlaveNodesRightBoundary.push_back(itNode->second);
            }
        }
    }

    sort(mSlaveNodesRightBoundary.begin(), mSlaveNodesRightBoundary.end(), less_YCoordinate2D());

    //calculate slave nodes top boundary
    mSlaveNodesTopBoundary.resize(0);
    if (mAngle>=45 && mAngle<135)
    {
        double crackShift;
        crackShift = length*tan((90-mAngle)*PI/180.);
        mSlaveNodesTopBoundary.reserve(mGroupTop->GetNumMembers()-1);
        for (Group<NodeBase>::const_iterator itNode=mGroupTop->begin(); itNode!=mGroupTop->end();itNode++)
        {
            if (itNode->second!=mLeftUpperCorner)
            {
                NodeBase* nodePtr = itNode->second;
                double coordinate = nodePtr->Get(Node::COORDINATES)[0];
                if (std::abs(coordinate-(length+crackShift)*0.5)>=mRadiusToCrackWithoutConstraints)
                    mSlaveNodesTopBoundary.push_back(nodePtr);
            }
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        mSlaveNodesTopBoundary.reserve(mGroupTop->GetNumMembers());
        for (Group<NodeBase>::const_iterator itNode=mGroupTop->begin(); itNode!=mGroupTop->end();itNode++)
        {
            NodeBase* nodePtr = itNode->second;
            double coordinate = nodePtr->Get(Node::COORDINATES)[0];
            double DeltaX(coordinate-(length+crackShift)*0.5);
            double DeltaY((length-crackShift)*0.5);

            if (DeltaX*DeltaX + DeltaY*DeltaY >=mRadiusToCrackWithoutConstraints * mRadiusToCrackWithoutConstraints)
                mSlaveNodesTopBoundary.push_back(itNode->second);
        }

    }
    sort(mSlaveNodesTopBoundary.begin(), mSlaveNodesTopBoundary.end(), less_XCoordinate2D());


    //Info about the nodes
    if (mStructure->GetVerboseLevel()>0)
    {
        for (int countBoundary=0; countBoundary<4; countBoundary++)
        {
            std::vector<NodeBase*>* nodeVectorPtr;
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
                Eigen::Matrix<double, 2, 1> coordinates = (*nodeVectorPtr)[countNodes]->Get(Node::COORDINATES);
                std::cout << "  " << mStructure->NodeGetId((*nodeVectorPtr)[countNodes]) << ": " << coordinates[0] << " " << coordinates[1] << std::endl;
            }
        }
    }
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
#define MIN_CONSTRAINT 1e-6
void NuTo::ConstraintLinearDisplacementsPeriodic2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrix<double>& rConstraintMatrix)const
{
    Eigen::Matrix<double, 2, 1> LeftUpperCoordinates = mLeftUpperCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> LeftLowerCoordinates = mLeftLowerCorner->Get(Node::COORDINATES);

    //calculate length of specimen
    double length = LeftUpperCoordinates[1]-LeftLowerCoordinates[1];

    //only du linear interpolation between nodes on the boundary, this is not exact for quadrativ elements, but the effort is much reduced
    if (mAngle>=45 && mAngle<135)
    {
        double crackShift;
        crackShift = length*tan((90-mAngle)*PI/180.);
       //**************************************************************
        //add  constraints for all the slave nodes of the right boundary
        //**************************************************************
        assert(mMasterNodesLeftBoundary.size()>1);
        unsigned int nextMasterNodecount(1);

        NodeBase* curMasterNodePtr(mMasterNodesLeftBoundary[0]);
        Eigen::Matrix<double, 2, 1> coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        NodeBase* nextMasterNodePtr(mMasterNodesLeftBoundary[nextMasterNodecount]);;
        Eigen::Matrix<double, 2, 1> coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        //double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
        	Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideY = coordinatesSlave[1];

            while (coordinatesNextMaster[1]<coordinatesSlaveonMasterSideY && nextMasterNodecount+1<mMasterNodesLeftBoundary.size())
            {
                 curMasterNodePtr = nextMasterNodePtr;
                 coordinatesCurMaster = coordinatesNextMaster;
                 nextMasterNodecount++;

                 assert (nextMasterNodecount<mMasterNodesLeftBoundary.size());
                 nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
                 coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }
            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlaveonMasterSideY);

            //deltaDisp[0] = length*mStrain[0] + mCrackOpening[0];
            //deltaDisp[1] = length*0.5*mStrain[2] + mCrackOpening[1];

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << 1-w << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 0),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 1),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the top boundary
        //**************************************************************

        assert(mMasterNodesBottomBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesBottomBoundary[0];
        coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
        coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double crackPosX((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
            Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideX = coordinatesSlave[0]-crackShift;
            double deltaRHS[2];
            if (coordinatesSlaveonMasterSideX<0)
            {
                coordinatesSlaveonMasterSideX+=length;
                deltaRHS[0] = mStrain[0]*length+mCrackOpening[0];
                deltaRHS[1] = 0.5 * mStrain[2]*length+mCrackOpening[1];
            }
            else
            {
                if (coordinatesSlaveonMasterSideX>length)
                {
                    coordinatesSlaveonMasterSideX-=length;
                    deltaRHS[0] = -mStrain[0]*length-mCrackOpening[0];
                    deltaRHS[1] = -0.5 * mStrain[2]*length-mCrackOpening[1];
                }
                else
                {
                    deltaRHS[0] = 0.;
                    deltaRHS[1] = 0.;
                }
            }

            while (coordinatesCurMaster[0]>coordinatesSlaveonMasterSideX || coordinatesNextMaster[0]<coordinatesSlaveonMasterSideX)
            {
                if (nextMasterNodePtr==mRightLowerCorner)
                {
                    //restart from the left corner
                    curMasterNodePtr = mMasterNodesBottomBoundary[0];
                    coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

                    nextMasterNodecount=1;
                }
                else
                {
                    curMasterNodePtr = nextMasterNodePtr;
                    coordinatesCurMaster[0] = coordinatesNextMaster[0];
                    coordinatesCurMaster[1] = coordinatesNextMaster[1];

                    nextMasterNodecount++;
                }
                nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
                coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }

            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlaveonMasterSideX);

            if (coordinatesCurMaster[0]>crackPosX)
            {
                deltaRHS[0]-=w*mCrackOpening[0];
                deltaRHS[1]-=w*mCrackOpening[1];
            }
            if (coordinatesNextMaster[0]>crackPosX)
            {
                deltaRHS[0]-=(1-w)*mCrackOpening[0];
                deltaRHS[1]-=(1-w)*mCrackOpening[1];
            }
            if (coordinatesSlaveonMasterSideX>crackPosX)
            {
                deltaRHS[0]+=mCrackOpening[0];
                deltaRHS[1]+=mCrackOpening[1];
            }

            //deltaDisp[0] = crackShift*mStrain[0] + length *0.5* mStrain[2] - deltaRHS[0];
            //deltaDisp[1] = crackShift*0.5*mStrain[2] + length*mStrain[1] - deltaRHS[1];

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 0),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 1),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        //**************************************************************
        // angle between 135 abnd 225 degrees
        //**************************************************************
        // add  constraints for all the slave nodes of the top boundary
        //**************************************************************
        assert(mMasterNodesBottomBoundary.size()>1);
        unsigned int nextMasterNodecount(1);

        NodeBase* curMasterNodePtr(mMasterNodesBottomBoundary[0]);
        Eigen::Matrix<double, 2, 1> coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        NodeBase* nextMasterNodePtr(mMasterNodesBottomBoundary[nextMasterNodecount]);;
        Eigen::Matrix<double, 2, 1> coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        //double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
        	Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideX;
            coordinatesSlaveonMasterSideX = coordinatesSlave[0];

            while (coordinatesNextMaster[0]<coordinatesSlaveonMasterSideX && nextMasterNodecount+1<mMasterNodesBottomBoundary.size())
            {
                 curMasterNodePtr = nextMasterNodePtr;
                 coordinatesCurMaster = coordinatesNextMaster;
                 nextMasterNodecount++;

                 assert (nextMasterNodecount<mMasterNodesBottomBoundary.size());
                 nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
                 coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }
            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlaveonMasterSideX);

            //deltaDisp[0] = length*0.5*mStrain[2] + mCrackOpening[0];;
            //deltaDisp[1] = length*mStrain[1] + mCrackOpening[1];;

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 0),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 1),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the right boundary
        //**************************************************************

        assert(mMasterNodesLeftBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesLeftBoundary[0];
        coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
        coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double crackPosY((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
            Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideY;
            coordinatesSlaveonMasterSideY = coordinatesSlave[1]-crackShift;
            double deltaRHS[2];
            if (coordinatesSlaveonMasterSideY<0)
            {
                coordinatesSlaveonMasterSideY+=length;
                deltaRHS[0] = 0.5 * mStrain[2]*length+mCrackOpening[0];
                deltaRHS[1] = mStrain[1]*length+mCrackOpening[1];
            }
            else
            {
                if (coordinatesSlaveonMasterSideY>length)
                {
                    coordinatesSlaveonMasterSideY-=length;
                    deltaRHS[0] = -0.5 * mStrain[2]*length-mCrackOpening[0];;
                    deltaRHS[1] = -mStrain[1]*length-mCrackOpening[1];;
                }
                else
                {
                    deltaRHS[0] = 0.;
                    deltaRHS[1] = 0.;
                }
            }

            while (coordinatesCurMaster[1]>coordinatesSlaveonMasterSideY || coordinatesNextMaster[1]<coordinatesSlaveonMasterSideY)
            {
                if (nextMasterNodePtr==mLeftUpperCorner)
                {
                    //restart from the left corner
                    curMasterNodePtr = mMasterNodesLeftBoundary[0];
                    coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

                    nextMasterNodecount=1;
                }
                else
                {
                    curMasterNodePtr = nextMasterNodePtr;
                    coordinatesCurMaster = coordinatesNextMaster;

                    nextMasterNodecount++;
                }
                nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
                coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }

            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlaveonMasterSideY);

            if (coordinatesCurMaster[1]>crackPosY)
            {
                deltaRHS[0]-=w*mCrackOpening[0];
                deltaRHS[1]-=w*mCrackOpening[1];
            }
            if (coordinatesNextMaster[1]>crackPosY)
            {
                deltaRHS[0]-=(1-w)*mCrackOpening[0];
                deltaRHS[1]-=(1-w)*mCrackOpening[1];
            }
            if (coordinatesSlaveonMasterSideY>crackPosY)
            {
                deltaRHS[0]+=mCrackOpening[0];
                deltaRHS[1]+=mCrackOpening[1];
            }

            //deltaDisp[0] = crackShift*0.5*mStrain[2] + length*mStrain[0] - deltaRHS[0];
            //deltaDisp[1] = crackShift*mStrain[1] + length *0.5* mStrain[2] - deltaRHS[1];

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 0),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 0),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddValue(curConstraintEquation,curSlaveNodePtr->GetDof(Node::DISPLACEMENTS, 1),1);
            if (std::abs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,curMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),-w);
            if (std::abs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddValue(curConstraintEquation,nextMasterNodePtr->GetDof(Node::DISPLACEMENTS, 1),w-1.);
            //rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }
    }
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearDisplacementsPeriodic2D::GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const
{
    Eigen::Matrix<double, 2, 1> LeftUpperCoordinates = mLeftUpperCorner->Get(Node::COORDINATES);
    Eigen::Matrix<double, 2, 1> LeftLowerCoordinates = mLeftLowerCorner->Get(Node::COORDINATES);

    //calculate length of specimen
    double length = LeftUpperCoordinates[1]-LeftLowerCoordinates[1];

    //only du linear interpolation between nodes on the boundary, this is not exact for quadrativ elements, but the effort is much reduced
    if (mAngle>=45 && mAngle<135)
    {
        double crackShift;
        crackShift = length*tan((90-mAngle)*PI/180.);
       //**************************************************************
        //add  constraints for all the slave nodes of the right boundary
        //**************************************************************
        assert(mMasterNodesLeftBoundary.size()>1);
        unsigned int nextMasterNodecount(1);

        NodeBase* curMasterNodePtr(mMasterNodesLeftBoundary[0]);
        Eigen::Matrix<double, 2, 1> coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        NodeBase* nextMasterNodePtr(mMasterNodesLeftBoundary[nextMasterNodecount]);;
        Eigen::Matrix<double, 2, 1> coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
        	Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideY;
            coordinatesSlaveonMasterSideY = coordinatesSlave[1];

            while (coordinatesNextMaster[1]<coordinatesSlaveonMasterSideY && nextMasterNodecount+1<mMasterNodesLeftBoundary.size())
            {
                 coordinatesCurMaster = coordinatesNextMaster;
                 nextMasterNodecount++;

                 assert (nextMasterNodecount<mMasterNodesLeftBoundary.size());
                 nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
                 coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }
            //slave is between two master nodes

            deltaDisp[0] = length*mStrain[0] + mCrackOpening[0];
            deltaDisp[1] = length*0.5*mStrain[2] + mCrackOpening[1];

            //constrain x direction
            rRHS(curConstraintEquation) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rRHS(curConstraintEquation) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the top boundary
        //**************************************************************

        assert(mMasterNodesBottomBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesBottomBoundary[0];
        coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
        coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double crackPosX((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
            Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideX;
            coordinatesSlaveonMasterSideX = coordinatesSlave[0]-crackShift;
            double deltaRHS[2];
            if (coordinatesSlaveonMasterSideX<0)
            {
                coordinatesSlaveonMasterSideX+=length;
                deltaRHS[0] = mStrain[0]*length+mCrackOpening[0];
                deltaRHS[1] = 0.5 * mStrain[2]*length+mCrackOpening[1];
            }
            else
            {
                if (coordinatesSlaveonMasterSideX>length)
                {
                    coordinatesSlaveonMasterSideX-=length;
                    deltaRHS[0] = -mStrain[0]*length-mCrackOpening[0];
                    deltaRHS[1] = -0.5 * mStrain[2]*length-mCrackOpening[1];
                }
                else
                {
                    deltaRHS[0] = 0.;
                    deltaRHS[1] = 0.;
                }
            }

            while (coordinatesCurMaster[0]>coordinatesSlaveonMasterSideX || coordinatesNextMaster[0]<coordinatesSlaveonMasterSideX)
            {
                if (nextMasterNodePtr==mRightLowerCorner)
                {
                    //restart from the left corner
                    curMasterNodePtr = mMasterNodesBottomBoundary[0];
                    coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

                    nextMasterNodecount = 1;
                }
                else
                {
                    coordinatesCurMaster = coordinatesNextMaster;
                    nextMasterNodecount++;
                }
                nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
                coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }

            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlaveonMasterSideX);

            if (coordinatesCurMaster[0]>crackPosX)
            {
                deltaRHS[0]-=w*mCrackOpening[0];
                deltaRHS[1]-=w*mCrackOpening[1];
            }
            if (coordinatesNextMaster[0]>crackPosX)
            {
                deltaRHS[0]-=(1-w)*mCrackOpening[0];
                deltaRHS[1]-=(1-w)*mCrackOpening[1];
            }
            if (coordinatesSlaveonMasterSideX>crackPosX)
            {
                deltaRHS[0]+=mCrackOpening[0];
                deltaRHS[1]+=mCrackOpening[1];
            }

            deltaDisp[0] = crackShift*mStrain[0] + length *0.5* mStrain[2] - deltaRHS[0];
            deltaDisp[1] = crackShift*0.5*mStrain[2] + length*mStrain[1] - deltaRHS[1];

            //constrain x direction
            rRHS(curConstraintEquation) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rRHS(curConstraintEquation) = deltaDisp[1];
            curConstraintEquation++;
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        //**************************************************************
        // angle between 135 abnd 225 degrees
        //**************************************************************
        // add  constraints for all the slave nodes of the top boundary
        //**************************************************************
        assert(mMasterNodesBottomBoundary.size()>1);
        unsigned int nextMasterNodecount(1);

        NodeBase* curMasterNodePtr(mMasterNodesBottomBoundary[0]);
        Eigen::Matrix<double, 2, 1> coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        NodeBase* nextMasterNodePtr(mMasterNodesBottomBoundary[nextMasterNodecount]);;
        Eigen::Matrix<double, 2, 1> coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
        	Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideX;
            coordinatesSlaveonMasterSideX = coordinatesSlave[0];

            while (coordinatesNextMaster[0]<coordinatesSlaveonMasterSideX && nextMasterNodecount+1<mMasterNodesBottomBoundary.size())
            {
                 coordinatesCurMaster = coordinatesNextMaster;
                 nextMasterNodecount++;

                 assert (nextMasterNodecount<mMasterNodesBottomBoundary.size());
                 nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
                 coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }
            //slave is between two master nodes

            deltaDisp[0] = length*0.5*mStrain[2] + mCrackOpening[0];;
            deltaDisp[1] = length*mStrain[1] + mCrackOpening[1];;

            //constrain x direction
            rRHS(curConstraintEquation) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rRHS(curConstraintEquation) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the right boundary
        //**************************************************************

        assert(mMasterNodesLeftBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesLeftBoundary[0];
        coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

        nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
        coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);

        double crackPosY((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
        	NodeBase* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
            Eigen::Matrix<double, 2, 1> coordinatesSlave = curSlaveNodePtr->Get(Node::COORDINATES);

            double coordinatesSlaveonMasterSideY;
            coordinatesSlaveonMasterSideY = coordinatesSlave[1]-crackShift;
            double deltaRHS[2];
            if (coordinatesSlaveonMasterSideY<0)
            {
                coordinatesSlaveonMasterSideY+=length;
                deltaRHS[0] = 0.5 * mStrain[2]*length+mCrackOpening[0];
                deltaRHS[1] = mStrain[1]*length+mCrackOpening[1];
            }
            else
            {
                if (coordinatesSlaveonMasterSideY>length)
                {
                    coordinatesSlaveonMasterSideY-=length;
                    deltaRHS[0] = -0.5 * mStrain[2]*length-mCrackOpening[0];;
                    deltaRHS[1] = -mStrain[1]*length-mCrackOpening[1];;
                }
                else
                {
                    deltaRHS[0] = 0.;
                    deltaRHS[1] = 0.;
                }
            }

            while (coordinatesCurMaster[1]>coordinatesSlaveonMasterSideY || coordinatesNextMaster[1]<coordinatesSlaveonMasterSideY)
            {
                if (nextMasterNodePtr==mLeftUpperCorner)
                {
                    //restart from the left corner
                    curMasterNodePtr = mMasterNodesLeftBoundary[0];
                    coordinatesCurMaster = curMasterNodePtr->Get(Node::COORDINATES);

                    nextMasterNodecount = 1;
                }
                else
                {
                    coordinatesCurMaster = coordinatesNextMaster;

                    nextMasterNodecount++;
                }
                nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
                coordinatesNextMaster = nextMasterNodePtr->Get(Node::COORDINATES);
            }

            //slave is between two master nodes

            //calculate weighting function for each master node
            double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlaveonMasterSideY);

            if (coordinatesCurMaster[1]>crackPosY)
            {
                deltaRHS[0]-=w*mCrackOpening[0];
                deltaRHS[1]-=w*mCrackOpening[1];
            }
            if (coordinatesNextMaster[1]>crackPosY)
            {
                deltaRHS[0]-=(1-w)*mCrackOpening[0];
                deltaRHS[1]-=(1-w)*mCrackOpening[1];
            }
            if (coordinatesSlaveonMasterSideY>crackPosY)
            {
                deltaRHS[0]+=mCrackOpening[0];
                deltaRHS[1]+=mCrackOpening[1];
            }

            deltaDisp[0] = crackShift*0.5*mStrain[2] + length*mStrain[0] - deltaRHS[0];
            deltaDisp[1] = crackShift*mStrain[1] + length *0.5* mStrain[2] - deltaRHS[1];

            //constrain x direction
            rRHS(curConstraintEquation) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rRHS(curConstraintEquation) = deltaDisp[1];
            curConstraintEquation++;
        }
    }
}

//calculate weighting function for each master node
double NuTo::ConstraintLinearDisplacementsPeriodic2D::CalculateWeightFunction(double rCoordinateCurMaster, double rCoordinateNextMaster, double rCoordinateSlave)const
{
    assert(rCoordinateNextMaster!=rCoordinateCurMaster);
    return 1.-(rCoordinateSlave-rCoordinateCurMaster)/(rCoordinateNextMaster-rCoordinateCurMaster);
}

//calculate delta displacement in x and y direction from the applied strain and the nodal position
void NuTo::ConstraintLinearDisplacementsPeriodic2D::CalculateDeltaDisp(double rCoordinates[2], double rDeltaDisp[2])const
{
    rDeltaDisp[0] = mStrain[0]*rCoordinates[0] + mStrain[2]*rCoordinates[1];
    rDeltaDisp[1] = mStrain[1]*rCoordinates[1] + mStrain[2]*rCoordinates[0];
}

//!@brief calculates all the nodes on the boundary
std::vector<NuTo::NodeBase*> NuTo::ConstraintLinearDisplacementsPeriodic2D::GetBoundaryNodes()
{
    Group<NodeBase> boundaryNodes;
    boundaryNodes.insert(mGroupTop->begin(),mGroupTop->end());
    boundaryNodes.insert(mGroupBottom->begin(),mGroupBottom->end());
    boundaryNodes.insert(mGroupLeft->begin(),mGroupLeft->end());
    boundaryNodes.insert(mGroupRight->begin(),mGroupRight->end());

    //write in a vector
    std::vector<NodeBase* > returnVector(boundaryNodes.GetNumMembers());
    int theNode(0);
    for (Group<NodeBase>::iterator itNode=boundaryNodes.begin(); itNode!=boundaryNodes.end();itNode++, theNode++)
    {
        returnVector[theNode] = itNode->second;
    }
    return returnVector;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearDisplacementsPeriodic2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearDisplacementsPeriodic2D::save(Archive & ar, const unsigned int version) const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearDisplacementsPeriodic2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);
    ar & BOOST_SERIALIZATION_NVP(mAngle);
    ar & BOOST_SERIALIZATION_NVP(mStrain);

    std::uintptr_t mLeftUpperCornerAddress = reinterpret_cast<std::uintptr_t>(mLeftUpperCorner);
    ar & boost::serialization::make_nvp("mLeftUpperCorner", mLeftUpperCornerAddress);

    std::uintptr_t mLeftLowerCornerAddress = reinterpret_cast<std::uintptr_t>(mLeftLowerCorner);
    ar & boost::serialization::make_nvp("mLeftLowerCorner", mLeftLowerCornerAddress);

    std::uintptr_t mRightUpperCornerAddress = reinterpret_cast<std::uintptr_t>(mRightUpperCorner);
    ar & boost::serialization::make_nvp("mRightUpperCorner", mRightUpperCornerAddress);

    std::uintptr_t mRightLowerCornerAddress = reinterpret_cast<std::uintptr_t>(mRightLowerCorner);
    ar & boost::serialization::make_nvp("mRightLowerCorner", mRightLowerCornerAddress);

    const std::uintptr_t* mSlaveNodesRightBoundaryAddress = reinterpret_cast<const std::uintptr_t*>(mSlaveNodesRightBoundary.data());
    int size = mSlaveNodesRightBoundary.size();
    ar & boost::serialization::make_nvp("mSlaveNodesRightBoundary_size", size);
    ar & boost::serialization::make_nvp("mSlaveNodesRightBoundary", boost::serialization::make_array(mSlaveNodesRightBoundaryAddress, size));

    const std::uintptr_t* mSlaveNodesTopBoundaryAddress = reinterpret_cast<const std::uintptr_t*>(mSlaveNodesTopBoundary.data());
    size = mSlaveNodesTopBoundary.size();
    ar & boost::serialization::make_nvp("mSlaveNodesTopBoundary_size", size);
    ar & boost::serialization::make_nvp("mSlaveNodesTopBoundary", boost::serialization::make_array(mSlaveNodesTopBoundaryAddress, size));

    const std::uintptr_t* mMasterNodesLeftBoundaryAddress = reinterpret_cast<const std::uintptr_t*>(mMasterNodesLeftBoundary.data());
    size = mMasterNodesLeftBoundary.size();
    ar & boost::serialization::make_nvp("mMasterNodesLeftBoundary_size", size);
    ar & boost::serialization::make_nvp("mMasterNodesLeftBoundary", boost::serialization::make_array(mMasterNodesLeftBoundaryAddress, size));

    const std::uintptr_t*  mMasterNodesBottomBoundaryAddress = reinterpret_cast<const std::uintptr_t*>( mMasterNodesBottomBoundary.data());
    size =  mMasterNodesBottomBoundary.size();
    ar & boost::serialization::make_nvp(" mMasterNodesBottomBoundary_size", size);
    ar & boost::serialization::make_nvp(" mMasterNodesBottomBoundary", boost::serialization::make_array( mMasterNodesBottomBoundaryAddress, size));

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearDisplacementsPeriodic2D" << std::endl;
#endif
}
template<class Archive>
void NuTo::ConstraintLinearDisplacementsPeriodic2D::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearDisplacementsPeriodic2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase);
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear);
    ar & BOOST_SERIALIZATION_NVP(mAngle);
    ar & BOOST_SERIALIZATION_NVP(mStrain);

    std::uintptr_t mLeftUpperCornerAddress;
    ar & boost::serialization::make_nvp("mLeftUpperCorner", mLeftUpperCornerAddress);
    mLeftUpperCorner = reinterpret_cast<const NodeBase*>(mLeftUpperCornerAddress);

    std::uintptr_t mLeftLowerCornerAddress;
    ar & boost::serialization::make_nvp("mLeftLowerCorner", mLeftLowerCornerAddress);
    mLeftLowerCorner = reinterpret_cast<const NodeBase*>(mLeftLowerCornerAddress);

    std::uintptr_t mRightUpperCornerAddress;
    ar & boost::serialization::make_nvp("mRightUpperCorner", mRightUpperCornerAddress);
    mRightUpperCorner = reinterpret_cast<const NodeBase*>(mRightUpperCornerAddress);

    std::uintptr_t mRightLowerCornerAddress;
    ar & boost::serialization::make_nvp("mRightLowerCorner", mRightLowerCornerAddress);
    mRightLowerCorner = reinterpret_cast<const NodeBase*>(mRightLowerCornerAddress);

    int size = 0;
    ar & boost::serialization::make_nvp("mSlaveNodesRightBoundary_size", size);
    std::uintptr_t* mSlaveNodesRightBoundaryAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mSlaveNodesRightBoundary", boost::serialization::make_array(mSlaveNodesRightBoundaryAddress, size));
    mSlaveNodesRightBoundary.assign(reinterpret_cast<NodeBase**>(&mSlaveNodesRightBoundaryAddress[0]), reinterpret_cast<NodeBase**>(&mSlaveNodesRightBoundaryAddress[size]));

    size = 0;
    ar & boost::serialization::make_nvp("mSlaveNodesTopBoundary_size", size);
    std::uintptr_t*  mSlaveNodesTopBoundaryAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mSlaveNodesTopBoundary", boost::serialization::make_array( mSlaveNodesTopBoundaryAddress, size));
     mSlaveNodesTopBoundary.assign(reinterpret_cast<NodeBase**>(& mSlaveNodesTopBoundaryAddress[0]), reinterpret_cast<NodeBase**>(& mSlaveNodesTopBoundaryAddress[size]));

    size = 0;
    ar & boost::serialization::make_nvp("mMasterNodesLeftBoundary_size", size);
    std::uintptr_t*  mMasterNodesLeftBoundaryAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mMasterNodesLeftBoundary", boost::serialization::make_array( mMasterNodesLeftBoundaryAddress, size));
     mMasterNodesLeftBoundary.assign(reinterpret_cast<NodeBase**>(& mMasterNodesLeftBoundaryAddress[0]), reinterpret_cast<NodeBase**>(& mMasterNodesLeftBoundaryAddress[size]));

    size = 0;
    ar & boost::serialization::make_nvp("mMasterNodesBottomBoundary_size", size);
    std::uintptr_t*  mMasterNodesBottomBoundaryAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mMasterNodesBottomBoundary", boost::serialization::make_array( mMasterNodesBottomBoundaryAddress, size));
     mMasterNodesBottomBoundary.assign(reinterpret_cast<NodeBase**>(& mMasterNodesBottomBoundaryAddress[0]), reinterpret_cast<NodeBase**>(& mMasterNodesBottomBoundaryAddress[size]));

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearDisplacementsPeriodic2D" << std::endl;
#endif
}

void NuTo::ConstraintLinearDisplacementsPeriodic2D::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    std::map<std::uintptr_t, std::uintptr_t>::const_iterator it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mLeftUpperCorner));
    if (it!=mNodeMapCast.end())
    {
        NodeBase** temp = const_cast<NodeBase**>(&mLeftUpperCorner);
        *temp = reinterpret_cast<NodeBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] The mLeftUpperCorner NodeBase-Pointer could not be updated.");

    it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mLeftLowerCorner));
    if (it!=mNodeMapCast.end())
    {
        NodeBase** temp = const_cast<NodeBase**>(&mLeftLowerCorner);
        *temp = reinterpret_cast<NodeBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] The mLeftLowerCorner NodeBase-Pointer could not be updated.");

    it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mRightUpperCorner));
    if (it!=mNodeMapCast.end())
    {
        NodeBase** temp = const_cast<NodeBase**>(&mRightUpperCorner);
        *temp = reinterpret_cast<NodeBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] The mRightUpperCorner NodeBase-Pointer could not be updated.");

    it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mRightLowerCorner));
    if (it!=mNodeMapCast.end())
    {
        NodeBase** temp = const_cast<NodeBase**>(&mRightLowerCorner);
        *temp = reinterpret_cast<NodeBase*>(it->second);
    }
    else
        throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] The mRightLowerCorner NodeBase-Pointer could not be updated.");

    for(std::vector<NodeBase*>::iterator it = mSlaveNodesRightBoundary.begin(); it != mSlaveNodesRightBoundary.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&(*it));
            *temp = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] A NodeBase-Pointer in mSlaveNodesRightBoundary could not be updated.");
    }

    for(std::vector<NodeBase*>::iterator it = mSlaveNodesTopBoundary.begin(); it != mSlaveNodesTopBoundary.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&(*it));
            *temp = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] A NodeBase-Pointer in mSlaveNodesTopBoundary could not be updated.");
    }

    for(std::vector<NodeBase*>::iterator it = mMasterNodesLeftBoundary.begin(); it != mMasterNodesLeftBoundary.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&(*it));
            *temp = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] A NodeBase-Pointer in mMasterNodesLeftBoundary could not be updated.");
    }

    for(std::vector<NodeBase*>::iterator it = mMasterNodesBottomBoundary.begin(); it != mMasterNodesBottomBoundary.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&(*it));
            *temp = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::ConstraintLinearDisplacementsPeriodic2D] A NodeBase-Pointer in mMasterNodesBottomBoundary could not be updated.");
    }

}


BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearDisplacementsPeriodic2D)
#endif // ENABLE_SERIALIZATION
