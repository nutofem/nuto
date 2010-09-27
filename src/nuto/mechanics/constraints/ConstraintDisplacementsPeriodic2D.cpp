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

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeDisplacements2D.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintDisplacementsPeriodic2D.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief constructor
NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D(const StructureBase* rStructure, double rAngle,
        NuTo::FullMatrix<double> rStrain,NuTo::FullMatrix<double> rCrackOpening, double rRadiusToCrackWithoutConstraints,
        const Group<NodeBase>* rGroupTop,const Group<NodeBase>* rGroupBottom,
        const Group<NodeBase>* rGroupLeft, const Group<NodeBase>* rGroupRight) :  ConstraintBase()
{
    mStructure = rStructure,

    mRadiusToCrackWithoutConstraints = rRadiusToCrackWithoutConstraints;

    SetCrackOpening(rCrackOpening);

    if (rStrain.GetNumRows()!=3 || rStrain.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::ConstraintNodeDisplacementsPeriodic2D::ConstraintNodeDisplacementsPeriodic2D] the strain is matrix (3,1) with (e_xx, e_yy, gamma_xy)");
    mStrain[0] = rStrain(0,0);
    mStrain[1] = rStrain(1,0);
    mStrain[2] = rStrain(2,0);

    mGroupTop = rGroupTop;
    mGroupBottom = rGroupBottom;
    mGroupLeft = rGroupLeft;
    mGroupRight = rGroupRight;

    //determine the corner nodes
    Group<NodeBase>* newGroup = dynamic_cast<Group<NodeBase>*>(rGroupTop->Intersection (rGroupLeft));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Top left corner node can not be determined.");
    mLeftUpperCorner = *(newGroup->begin());
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupTop->Intersection (rGroupRight));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Top right corner node can not be determined.");
    mRightUpperCorner = *(newGroup->begin());
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupBottom->Intersection (rGroupLeft));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Bottom left corner node can not be determined.");
    mLeftLowerCorner = *(newGroup->begin());
    delete newGroup;

    newGroup = dynamic_cast<Group<NodeBase>*>(rGroupBottom->Intersection (rGroupRight));
    if (newGroup->GetNumMembers()!=1)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Bottom right corner node can not be determined.");
    mRightLowerCorner = *(newGroup->begin());
    delete newGroup;

    //get box coordinates
    double LeftUpperCoordinates[2], RightUpperCoordinates[2], LeftLowerCoordinates[2], RightLowerCoordinates[2];
    if (mLeftUpperCorner->GetNumCoordinates()!=2)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Upper left node does not have 2 coordinates.");
    mLeftUpperCorner->GetCoordinates2D(LeftUpperCoordinates);

    if (mRightUpperCorner->GetNumCoordinates()!=2)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Upper right node does not have 2 coordinates.");
    mRightUpperCorner->GetCoordinates2D(RightUpperCoordinates);

    if (mLeftLowerCorner->GetNumCoordinates()!=2)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Lower left node does not have 2 coordinates.");
    mLeftLowerCorner->GetCoordinates2D(LeftLowerCoordinates);

    if (mRightLowerCorner->GetNumCoordinates()!=2)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Lower right node does not have 2 coordinates.");
    mRightLowerCorner->GetCoordinates2D(RightLowerCoordinates);

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
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Left boundary is not correct - check your node groups.");
    }
    if (RightUpperCoordinates[0]!=RightLowerCoordinates[0])
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] Right boundary is not correct - check your node groups.");

    if (RightUpperCoordinates[0]<=LeftUpperCoordinates[0])
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] left boundary coordinate is larger than right boundary coordinate.");

    if (LeftUpperCoordinates[1]!=RightUpperCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] upper boundary is not correct - check your node groups.");

    if (LeftLowerCoordinates[0]!=RightLowerCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] lower boundary is not correct - check your node groups.");

    if (LeftUpperCoordinates[1]<=LeftLowerCoordinates[1])
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::ConstraintDisplacementsPeriodic2D] upper boundary coordinate is larger than lower boundary coordinate.");


    SetAngle(rAngle);
}

//!@brief set the angle of the periodic boundary conditions
//!@param rRHS new right hand side
void NuTo::ConstraintDisplacementsPeriodic2D::SetAngle(double rAngle)
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
void NuTo::ConstraintDisplacementsPeriodic2D::SetCrackOpening(const NuTo::FullMatrix<double>& rCrackOpening)
{
    if (rCrackOpening.GetNumRows()!=2 || rCrackOpening.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::ConstraintDisplacementsPeriodic2D::SetCrackOpening] crack opening should be a (2,1) matrix.");
    mCrackOpening[0] = rCrackOpening(0,0);
    mCrackOpening[1] = rCrackOpening(1,0);
}

//!@brief set the strain of the periodic boundary conditions
//!@param rStrain strain (e_xx,e_yy,gamma_xy)
void NuTo::ConstraintDisplacementsPeriodic2D::SetStrain(const NuTo::FullMatrix<double>& rStrain)
{
    if (rStrain.GetNumRows()!=3 || rStrain.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::ConstraintNodeDisplacementsPeriodic2D::ConstraintNodeDisplacementsPeriodic2D] the strain is matrix (3,1) with (e_xx, e_yy, gamma_xy)");
    mStrain[0] = rStrain(0,0);
    mStrain[1] = rStrain(1,0);
    mStrain[2] = rStrain(2,0);
}

#define PI 3.14159265359

//!@brief calculate the border vectors in counterclockwise direction
void NuTo::ConstraintDisplacementsPeriodic2D::SetBoundaryVectors()
{
    double LeftUpperCoordinates[2], LeftLowerCoordinates[2];
    mLeftUpperCorner->GetCoordinates2D(LeftUpperCoordinates);
    mLeftLowerCorner->GetCoordinates2D(LeftLowerCoordinates);

    //calculate length of specimen
    double length = LeftUpperCoordinates[1]-LeftLowerCoordinates[1];
    //calculate master nodes left boundary (green)
    mMasterNodesLeftBoundary.resize(0);
    mMasterNodesLeftBoundary.reserve(mGroupLeft->GetNumMembers());
    for (Group<NodeBase>::iterator itNode=mGroupLeft->begin(); itNode!=mGroupLeft->end();itNode++)
    {
        mMasterNodesLeftBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode));
    }
    sort(mMasterNodesLeftBoundary.begin(), mMasterNodesLeftBoundary.end(), less_YCoordinate2D());

    //calculate master nodes bottom boundary (yellow-orange-blue)
    mMasterNodesBottomBoundary.resize(0);
    mMasterNodesBottomBoundary.reserve(mGroupBottom->GetNumMembers());
    for (Group<NodeBase>::iterator itNode=mGroupBottom->begin(); itNode!=mGroupBottom->end();itNode++)
    {
        mMasterNodesBottomBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode));
    }
    sort(mMasterNodesBottomBoundary.begin(), mMasterNodesBottomBoundary.end(), less_XCoordinate2D());

    //calculate slave nodes right boundary (green)
    mSlaveNodesRightBoundary.resize(0);
    if (mAngle>=45 && mAngle<135)
    {
        double crackShift;
        crackShift = length*tan((90-mAngle)*PI/180.);
        mSlaveNodesRightBoundary.reserve(mGroupRight->GetNumMembers());
        double coordinates[3];
        for (Group<NodeBase>::iterator itNode=mGroupRight->begin(); itNode!=mGroupRight->end();itNode++)
        {
            NodeCoordinatesDisplacements2D* nodePtr = dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode);
            nodePtr->GetCoordinates2D(coordinates);
            double DeltaX((length-crackShift)*0.5);
            double DeltaY(length-coordinates[1]);

            if (DeltaX*DeltaX + DeltaY*DeltaY >=mRadiusToCrackWithoutConstraints * mRadiusToCrackWithoutConstraints)
                mSlaveNodesRightBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode));
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        mSlaveNodesRightBoundary.reserve(mGroupRight->GetNumMembers()-1);
        double coordinates[3];
        for (Group<NodeBase>::iterator itNode=mGroupRight->begin(); itNode!=mGroupRight->end();itNode++)
        {
            if ((*itNode)!=mRightLowerCorner)
            {
                NodeCoordinatesDisplacements2D* nodePtr = dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode);
                nodePtr->GetCoordinates2D(coordinates);
                if (fabs(coordinates[1]-(length+crackShift)*0.5)>=mRadiusToCrackWithoutConstraints)
                    mSlaveNodesRightBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode));
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
        double coordinates[3];
        for (Group<NodeBase>::iterator itNode=mGroupTop->begin(); itNode!=mGroupTop->end();itNode++)
        {
            if ((*itNode)!=mLeftUpperCorner)
            {
                NodeCoordinatesDisplacements2D* nodePtr = dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode);
                nodePtr->GetCoordinates2D(coordinates);
                if (fabs(coordinates[0]-(length+crackShift)*0.5)>=mRadiusToCrackWithoutConstraints)
                    mSlaveNodesTopBoundary.push_back(nodePtr);
            }
        }
    }
    else
    {
        double crackShift;
        crackShift = length*tan((mAngle)*PI/180.);
        mSlaveNodesTopBoundary.reserve(mGroupTop->GetNumMembers());
        double coordinates[3];
        for (Group<NodeBase>::iterator itNode=mGroupTop->begin(); itNode!=mGroupTop->end();itNode++)
        {
            NodeCoordinatesDisplacements2D* nodePtr = dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode);
            nodePtr->GetCoordinates2D(coordinates);
            double DeltaX(coordinates[0]-(length+crackShift)*0.5);
            double DeltaY((length-crackShift)*0.5);

            if (DeltaX*DeltaX + DeltaY*DeltaY >=mRadiusToCrackWithoutConstraints * mRadiusToCrackWithoutConstraints)
                mSlaveNodesTopBoundary.push_back(dynamic_cast<NodeCoordinatesDisplacements2D*>(*itNode));
        }

    }
    sort(mSlaveNodesTopBoundary.begin(), mSlaveNodesTopBoundary.end(), less_XCoordinate2D());


    //Info about the nodes
    if (mStructure->GetVerboseLevel()>0)
    {
        for (int countBoundary=0; countBoundary<4; countBoundary++)
        {
            std::vector<NodeCoordinatesDisplacements2D*>* nodeVectorPtr;
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
                std::cout << "  " << mStructure->NodeGetId((*nodeVectorPtr)[countNodes]) << ": " << coordinates[0] << " " << coordinates[1] << std::endl;
            }
        }
    }
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
#define MIN_CONSTRAINT 1e-6
void NuTo::ConstraintDisplacementsPeriodic2D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    double LeftUpperCoordinates[2], LeftLowerCoordinates[2];
    mLeftUpperCorner->GetCoordinates2D(LeftUpperCoordinates);
    mLeftLowerCorner->GetCoordinates2D(LeftLowerCoordinates);

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

        NodeCoordinatesDisplacements2D* curMasterNodePtr(mMasterNodesLeftBoundary[0]);
        double coordinatesCurMaster[2];
        curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

        NodeCoordinatesDisplacements2D* nextMasterNodePtr(mMasterNodesLeftBoundary[nextMasterNodecount]);;
        double coordinatesNextMaster[2];
        nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

        double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
            NodeCoordinatesDisplacements2D* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
            double coordinatesSlave[2];
            curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

            double coordinatesSlaveonMasterSideY;
            coordinatesSlaveonMasterSideY = coordinatesSlave[1];

            while (coordinatesNextMaster[1]<coordinatesSlaveonMasterSideY && nextMasterNodecount+1<mMasterNodesLeftBoundary.size())
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
            double w = this->CalculateWeightFunction(coordinatesCurMaster[1],coordinatesNextMaster[1],coordinatesSlaveonMasterSideY);

            deltaDisp[0] = length*mStrain[0] + mCrackOpening[0];
            deltaDisp[1] = length*0.5*mStrain[2] + mCrackOpening[1];

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << 1-w << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(0),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(0),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(0),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(1),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(1),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(1),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the top boundary
        //**************************************************************

        assert(mMasterNodesBottomBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesBottomBoundary[0];
        curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

        nextMasterNodePtr = mMasterNodesBottomBoundary[nextMasterNodecount];
        nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

        double crackPosX((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
            NodeCoordinatesDisplacements2D* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
            double coordinatesSlave[2];
            curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

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
                    curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

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
                nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
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

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(0),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(0),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(0),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(1),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(1),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(1),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[1];
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

        NodeCoordinatesDisplacements2D* curMasterNodePtr(mMasterNodesBottomBoundary[0]);
        double coordinatesCurMaster[2];
        curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

        NodeCoordinatesDisplacements2D* nextMasterNodePtr(mMasterNodesBottomBoundary[nextMasterNodecount]);;
        double coordinatesNextMaster[2];
        nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

        double deltaDisp[2];
        for (unsigned int countNode=0; countNode<mSlaveNodesTopBoundary.size(); countNode++)
        {
            NodeCoordinatesDisplacements2D* curSlaveNodePtr(mSlaveNodesTopBoundary[countNode]);
            double coordinatesSlave[2];
            curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

            double coordinatesSlaveonMasterSideX;
            coordinatesSlaveonMasterSideX = coordinatesSlave[0];

            while (coordinatesNextMaster[0]<coordinatesSlaveonMasterSideX && nextMasterNodecount+1<mMasterNodesBottomBoundary.size())
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
            double w = this->CalculateWeightFunction(coordinatesCurMaster[0],coordinatesNextMaster[0],coordinatesSlaveonMasterSideX);

            deltaDisp[0] = length*0.5*mStrain[2] + mCrackOpening[0];;
            deltaDisp[1] = length*mStrain[1] + mCrackOpening[1];;

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(0),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(0),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(0),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(1),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(1),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(1),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }

        //**************************************************************
        //add  constraints for all the slave nodes of the right boundary
        //**************************************************************

        assert(mMasterNodesLeftBoundary.size()>1);
        nextMasterNodecount = 1;

        curMasterNodePtr = mMasterNodesLeftBoundary[0];
        curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

        nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
        nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);

        double crackPosY((length-crackShift)*0.5);

        for (unsigned int countNode=0; countNode<mSlaveNodesRightBoundary.size(); countNode++)
        {
            NodeCoordinatesDisplacements2D* curSlaveNodePtr(mSlaveNodesRightBoundary[countNode]);
            double coordinatesSlave[2];
            curSlaveNodePtr->GetCoordinates2D(coordinatesSlave);

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
                    curMasterNodePtr->GetCoordinates2D(coordinatesCurMaster);

                    nextMasterNodecount=1;
                }
                else
                {
                    curMasterNodePtr = nextMasterNodePtr;
                    coordinatesCurMaster[0] = coordinatesNextMaster[0];
                    coordinatesCurMaster[1] = coordinatesNextMaster[1];

                    nextMasterNodecount++;
                }
                nextMasterNodePtr = mMasterNodesLeftBoundary[nextMasterNodecount];
                nextMasterNodePtr->GetCoordinates2D(coordinatesNextMaster);
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

/*            std::cout << "constraint equation " << curConstraintEquation
                    << ": node " << mStructure->NodeGetId(curSlaveNodePtr) << " + "
                    << -w << " node " << mStructure->NodeGetId(curMasterNodePtr) << " + "
                    << w-1 << " node " << mStructure->NodeGetId(nextMasterNodePtr)
                    << " = (" << deltaDisp[0] << ", " << deltaDisp[1] << ")"
                    << std::endl;
*/
            //constrain x direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(0),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(0),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(0),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[0];
            curConstraintEquation++;

            //constrain y direction
            rConstraintMatrix.AddEntry(curConstraintEquation,curSlaveNodePtr->GetDofDisplacement(1),1);
            if (fabs(w)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,curMasterNodePtr->GetDofDisplacement(1),-w);
            if (fabs(w-1.)>MIN_CONSTRAINT)
                rConstraintMatrix.AddEntry(curConstraintEquation,nextMasterNodePtr->GetDofDisplacement(1),w-1.);
            rRHS(curConstraintEquation,0) = deltaDisp[1];
            curConstraintEquation++;
        }
    }
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintDisplacementsPeriodic2D::GetNumConstraintEquations()const
{
    return 2*(mSlaveNodesRightBoundary.size()+mSlaveNodesTopBoundary.size());
}

//calculate weighting function for each master node
double NuTo::ConstraintDisplacementsPeriodic2D::CalculateWeightFunction(double rCoordinateCurMaster, double rCoordinateNextMaster, double rCoordinateSlave)const
{
    assert(rCoordinateNextMaster!=rCoordinateCurMaster);
    return 1.-(rCoordinateSlave-rCoordinateCurMaster)/(rCoordinateNextMaster-rCoordinateCurMaster);
}

//calculate delta displacement in x and y direction from the applied strain and the nodal position
void NuTo::ConstraintDisplacementsPeriodic2D::CalculateDeltaDisp(double rCoordinates[2], double rDeltaDisp[2])const
{
    rDeltaDisp[0] = mStrain[0]*rCoordinates[0] + mStrain[2]*rCoordinates[1];
    rDeltaDisp[1] = mStrain[1]*rCoordinates[1] + mStrain[2]*rCoordinates[0];
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintDisplacementsPeriodic2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintDisplacementsPeriodic2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintDisplacementsPeriodic2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_NVP(mAngle)
       & BOOST_SERIALIZATION_NVP(mStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintDisplacementsPeriodic2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintDisplacementsPeriodic2D)
#endif // ENABLE_SERIALIZATION
