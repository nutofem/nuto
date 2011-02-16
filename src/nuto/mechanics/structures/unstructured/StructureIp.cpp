// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION
#include <cmath>

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalCrackAngle.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalCrackOpening.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalTotalStrain.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackOpening2D.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackAngle2D.h"
#include "nuto/mechanics/constraints/ConstraintNonlinearGlobalCrackAngle2D.h"
#include "nuto/mechanics/constraints/ConstraintNonlinearGlobalCrackOpeningTangential2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsMultiscale2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/groups/GroupEnum.h"

#include <ANN/ANN.h>
#include <set>

#define PRINTRESULT true
//#define MYDEBUG
#define tolCrackOpeningConstraint 1e-5
//! @brief constructor
//! @param mDimension number of nodes
NuTo::StructureIp::StructureIp ( int rDimension )  : Structure ( rDimension )
{
    if (rDimension!=2)
        throw MechanicsException("[NuTo::StructureIp::StructureIp] The concurrent multiscale model is only implemented for 2D.");
    mCrackAngle = M_PI*(0.7000);
    mPrevCrackAngle = 0.5*M_PI;
    mDOFCrackAngle = -1;
    mCrackOpening[0] = 0.0;
    mCrackOpening[1] = 0.0;
    mDOFCrackOpening[0] = -1;
    mDOFCrackOpening[1] = -1;
    mEpsilonHom.mEngineeringStrain[0] = 0.;
    mEpsilonHom.mEngineeringStrain[1] = 0.;
    mEpsilonHom.mEngineeringStrain[2] = 0.;
    mDOFGlobalTotalStrain[0] = -1;
    mDOFGlobalTotalStrain[1] = -1;
    mDOFGlobalTotalStrain[2] = -1;
    mConstraintFineScaleX = -1;
    mConstraintFineScaleY = -1;
    mConstraintCrackAngle = -1;
    mConstraintNormalCrackOpening = -1;
    mConstraintTangentialCrackOpening = -1;
     mlX = 0.;
    mlY = 0.;
    mCrackTransitionZone = 0;
    mIPName = std::string("fineScaleIp");
    mGroupBoundaryNodes = -1;
    mBoundaryNodesAssigned = false;

    mToleranceElasticCrackAngleLow = 0.0;
    mToleranceElasticCrackAngleHigh = 0.0;

    // generate the constrain equation
    mConstraintTotalStrain = ConstraintLinearGlobalTotalStrain(NuTo::EngineeringStrain2D());

    //add constraint equation to avoid negative normal Crack opening
    //find unused integer id
/*    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    mConstraintCrackOpening = id;
    ConstraintLagrangeGlobalCrackOpening2D *mConst = new NuTo::ConstraintLagrangeGlobalCrackOpening2D(this, 0);
    mConstraintMap.insert(id, mConst);
*/
/*    NuTo::FullMatrix<double> dir(2,1);
    dir(0,0)=1;
    dir(1,0)=0;
    ConstraintLinearGlobalCrackOpening *mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this,dir,0);
    mConstraintMap.insert(id, mConst);
    dir(0,0)=0;
    dir(1,0)=1;
    ConstraintLinearGlobalCrackOpening *mConstb = new NuTo::ConstraintLinearGlobalCrackOpening(this,dir,1);
    id++;
    mConstraintMap.insert(id, mConstb);
*/
/*
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0) = 1;
    direction(1,0) = 0;
    ConstraintLinearGlobalCrackOpening *mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction,0.);
    mConstraintMap.insert(id, mConst);
*/
/*    direction(0,0) = 1;
    direction(1,0) = 0;
    mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction,0.1);
    id++;
    mConstraintMap.insert(id, mConst);
*/
    //add constraint equation to avoid singular stiffness for zero crack opening (because then the angle has no influence)
    //find unused integer id
/*    id = 1;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    ConstraintLagrangeGlobalCrackAngle2D *mConst2 = new NuTo::ConstraintLagrangeGlobalCrackAngle2D(this);
    //ConstraintLinearGlobalCrackAngle *mConst2 = new NuTo::ConstraintLinearGlobalCrackAngle(this);

    mConstraintMap.insert(id, mConst2);
    mConstraintCrackAngle = id;
*/
    std::cout << "number of constraints after reading in conversion " << mConstraintMap.size() << std::endl;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureIp::Info()const
{
    Structure::Info();
    //add crack angle and crack orientation

}

#ifdef ENABLE_SERIALIZATION
template void NuTo::StructureIp::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::StructureIp::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of structureIp" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Structure);
    ar & BOOST_SERIALIZATION_NVP(mCrackAngle)
       & BOOST_SERIALIZATION_NVP(mDOFCrackAngle)
       & BOOST_SERIALIZATION_NVP(mCrackOpening)
       & BOOST_SERIALIZATION_NVP(mDOFCrackOpening)
       & BOOST_SERIALIZATION_NVP(mEpsilonHom)
       & BOOST_SERIALIZATION_NVP(mDOFGlobalTotalStrain)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleX)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleY)
       & BOOST_SERIALIZATION_NVP(mConstraintNormalCrackOpening)
       & BOOST_SERIALIZATION_NVP(mConstraintTangentialCrackOpening)
       & BOOST_SERIALIZATION_NVP(mConstraintCrackAngle)
       & BOOST_SERIALIZATION_NVP(mConstraintTotalStrain)
       & BOOST_SERIALIZATION_NVP(mlX)
       & BOOST_SERIALIZATION_NVP(mlY)
       & BOOST_SERIALIZATION_NVP(mCrackTransitionZone)
       & BOOST_SERIALIZATION_NVP(mBoundaryNodesAssigned)
       & BOOST_SERIALIZATION_NVP(mGroupBoundaryNodes)
       & BOOST_SERIALIZATION_NVP(mIPName)
       & BOOST_SERIALIZATION_NVP(mToleranceElasticCrackAngleLow)
       & BOOST_SERIALIZATION_NVP(mToleranceElasticCrackAngleHigh);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structureIp" << std::endl;
#endif
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureIp::Save (const std::string &filename, std::string rType )const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::StructureIp::Save] Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ("Object_type", typeIdString );
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            std::string tmpString(this->GetTypeId());
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::StructureIp::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureIp::Save]File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureIp::Restore (const std::string &filename, std::string rType )
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::StructureIp::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::StructureIp::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureIp::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureIp::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
#endif// ENABLE_SERIALIZATION


void NuTo::StructureIp::SaveStructure()const
{
#ifdef ENABLE_SERIALIZATION
    if (mSavedToStringStream==false)
    {
        boost::archive::binary_oarchive oba(mSaveStringStream);
        oba << (*this);
        mSavedToStringStream = true;
    }
#else
    throw MechanicsException("[NuTo::StructureIp::Save] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

void NuTo::StructureIp::RestoreStructure()
{
#ifdef ENABLE_SERIALIZATION
    if (mSavedToStringStream==true)
    {
        boost::archive::binary_iarchive iba(mSaveStringStream);
        iba >> (*this);
        //delete string stream afterwards
        mSaveStringStream.str("");
    }
#else
    throw MechanicsException("[NuTo::StructureIp::Restore] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

//! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
//! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
void NuTo::StructureIp::TransformBoundaryNodes()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (mBoundaryNodesAssigned==false)
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] Group of boundary nodes has not been assigned yet.");
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(mGroupBoundaryNodes);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] Group for boundary nodes is not a node group.");

    Group<NodeBase>* nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);

    //copy the boundary nodes first, since the iterators are invalidated in NodeExchangePtr
    std::vector<NodeBase*> nodeVec(nodeGroup->size());
    unsigned int countNode=0;
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end(); itNode++, countNode++)
    {
        nodeVec[countNode] = (*itNode);
    }
    double minX(0.), minY(0.), maxX(0.), maxY(0.);

    for (countNode=0; countNode<nodeVec.size(); countNode++)
    {
        Node::eNodeType nodeType = nodeVec[countNode]->GetNodeType();
        NodeBase* newNode(0);

        switch (nodeType)
        {
        case Node::NodeCoordinatesDisplacements2D:
            newNode = new NodeCoordinatesDisplacementsMultiscale2D(this);
            double coordinates[2];
            nodeVec[countNode]->GetCoordinates2D(coordinates);
            if (countNode==0)
            {
                minX = coordinates[0];
                maxX = coordinates[0];
                minY = coordinates[1];
                maxY = coordinates[1];
            }
            else
            {
                if (minX > coordinates[0])
                    minX = coordinates[0];
                if (maxX < coordinates[0])
                    maxX = coordinates[0];
                if (minY > coordinates[1])
                    minY = coordinates[1];
                if (maxY < coordinates[1])
                    maxY = coordinates[1];
            }


            newNode->SetCoordinates2D(coordinates);
        case Node::NodeCoordinatesDisplacementsMultiscale2D:
            break;
        break;
        default:
            throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] node type not implemented.");
        }

        //old node is deleted in Exchange routine when being removed from node table
        NodeExchangePtr(nodeVec[countNode],newNode);
    }
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0) = 1.;
    direction(1,0) = 0.;
    mConstraintFineScaleX = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    direction(0,0) = 0.;
    direction(1,0) = 1.;
    mConstraintFineScaleY = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    ConstraintInfo(10);

    if (minX!=0 || minY!=0)
    {
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] left lower corner has to be at the origin.");
    }
    if (minX==maxX || minY==maxY)
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] structure has zero width or height, either check your boundary group or the full structure.");
    mlX = maxX;
    mlY = maxY;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureIp::TransformBoundaryNodes] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief numbers non standard DOFs' e.g. in StructureIp, for standard structures this routine is empty
void NuTo::StructureIp::NumberAdditionalGlobalDofs()
{
    // DOFs related to the crack angle and crack opening
    if (mDimension==2)
    {
        mDOFCrackAngle = mNumDofs++;
        mDOFCrackOpening[0] = mNumDofs++;
        mDOFCrackOpening[1] = mNumDofs++;
        mDOFGlobalTotalStrain[0] = mNumDofs++;
        mDOFGlobalTotalStrain[1] = mNumDofs++;
        mDOFGlobalTotalStrain[2] = mNumDofs++;
    }
    else
        throw MechanicsException("[NuTo::StructureIp::SetAdditionalGlobalDofs] Only implemented for 2D.");
}

// merge dof values
void NuTo::StructureIp::NodeMergeAdditionalGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    //merge alpha
    if (mDimension==2)
    {
        double value;
        int dof = this->mDOFCrackAngle;
        if (mDOFCrackAngle >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof,0);
        }
        else
        {
            value = rActiveDofValues(dof,0);
        }
        while (value<0)
        {
            value+=2.*M_PI;
        }
        while (value>2.*M_PI)
        {
            value-=2.*M_PI;
        }
        this->mCrackAngle = value;


        //merge crack opening
        for (int count=0; count<2; count++)
        {
            int dof = this->mDOFCrackOpening[count];
            double value;
            if (dof >= rActiveDofValues.GetNumRows())
            {
                dof -= rActiveDofValues.GetNumRows();
                assert(dof < rDependentDofValues.GetNumRows());
                value = rDependentDofValues(dof,0);
            }
            else
            {
                value = rActiveDofValues(dof,0);
            }
            this->mCrackOpening[count] = value;
        }

        //merge total strain
        for (int count=0; count<3; count++)
        {
            int dof = this->mDOFGlobalTotalStrain[count];
            double value;
            if (dof >= rActiveDofValues.GetNumRows())
            {
                dof -= rActiveDofValues.GetNumRows();
                assert(dof < rDependentDofValues.GetNumRows());
                value = rDependentDofValues(dof,0);
            }
            else
            {
                value = rActiveDofValues(dof,0);
            }
            this->mEpsilonTot.mEngineeringStrain[count] = value;
        }

    }
    else
        throw MechanicsException("[NuTo::StructureIp::NodeMergeActiveDofValues] Only implemented for 2D");

    //update the homogeneous strain
    CalculateHomogeneousEngineeringStrain();
}

//! @brief extract dof values additional dof values
//! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
//! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
void NuTo::StructureIp::NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const
{
    if (mDimension==2)
    {
        int dof = this->mDOFCrackAngle;
        double value = this->mCrackAngle;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof,0) = value;
        }
        else
        {
            rActiveDofValues(dof,0) = value;
        }
        for (int count=0; count<2; count++)
        {
            int dof = this->mDOFCrackOpening[count];
            double value = this->mCrackOpening[count];
            if (dof >= rActiveDofValues.GetNumRows())
            {
                dof -= rActiveDofValues.GetNumRows();
                assert(dof < rDependentDofValues.GetNumRows());
                rDependentDofValues(dof,0) = value;
            }
            else
            {
                rActiveDofValues(dof,0) = value;
            }
        }
        for (int count=0; count<3; count++)
        {
            int dof = this->mDOFGlobalTotalStrain[count];
            double value = this->mEpsilonTot.mEngineeringStrain[count];
            if (dof >= rActiveDofValues.GetNumRows())
            {
                dof -= rActiveDofValues.GetNumRows();
                assert(dof < rDependentDofValues.GetNumRows());
                rDependentDofValues(dof,0) = value;
            }
            else
            {
                rActiveDofValues(dof,0) = value;
            }
        }
    }
}

void NuTo::StructureIp::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
{
    // initialize vectors
    assert(rActiveDofGradientVector.GetNumRows() == this->mNumActiveDofs);
    assert(rActiveDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rActiveDofGradientVector.GetNumRows(); row ++)
    {
        rActiveDofGradientVector(row,0) = 0.0;
    }
    assert(rDependentDofGradientVector.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rDependentDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rDependentDofGradientVector.GetNumRows(); row ++)
    {
        rDependentDofGradientVector(row,0) = 0.0;
    }

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,6> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);
    double bAlphaRow, bURow[2], bMRow[3];

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        elementIter->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
        assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());
        assert(static_cast<unsigned int>(elementVector.GetNumColumns()) == 1);

        // write element contribution to global vectors
        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
        {
            int globalRowDof = elementVectorGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                rActiveDofGradientVector(globalRowDof,0) += elementVector(rowCount,0);
            }
            else
            {
                assert(globalRowDof-this->mNumActiveDofs < this->mNumDofs - this->mNumActiveDofs);
                rDependentDofGradientVector(globalRowDof - this->mNumActiveDofs,0) += elementVector(rowCount,0);
            }
            //add contribution of the global degrees of freedom (crack opening and crack orientation)
            if (mappingDofMultiscaleNode[globalRowDof]!=-1)
            {
                //influence of alpha
                int theDofMapRow = mappingDofMultiscaleNode[globalRowDof];
                bAlphaRow =dDOF[theDofMapRow][0];

                //influence of u
                bURow[0] = dDOF[theDofMapRow][1];
                bURow[1] = dDOF[theDofMapRow][2];

                //influence of epsilon_tot = epsilon_hom
                bMRow[0] = dDOF[theDofMapRow][3];
                bMRow[1] = dDOF[theDofMapRow][4];
                bMRow[2] = dDOF[theDofMapRow][5];

                if (mDOFCrackAngle< this->mNumActiveDofs)
                    rActiveDofGradientVector(mDOFCrackAngle,0) += bAlphaRow * elementVector(rowCount,0);
                else
                    rDependentDofGradientVector(mDOFCrackAngle - this->mNumActiveDofs,0) += bAlphaRow * elementVector(rowCount,0);

                if (mDOFCrackOpening[0]< this->mNumActiveDofs)
                    rActiveDofGradientVector(mDOFCrackOpening[0],0) += bURow[0] * elementVector(rowCount,0);
                else
                    rDependentDofGradientVector(mDOFCrackOpening[0] - this->mNumActiveDofs,0) += bURow[0] * elementVector(rowCount,0);

                if (mDOFCrackOpening[1]< this->mNumActiveDofs)
                    rActiveDofGradientVector(mDOFCrackOpening[1],0) += bURow[1] * elementVector(rowCount,0);
                else
                    rDependentDofGradientVector(mDOFCrackOpening[1] - this->mNumActiveDofs,0) += bURow[1] * elementVector(rowCount,0);

                if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                {
                    assert(mDOFGlobalTotalStrain[1]< this->mNumActiveDofs);
                    assert(mDOFGlobalTotalStrain[2]< this->mNumActiveDofs);
                    rActiveDofGradientVector(mDOFGlobalTotalStrain[0],0) += bMRow[0] * elementVector(rowCount,0);
                    rActiveDofGradientVector(mDOFGlobalTotalStrain[1],0) += bMRow[1] * elementVector(rowCount,0);
                    rActiveDofGradientVector(mDOFGlobalTotalStrain[2],0) += bMRow[2] * elementVector(rowCount,0);
                }
                else
                {
                    assert(mDOFGlobalTotalStrain[1]>= this->mNumActiveDofs);
                    assert(mDOFGlobalTotalStrain[2]>= this->mNumActiveDofs);
                    rDependentDofGradientVector(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs,0) += bMRow[0] * elementVector(rowCount,0);
                    rDependentDofGradientVector(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs,0) += bMRow[1] * elementVector(rowCount,0);
                    rDependentDofGradientVector(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs,0) += bMRow[2] * elementVector(rowCount,0);

                }
            }
        }

        elementIter++;
    }

    //write contribution of Lagrange Multipliers
    ConstraintBuildGlobalGradientInternalPotentialSubVectors(rActiveDofGradientVector, rDependentDofGradientVector);
}




// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
//    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementMatrix;
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,6> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        elementIter->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);

        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

        this->AddElementMatrixToGlobalSubMatricesGeneral(
                elementMatrix,
                elementMatrixGlobalDofsRow,
                elementMatrixGlobalDofsColumn,
                elementVector,
                mappingDofMultiscaleNode,
                dDOF,
                dDOF2,
                &rMatrixJJ,
                &rMatrixJK,
                0,
                0,
                false);

        elementIter++;
    }

/*
    {
        std::cout << "JJ before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
        rMatrixJJFull.Info(12,5);
        std::cout << "JK before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
        rMatrixJKFull.Info(12,5);

    }
*/
    //add contribution of Lagrange Multipliers
    ConstraintsBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK);
/*
    {
        std::cout << "JJ after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
        rMatrixJJFull.Info(12,5);
        std::cout << "JK after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
        rMatrixJKFull.Info(12,5);

    }
*/
}

// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    //assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);
    assert(rMatrixKJ.IsSymmetric() == false);
    assert(rMatrixKJ.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixKJ.GetNumEntries() == 0);
    assert(rMatrixKK.IsSymmetric() == false);
    assert(rMatrixKK.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumEntries() == 0);

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementMatrix;
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,6> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        elementIter->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);

        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

        this->AddElementMatrixToGlobalSubMatricesGeneral(
                elementMatrix,
                elementMatrixGlobalDofsRow,
                elementMatrixGlobalDofsColumn,
                elementVector,
                mappingDofMultiscaleNode,
                dDOF,
                dDOF2,
                &rMatrixJJ,
                &rMatrixJK,
                &rMatrixKJ,
                &rMatrixKK,
                true);

        elementIter++;
    }

    /*    {
        std::cout << "JJ before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
        rMatrixJJFull.Info(12,5);
        std::cout << "JK before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
        rMatrixJKFull.Info(12,5);
        std::cout << "JJ before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixKJFull(rMatrixKJ);
        rMatrixKJFull.Info(12,5);
        std::cout << "KK before "<< std::endl;
        NuTo::FullMatrix<double> rMatrixKKFull(rMatrixKK);
        rMatrixKKFull.Info(12,5);

    }
*/
    //add contribution of Lagrange Multipliers
    ConstraintBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK, rMatrixKJ, rMatrixKK);
/*    {
        std::cout << "JJ after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
        rMatrixJJFull.Info(12,5);
        std::cout << "JK after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
        rMatrixJKFull.Info(12,5);
        std::cout << "JK after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixKJFull(rMatrixKJ);
        rMatrixKJFull.Info(12,5);
        std::cout << "KK after "<< std::endl;
        NuTo::FullMatrix<double> rMatrixKKFull(rMatrixKK);
        rMatrixKKFull.Info(12,5);

    }
*/
}

void NuTo::StructureIp::AddElementMatrixToGlobalSubMatricesGeneral(
        NuTo::FullMatrix<double>& elementMatrix,
        std::vector<int>& elementMatrixGlobalDofsRow,
        std::vector<int>& elementMatrixGlobalDofsColumn,
        NuTo::FullMatrix<double>& elementVector,
        std::vector<int>& mappingDofMultiscaleNode,
        std::vector<std::array<double,6> >&dDOF,
        std::vector<std::array<double,3> >&dDOF2,
        SparseMatrix<double>* rMatrixJJ,
        SparseMatrix<double>* rMatrixJK,
        SparseMatrix<double>* rMatrixKJ,
        SparseMatrix<double>* rMatrixKK,
        bool rCalcMatrixKJ_KK)const
{
    assert(rMatrixJJ!=0);
    assert(rMatrixJK!=0);
    //either both are zero or both are nonzero
    assert(!rCalcMatrixKJ_KK || ((rMatrixKJ!=0) && (rMatrixKK!=0)));
    double bAlphaRow, bURow[2], bMRow[3], bAlphaCol, bUCol[2], bMCol[3];

    // write element contribution to global matrix
    for (unsigned int rowCount = 0; rowCount < elementMatrixGlobalDofsRow.size(); rowCount++)
    {
        int globalRowDof = elementMatrixGlobalDofsRow[rowCount];

        if (globalRowDof < this->mNumActiveDofs)
        {
            for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
            {
                int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
                if (globalColumnDof < this->mNumActiveDofs)
                {
                    rMatrixJJ->AddEntry(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));

                }
                else
                {
                    rMatrixJK->AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                }
            }
        }
        else
        {
            if (rCalcMatrixKJ_KK)
            {
                for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
                {
                    int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
                    if (globalColumnDof < this->mNumActiveDofs)
                    {
                         rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                         rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                    }
                }
            }
        }

        // add influence of global variables
        //include influence of global dofs like crack orientation (alpha) and crackopening (ux,uy)
        if (mappingDofMultiscaleNode[globalRowDof]!=-1)
        {
            //influence of alpha
            int theDofMapRow = mappingDofMultiscaleNode[globalRowDof];
            bAlphaRow =dDOF[theDofMapRow][0];

            //influence of u
            bURow[0] = dDOF[theDofMapRow][1];
            bURow[1] = dDOF[theDofMapRow][2];

            //influence of e_hom = e_tot
            bMRow[0] = dDOF[theDofMapRow][3];
            bMRow[1] = dDOF[theDofMapRow][4];
            bMRow[2] = dDOF[theDofMapRow][5];

            //add the influence of the second order derivatives
            if (mDOFCrackAngle < this->mNumActiveDofs)
            {
                rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                {
                    rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                }
                else
                {
                    rMatrixJK->AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    if (rCalcMatrixKJ_KK)
                        rMatrixKJ->AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                }
                if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                {
                    rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                }
                else
                {
                    rMatrixJK->AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    if (rCalcMatrixKJ_KK)
                        rMatrixKJ->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                }
            }
            else
            {
                if (rCalcMatrixKJ_KK)
                    rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                {
                    if (rCalcMatrixKJ_KK)
                        rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                }
                else
                {
                    if (rCalcMatrixKJ_KK)
                    {
                        rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixKK->AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                }
                if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                {
                    if (rCalcMatrixKJ_KK)
                        rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                }
                else
                {
                    if (rCalcMatrixKJ_KK)
                    {
                        rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixKK->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                }
            }
        }

        for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
        {
            int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
            if (mappingDofMultiscaleNode[globalColumnDof]!=-1)
            {
                //include influence of global dofs like crack orientation (alpha) and crackopening (ux,uy)
                //influence of alpha
                int theDofMapCol = mappingDofMultiscaleNode[globalColumnDof];
                bAlphaCol = dDOF[theDofMapCol][0];

                //influence of u
                bUCol[0] = dDOF[theDofMapCol][1];
                bUCol[1] = dDOF[theDofMapCol][2];

                //influence of e_hom and e_tot
                bMCol[0] = dDOF[theDofMapCol][3];
                bMCol[1] = dDOF[theDofMapCol][4];
                bMCol[2] = dDOF[theDofMapCol][5];

                if (globalRowDof < this->mNumActiveDofs)
                {
                    if (mDOFCrackAngle < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(globalRowDof, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                    }
                    else
                    {
                        rMatrixJK->AddEntry(globalRowDof, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                    }

                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(globalRowDof, mDOFCrackOpening[0], elementMatrix(rowCount, colCount) * bUCol[0]);
                    }
                    else
                    {
                        rMatrixJK->AddEntry(globalRowDof, mDOFCrackOpening[0] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[0]);
                    }

                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(globalRowDof, mDOFCrackOpening[1], elementMatrix(rowCount, colCount) * bUCol[1]);
                    }
                    else
                    {
                        rMatrixJK->AddEntry(globalRowDof, mDOFCrackOpening[1] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[1]);
                    }
                    if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                    {
                        assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                        rMatrixJJ->AddEntry(globalRowDof, mDOFGlobalTotalStrain[0],  elementMatrix(rowCount, colCount) * bMCol[0]);
                        rMatrixJJ->AddEntry(globalRowDof, mDOFGlobalTotalStrain[1],  elementMatrix(rowCount, colCount) * bMCol[1]);
                        rMatrixJJ->AddEntry(globalRowDof, mDOFGlobalTotalStrain[2],  elementMatrix(rowCount, colCount) * bMCol[2]);
                    }
                    else
                    {
                        assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                        rMatrixJK->AddEntry(globalRowDof, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[0]);
                        rMatrixJK->AddEntry(globalRowDof, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[1]);
                        rMatrixJK->AddEntry(globalRowDof, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[2]);
                    }

                }
                else
                {
                    if (rCalcMatrixKJ_KK)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[0], elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[0]);
                        }

                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[1], elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                        {
                            assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[0], elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[1], elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[2], elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                        else
                        {
                            assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                    }
                }
                if (mappingDofMultiscaleNode[globalRowDof]!=-1)
                {
                    if (mDOFCrackAngle < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(mDOFCrackAngle,      mDOFCrackAngle,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }

                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                        {
                            assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                            rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[0], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[1], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[2], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                        else
                        {
                            assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                            rMatrixJK->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJK->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJK->AddEntry(mDOFCrackAngle, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                        {
                            rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs,      mDOFCrackAngle - this->mNumActiveDofs,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }

                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                            {
                                assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[0], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[1], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[2], bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }
                            else
                            {
                                assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }

                        }
                    }

                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[0], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1]-this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                        {
                            assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[0], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[1], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[2], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                        else
                        {
                            assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                            rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[0] -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[1], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[1]-this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                            {
                                assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[0], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[1], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[2], bURow[0] * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }
                            else
                            {
                                assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                                rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bURow[0] * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }
                        }
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[1], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                        {
                            assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[0], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[1], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[2], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                        else
                        {
                            assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                            assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                            rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                            }
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackOpening[0], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFCrackOpening[1]- this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            rMatrixKK->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                            {
                                assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[0], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[1], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKJ->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[2], bURow[1] * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }
                            else
                            {
                                assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                                assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                                rMatrixKK->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[0]);
                                rMatrixKK->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[1]);
                                rMatrixKK->AddEntry(mDOFCrackOpening[1] -this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bURow[1] * elementMatrix(rowCount, colCount) * bMCol[2]);
                            }

                        }
                    }
                    if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                    {
                       assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                       assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                       if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackAngle, bMRow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackAngle, bMRow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackAngle, bMRow[2]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackAngle - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackAngle - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackAngle - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackOpening[0], bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackOpening[0], bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackOpening[0], bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackOpening[1], bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackOpening[1], bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackOpening[1], bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[0], mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[1], mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            rMatrixJK->AddEntry(mDOFGlobalTotalStrain[2], mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFGlobalTotalStrain[0], bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFGlobalTotalStrain[1], bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], mDOFGlobalTotalStrain[2], bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFGlobalTotalStrain[0], bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFGlobalTotalStrain[1], bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], mDOFGlobalTotalStrain[2], bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFGlobalTotalStrain[0], bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFGlobalTotalStrain[1], bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], mDOFGlobalTotalStrain[2], bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                    }
                    else
                    {
                       assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                       assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                       if (rCalcMatrixKJ_KK)
                       {
                           if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackAngle, bMRow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackAngle, bMRow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackAngle, bMRow[2]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackOpening[0], bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackOpening[0], bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackOpening[0], bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackOpening[1], bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackOpening[1], bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                                rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackOpening[1], bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                                rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFGlobalTotalStrain[0]- this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[0]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFGlobalTotalStrain[1]- this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[1]);
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, mDOFGlobalTotalStrain[2]- this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount) * bMCol[2]);
                        }
                    }
                }
            }
            if (mappingDofMultiscaleNode[globalRowDof]!=-1)
            {
                if (globalColumnDof < this->mNumActiveDofs)
                {
                    if (mDOFCrackAngle < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(mDOFCrackAngle,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(mDOFCrackOpening[0], globalColumnDof, bURow[0]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKJ->AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, globalColumnDof, bURow[0]  * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixJJ->AddEntry(mDOFCrackOpening[1], globalColumnDof, bURow[1]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKJ->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, globalColumnDof, bURow[1]  * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                    {
                        assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[0], globalColumnDof, bMRow[0]  * elementMatrix(rowCount, colCount));
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[1], globalColumnDof, bMRow[1]  * elementMatrix(rowCount, colCount));
                        rMatrixJJ->AddEntry(mDOFGlobalTotalStrain[2], globalColumnDof, bMRow[2]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                        if (rCalcMatrixKJ_KK)
                        {
                            rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[0] - this->mNumActiveDofs, globalColumnDof, bMRow[0]  * elementMatrix(rowCount, colCount));
                            rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[1] - this->mNumActiveDofs, globalColumnDof, bMRow[1]  * elementMatrix(rowCount, colCount));
                            rMatrixKJ->AddEntry(mDOFGlobalTotalStrain[2] - this->mNumActiveDofs, globalColumnDof, bMRow[2]  * elementMatrix(rowCount, colCount));
                        }
                    }
                }
                else
                {
                    if (mDOFCrackAngle < this->mNumActiveDofs)
                    {
                        rMatrixJK->AddEntry(mDOFCrackAngle, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixJK->AddEntry(mDOFCrackOpening[0], globalColumnDof - this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKK->AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixJK->AddEntry(mDOFCrackOpening[1], globalColumnDof - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        if (rCalcMatrixKJ_KK)
                            rMatrixKK->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount));
                    }
                    if (mDOFGlobalTotalStrain[0] < this->mNumActiveDofs)
                    {
                        assert(mDOFGlobalTotalStrain[1]<this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]<this->mNumActiveDofs);
                        rMatrixJK->AddEntry(mDOFGlobalTotalStrain[0], globalColumnDof - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount));
                        rMatrixJK->AddEntry(mDOFGlobalTotalStrain[1], globalColumnDof - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount));
                        rMatrixJK->AddEntry(mDOFGlobalTotalStrain[2], globalColumnDof - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                        assert(mDOFGlobalTotalStrain[1]>=this->mNumActiveDofs);
                        assert(mDOFGlobalTotalStrain[2]>=this->mNumActiveDofs);
                        if (rCalcMatrixKJ_KK)
                        {
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[0] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bMRow[0]  * elementMatrix(rowCount, colCount));
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[1] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bMRow[1]  * elementMatrix(rowCount, colCount));
                            rMatrixKK->AddEntry(mDOFGlobalTotalStrain[2] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bMRow[2]  * elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
    }
}

//! @brief calculate the derivative of the displacements at the nodes with respect to crack opening and crack orientation
//! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
//! @param rDOF return value, for each dof, the corresponding derivatives (alpha, ux, uy, exx, eyy, gamma_xy)
//! @param r2DOF return value, for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy)
void NuTo::StructureIp::CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<std::array<double,6> >& rDOF, std::vector<std::array<double,3> >& rDOF2)const
{
    rMappingDofMultiscaleNode.resize(mNumDofs,-1);

    //first loop, just count the DOFS
    int numMultiscaleDofs(0);
    for (boost::ptr_map<int,NodeBase>::const_iterator itNode=mNodeMap.begin(); itNode!=mNodeMap.end(); itNode++)
    {
        if (itNode->second->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
        {
            numMultiscaleDofs+=2;
        }
        else
        {
            if (itNode->second->GetNumFineScaleDisplacements()!=0)
                throw MechanicsException("[NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0General] This multiscale node type has not been implemented.");
        }
    }

    rDOF.resize(numMultiscaleDofs);
    rDOF2.resize(numMultiscaleDofs);
    int countMultiscaleDofs(0);
    double bHomAlpha[3], bHomU[6]; //bHomU[0-2] for ux [3-5] for uy
    double bHessian[9]; //depsilondalpha2[0-2], depsilondalphadux[3-5], depsilondalphadux[6-8]
    GetdEpsilonHomdCrack(bHomAlpha,bHomU,bHessian);

    //check the derivatives
#ifdef MYDEBUG
    double interval(1e-8);
    EngineeringStrain2D strain1 = mEpsilonHom;
    std::cout << "algo alpha " << bHomAlpha[0] << "  " << bHomAlpha[1] << "  " <<  bHomAlpha[2] << std::endl;
    const_cast<StructureIp*>(&*this)->mCrackAngle += interval;
    const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
    EngineeringStrain2D strain2 = mEpsilonHom;
    std::cout << "cdf alpha " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
                        << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
                        << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< std::endl;
    const_cast<StructureIp*>(&*this)->mCrackAngle -= interval;
    const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();

    for (int count=0; count<2; count++)
    {
        std::cout << "algo u" << bHomU[0+3*count] << "  " << bHomU[1+3*count] << "  " <<  bHomU[2+3*count] << std::endl;
        strain1 = mEpsilonHom;
        const_cast<StructureIp*>(&*this)->mCrackOpening[count]+=interval;
        const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
        EngineeringStrain2D strain2 = mEpsilonHom;
        std::cout << "cdf " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
                            << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
                            << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< std::endl;
        const_cast<StructureIp*>(&*this)->mCrackOpening[count] -= interval;
        const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
    }

    //check hessian
    double hessianCDF[9],bHomAlpha2[3],bHomU2[6];
    std::cout << "algo hessian " << std::endl;
    std::cout << bHessian[0] << "  " << bHessian[3] << "  " <<  bHessian[6] << std::endl;
    std::cout << bHessian[1] << "  " << bHessian[4] << "  " <<  bHessian[7] << std::endl;
    std::cout << bHessian[2] << "  " << bHessian[5] << "  " <<  bHessian[8] << std::endl;
    const_cast<StructureIp*>(&*this)->mCrackAngle += interval;
    const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
    GetdEpsilonHomdCrack(bHomAlpha2,bHomU2,0);
    hessianCDF[0] = (bHomAlpha2[0]-bHomAlpha[0])/interval;
    hessianCDF[1] = (bHomAlpha2[1]-bHomAlpha[1])/interval;
    hessianCDF[2] = (bHomAlpha2[2]-bHomAlpha[2])/interval;

    const_cast<StructureIp*>(&*this)->mCrackAngle -= interval;
    const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();

    for (int count=0; count<2; count++)
    {
        const_cast<StructureIp*>(&*this)->mCrackOpening[count]+=interval;
        const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
        GetdEpsilonHomdCrack(bHomAlpha2,bHomU2,0);
        hessianCDF[3+3*count] = (bHomAlpha2[0]-bHomAlpha[0])/interval;
        hessianCDF[4+3*count] = (bHomAlpha2[1]-bHomAlpha[1])/interval;
        hessianCDF[5+3*count] = (bHomAlpha2[2]-bHomAlpha[2])/interval;

        const_cast<StructureIp*>(&*this)->mCrackOpening[count] -= interval;
        const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
    }
    std::cout << "cdf hessian " << std::endl;
    std::cout << hessianCDF[0] << "  " << hessianCDF[3] << "  " <<  hessianCDF[6] << std::endl;
    std::cout << hessianCDF[1] << "  " << hessianCDF[4] << "  " <<  hessianCDF[7] << std::endl;
    std::cout << hessianCDF[2] << "  " << hessianCDF[5] << "  " <<  hessianCDF[8] << std::endl;

    std::cout << std::endl;
#endif

    //second loop, calculate derivates
    for (boost::ptr_map<int,NodeBase>::const_iterator itNode=mNodeMap.begin(); itNode!=mNodeMap.end(); itNode++)
    {
        if (itNode->second->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
        {
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(0)] = countMultiscaleDofs;
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(1)] = countMultiscaleDofs+1;

            double coord[2];
            itNode->second->GetCoordinates2D(coord);
            //derivative of displacement with respect to homogeneous strain
            double dDOFdEpsilonHomX[3];
            double dDOFdEpsilonHomY[3];
            GetdDisplacementdEpsilonHom(coord,dDOFdEpsilonHomX,dDOFdEpsilonHomY);
#ifdef MYDEBUG
            {
            double interval(1e-3);
            double disp1[2],disp2[2];
            double dDOFdEpsilonHomXcdf[3];
            double dDOFdEpsilonHomYcdf[3];
            std::cout << "dDofxdEpsilonHom algo " << dDOFdEpsilonHomX[0] << "  " << dDOFdEpsilonHomX[1] << "  " << dDOFdEpsilonHomX[2] << std::endl;
            std::cout << "dDofydEpsilonHom algo " << dDOFdEpsilonHomY[0] << "  " << dDOFdEpsilonHomY[1] << "  " << dDOFdEpsilonHomY[2] << std::endl;
            GetDisplacementsEpsilonHom2D(coord, disp1);
            for (int count=0; count<3; count++)
            {
                const_cast<StructureIp*>(&*this)->mEpsilonHom.mEngineeringStrain[count] += interval;
                GetDisplacementsEpsilonHom2D(coord, disp2);
                dDOFdEpsilonHomXcdf[count] = (disp2[0]-disp1[0])/interval;
                dDOFdEpsilonHomYcdf[count] = (disp2[1]-disp1[1])/interval;
                const_cast<StructureIp*>(&*this)->mEpsilonHom.mEngineeringStrain[count] -= interval;
            }
            std::cout << "dDofxdEpsilonHom cdf " << dDOFdEpsilonHomXcdf[0] << "  " << dDOFdEpsilonHomXcdf[1] << "  " << dDOFdEpsilonHomXcdf[2] << std::endl;
            std::cout << "dDofydEpsilonHom cdf " << dDOFdEpsilonHomYcdf[0] << "  " << dDOFdEpsilonHomYcdf[1] << "  " << dDOFdEpsilonHomYcdf[2] << std::endl;
            std::cout << std::endl;
            }
#endif

            //derivative of displacement with respect to discontinuity (crack orientation)
            GetdDisplacementdCrackOrientation(coord,&(rDOF[countMultiscaleDofs][0]),&(rDOF[countMultiscaleDofs+1][0]));

            rDOF[countMultiscaleDofs][0]  +=  dDOFdEpsilonHomX[0]*bHomAlpha[0] +
                                              dDOFdEpsilonHomX[1]*bHomAlpha[1] +
                                              dDOFdEpsilonHomX[2]*bHomAlpha[2];
            rDOF[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHomAlpha[0] +
                                              dDOFdEpsilonHomY[1]*bHomAlpha[1] +
                                              dDOFdEpsilonHomY[2]*bHomAlpha[2];
#ifdef MYDEBUG
            {
            double interval(1e-5);
            double disp1[2],disp2[2],disptmp[2];
            std::cout << "dDofdAlpha algo " << rDOF[countMultiscaleDofs][0] << "  " << rDOF[countMultiscaleDofs+1][0] << std::endl;
            GetDisplacementsCrack2D(coord, disp1);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp1[0] += disptmp[0];
            disp1[1] += disptmp[1];
            const_cast<StructureIp*>(&*this)->mCrackAngle += interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            GetDisplacementsCrack2D(coord, disp2);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp2[0] += disptmp[0];
            disp2[1] += disptmp[1];
            std::cout << "dDofdAlpha cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][0])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][0])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureIp*>(&*this)->mCrackAngle -= interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            }
#endif
            //derivative of displacement with respect to discontinuity (crack opening)
            GetdDisplacementdCrackOpening(coord,&(rDOF[countMultiscaleDofs][1]),&(rDOF[countMultiscaleDofs+1][1]));
            rDOF[countMultiscaleDofs][1]   += dDOFdEpsilonHomX[0]*bHomU[0] +
                                              dDOFdEpsilonHomX[1]*bHomU[1] +
                                              dDOFdEpsilonHomX[2]*bHomU[2];
            rDOF[countMultiscaleDofs+1][1] += dDOFdEpsilonHomY[0]*bHomU[0] +
                                              dDOFdEpsilonHomY[1]*bHomU[1] +
                                              dDOFdEpsilonHomY[2]*bHomU[2];
#ifdef MYDEBUG
            {
            double interval(1e-5);
            double disp1[2],disp2[2],disptmp[2];
            std::cout << "dDofdU1 algo " << rDOF[countMultiscaleDofs][1] << "  " << rDOF[countMultiscaleDofs+1][1] << std::endl;
            GetDisplacementsCrack2D(coord, disp1);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp1[0] += disptmp[0];
            disp1[1] += disptmp[1];
            const_cast<StructureIp*>(&*this)->mCrackOpening[0] += interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            GetDisplacementsCrack2D(coord, disp2);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp2[0] += disptmp[0];
            disp2[1] += disptmp[1];
            std::cout << "dDofdU1 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][1])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][1])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureIp*>(&*this)->mCrackOpening[0] -= interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            }
#endif
            rDOF[countMultiscaleDofs][2]   += dDOFdEpsilonHomX[0]*bHomU[3] +
                                              dDOFdEpsilonHomX[1]*bHomU[4] +
                                              dDOFdEpsilonHomX[2]*bHomU[5];
            rDOF[countMultiscaleDofs+1][2] += dDOFdEpsilonHomY[0]*bHomU[3] +
                                              dDOFdEpsilonHomY[1]*bHomU[4] +
                                              dDOFdEpsilonHomY[2]*bHomU[5];
#ifdef MYDEBUG
            {
            double interval(1e-5);
            double disp1[2],disp2[2],disptmp[2];
            std::cout << "dDofdU2 algo " << rDOF[countMultiscaleDofs][2] << "  " << rDOF[countMultiscaleDofs+1][2] << std::endl;
            GetDisplacementsCrack2D(coord, disp1);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp1[0] += disptmp[0];
            disp1[1] += disptmp[1];
            const_cast<StructureIp*>(&*this)->mCrackOpening[1] += interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            GetDisplacementsCrack2D(coord, disp2);
            GetDisplacementsEpsilonHom2D(coord, disptmp);
            disp2[0] += disptmp[0];
            disp2[1] += disptmp[1];
            std::cout << "dDofdU2 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][2])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][2])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureIp*>(&*this)->mCrackOpening[1] -= interval;
            const_cast<StructureIp*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            }
#endif
            //derivative of displacement with respect to homogeneous strain (is equal to derivative w.r.t. total strain)
            rDOF[countMultiscaleDofs][3]   = dDOFdEpsilonHomX[0];
            rDOF[countMultiscaleDofs][4]   = dDOFdEpsilonHomX[1];
            rDOF[countMultiscaleDofs][5]   = dDOFdEpsilonHomX[2];

            rDOF[countMultiscaleDofs+1][3]   = dDOFdEpsilonHomY[0];
            rDOF[countMultiscaleDofs+1][4]   = dDOFdEpsilonHomY[1];
            rDOF[countMultiscaleDofs+1][5]   = dDOFdEpsilonHomY[2];

            //calculate second order derivatives
            //second derivative of displacement with respect to discontinuity (crack orientation)
            Getd2Displacementd2CrackOrientation(coord,&(rDOF2[countMultiscaleDofs][0]),&(rDOF2[countMultiscaleDofs+1][0]));

            rDOF2[countMultiscaleDofs][0]   += dDOFdEpsilonHomX[0]*bHessian[0] +
                                               dDOFdEpsilonHomX[1]*bHessian[1] +
                                               dDOFdEpsilonHomX[2]*bHessian[2];
            rDOF2[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHessian[0] +
                                               dDOFdEpsilonHomY[1]*bHessian[1] +
                                               dDOFdEpsilonHomY[2]*bHessian[2];

            //derivative of displacement with respect to alpha and discontinuity (crack opening)
            Getd2Displacementd2CrackOpening(coord, &(rDOF2[countMultiscaleDofs][1]), &(rDOF2[countMultiscaleDofs+1][1]));
            rDOF2[countMultiscaleDofs][1]   += dDOFdEpsilonHomX[0]*bHessian[3] +
                                               dDOFdEpsilonHomX[1]*bHessian[4] +
                                               dDOFdEpsilonHomX[2]*bHessian[5];
            rDOF2[countMultiscaleDofs][2]   += dDOFdEpsilonHomX[0]*bHessian[6] +
                                               dDOFdEpsilonHomX[1]*bHessian[7] +
                                               dDOFdEpsilonHomX[2]*bHessian[8];
            rDOF2[countMultiscaleDofs+1][1] += dDOFdEpsilonHomY[0]*bHessian[3] +
                                               dDOFdEpsilonHomY[1]*bHessian[4] +
                                               dDOFdEpsilonHomY[2]*bHessian[5];
            rDOF2[countMultiscaleDofs+1][2] += dDOFdEpsilonHomY[0]*bHessian[6] +
                                               dDOFdEpsilonHomY[1]*bHessian[7] +
                                               dDOFdEpsilonHomY[2]*bHessian[8];

            countMultiscaleDofs+=2;
            if (countMultiscaleDofs>numMultiscaleDofs)
            {
                throw MechanicsException("[NuTo::StructureIp::CalculatedDispdGlobalDofs] countMultiscaleDofs>numMultiscaleDofs internal error.");
            }
        }
    }
#ifdef MYDEBUG
    std::cout << "End of routine" << std::endl << std::endl;
#endif

    if (countMultiscaleDofs!=numMultiscaleDofs)
        throw MechanicsException("[NuTo::StructureIp::CalculatedDispdGlobalDofs] countMultiscaleDofs!=numMultiscaleDofs internal error.");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
void NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const
{
    throw MechanicsException("[NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0Symmetric] not implemented for StructureIp");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
//! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
void NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const
{
    throw MechanicsException("[NuTo::StructureIp::BuildGlobalCoefficientSubMatrices0Symmetric] not implemented for StructureIp");
}

//! @brief ... calculate the displacement based on the homogeneous strain
//! @param rCoordinates ... coordinates of the point
//! @param rCoarseDisplacements ... return array of displacements
void NuTo::StructureIp::GetDisplacementsEpsilonHom2D(double rCoordinates[2], double rDisplacements[2])const
{
    rDisplacements[0] = mEpsilonHom.mEngineeringStrain[0] * (rCoordinates[0]-mlX*0.5) + 0.5 * mEpsilonHom.mEngineeringStrain[2] * (rCoordinates[1]-mlX*0.5);
    rDisplacements[1] = mEpsilonHom.mEngineeringStrain[1] * (rCoordinates[1]-mlX*0.5) + 0.5 * mEpsilonHom.mEngineeringStrain[2] * (rCoordinates[0]-mlX*0.5);
}

//! @brief ... calculate the displacement based on the crack opening
//! @param rCoordinates ... coordinates of the point
//! @param rCoarseDisplacements ... return array of displacements
void NuTo::StructureIp::GetDisplacementsCrack2D(double rCoordinates[2], double rDisplacements[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    double factor;
    if (d<-mCrackTransitionZone)
    {
        factor = -0.5;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            factor = 0.5;
        }
        else
        {
            //smooth transition from cracking to none cracking
            factor = 0.5*sin(0.5*M_PI*d/mCrackTransitionZone);
        }
    }

    rDisplacements[0] = factor*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
    rDisplacements[1] = factor*(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1]);

}

//! @brief derivative of displacement with respect to homogeneous strain
//! @param rdX_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
//! @param rdY_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
void NuTo::StructureIp::GetdDisplacementdEpsilonHom(double rCoordinates[2], double rdX_dEpsilonHom[3], double rdY_dEpsilonHom[3])const
{
    rdX_dEpsilonHom[0] = rCoordinates[0]-mlX*0.5;
    rdX_dEpsilonHom[1] = 0.;
    rdX_dEpsilonHom[2] = 0.5 * (rCoordinates[1]-mlX*0.5);

    rdY_dEpsilonHom[0] = 0.;
    rdY_dEpsilonHom[1] = rCoordinates[1];
    rdY_dEpsilonHom[2] = 0.5 * (rCoordinates[0]-mlX*0.5);
}

//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
//! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to crack opening (ux, uy)
void NuTo::StructureIp::GetdDisplacementdCrackOpening(double rCoordinates[2], double rdX_dCrackOpening[2], double rdY_dCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double factor;

    if (d<-mCrackTransitionZone)
    {
        factor = -0.5;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            factor = 0.5;
        }
        else
        {
            //smooth transition from cracking to none cracking
            factor=0.5*sin(0.5*M_PI*d/mCrackTransitionZone);
        }
    }
    rdX_dCrackOpening[0] =  factor*cosAlpha;
    rdX_dCrackOpening[1] = -factor*sinAlpha;
    rdY_dCrackOpening[0] =  factor*sinAlpha;
    rdY_dCrackOpening[1] =  factor*cosAlpha;
}

//! @brief second derivative of displacement with respect to alpha and discontinuity (crack opening)
//! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to alpha and crack opening (ux, uy)
//! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to alpha and crack opening (ux, uy)
void NuTo::StructureIp::Getd2Displacementd2CrackOpening(double rCoordinates[2], double rdX_dAlphaCrackOpening[2], double rdY_dAlphaCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double factor;
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionZone)
    {
        factor=-0.5;
        rdX_dAlphaCrackOpening[0] = -factor*sinAlpha;
        rdX_dAlphaCrackOpening[1] = -factor*cosAlpha;
        rdY_dAlphaCrackOpening[0] =  factor*cosAlpha;
        rdY_dAlphaCrackOpening[1] = -factor*sinAlpha;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            factor=0.5;
            rdX_dAlphaCrackOpening[0] = -factor*sinAlpha;
            rdX_dAlphaCrackOpening[1] = -factor*cosAlpha;
            rdY_dAlphaCrackOpening[0] =  factor*cosAlpha;
            rdY_dAlphaCrackOpening[1] = -factor*sinAlpha;
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double dFactor(0.25*M_PI/mCrackTransitionZone*cos(0.5*M_PI*d/mCrackTransitionZone)*dDdAlpha);
            factor = 0.5*sin(0.5*M_PI*d/mCrackTransitionZone);

            rdX_dAlphaCrackOpening[0] = -factor*sinAlpha+dFactor*cosAlpha;
            rdX_dAlphaCrackOpening[1] = -factor*cosAlpha-dFactor*sinAlpha;
            rdY_dAlphaCrackOpening[0] =  factor*cosAlpha+dFactor*sinAlpha;
            rdY_dAlphaCrackOpening[1] = -factor*sinAlpha+dFactor*cosAlpha;
        }
    }
}


//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, derivative of y-displacement with respect to crack orientation (alpha)
void NuTo::StructureIp::GetdDisplacementdCrackOrientation(double rCoordinates[2], double rdX_dAlpha[1], double rdY_dAlpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionZone)
    {
        rdX_dAlpha[0] = -0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
        rdY_dAlpha[0] = -0.5*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rdX_dAlpha[0] = 0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
            rdY_dAlpha[0] = 0.5*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double factor(0.5*sin(0.5*M_PI*d/mCrackTransitionZone));
            double dFactor(0.25*M_PI/mCrackTransitionZone*cos(0.5*M_PI*d/mCrackTransitionZone)*dDdAlpha);

            rdX_dAlpha[0] = factor*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1])+
                           dFactor*(cosAlpha*mCrackOpening[0] - sinAlpha*mCrackOpening[1]);
            rdY_dAlpha[0] = factor*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1])+
                           dFactor*(sinAlpha*mCrackOpening[0] + cosAlpha*mCrackOpening[1]);
        }
    }
#ifdef MYDEBUG

    // test using central differences
    double displacements1[2],displacements2[2];
    double interval(1e-8);
    GetDisplacementsCrack2D(rCoordinates, displacements1);
    const_cast<NuTo::StructureIp*> (this)->mCrackAngle+=interval;
    GetDisplacementsCrack2D(rCoordinates, displacements2);
    const_cast<NuTo::StructureIp*> (this)->mCrackAngle-=interval;
    std::cout << "rdX_dAlpha "<< rdX_dAlpha[0] << " " << (displacements2[0]-displacements1[0])/interval << std::endl;
    std::cout << "rdY_dAlpha "<< rdY_dAlpha[0] << " " << (displacements2[1]-displacements1[1])/interval << std::endl;
    std::cout << d << std::endl;
#endif
}
//! @brief second derivative of displacement with respect to orientation of discontinuity (crack angle)
//! @param rdX_dAlpha[2] return value, second derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, second derivative of y-displacement with respect to crack orientation (alpha)
void NuTo::StructureIp::Getd2Displacementd2CrackOrientation(double rCoordinates[2], double rd2X_d2Alpha[1], double rd2Y_d2Alpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionZone)
    {
        rd2X_d2Alpha[0] = -0.5*(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
        rd2Y_d2Alpha[0] = -0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rd2X_d2Alpha[0] = 0.5*(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
            rd2Y_d2Alpha[0] = 0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double d2DdAlpha = Calculated2DistanceToCrack2Dd2Alpha(rCoordinates);
            double factor(0.5*sin(0.5*M_PI*d/mCrackTransitionZone));
            double dFactor(0.25*M_PI/mCrackTransitionZone*cos(0.5*M_PI*d/mCrackTransitionZone)*dDdAlpha);
            double dFactor2(0.25*M_PI/mCrackTransitionZone*(-0.5*M_PI/mCrackTransitionZone*sin(0.5*M_PI*d/mCrackTransitionZone)*dDdAlpha*dDdAlpha
                                      +cos(0.5*M_PI*d/mCrackTransitionZone)*d2DdAlpha));

            rd2X_d2Alpha[0] = (dFactor2-factor)*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1])+2.*dFactor*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
            rd2Y_d2Alpha[0] = (dFactor2-factor)*(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1])+2.*dFactor*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
        }
    }
#ifdef MYDEBUG

    //check dd2d2alpha
    double interval(1e-8);
    double dX_dAlpha1,dY_dAlpha1,dX_dAlpha2,dY_dAlpha2;
    GetdDisplacementdCrackOrientation(rCoordinates, &dX_dAlpha1, &dY_dAlpha1);
    const_cast<NuTo::StructureIp*> (this)->mCrackAngle+=interval;
    GetdDisplacementdCrackOrientation(rCoordinates, &dX_dAlpha2, &dY_dAlpha2);
    const_cast<NuTo::StructureIp*> (this)->mCrackAngle-=interval;

    std::cout << "dd2d2alpha algo " << rd2X_d2Alpha[0] << " " << rd2Y_d2Alpha[0] << std::endl;
    std::cout << "dd2d2alpha cdf " << (dX_dAlpha2-dX_dAlpha1)/interval << " " << (dY_dAlpha2-dY_dAlpha1)/interval << std::endl;

    if (fabs(rd2X_d2Alpha[0]-(dX_dAlpha2-dX_dAlpha1)/interval)>1e-2)
        exit(0);
    if (fabs(rd2Y_d2Alpha[0]-(dY_dAlpha2-dY_dAlpha1)/interval)>1e-2)
        exit(0);
#endif
}

//! @brief calculate the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return distance to crack
double NuTo::StructureIp::CalculateDistanceToCrack2D(double rCoordinates[2])const
{
    return sin(mCrackAngle)*(mlX/2-rCoordinates[0]) - cos(mCrackAngle)*(mlY/2-rCoordinates[1]);
}

//! @brief calculate the derivative of the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return derivative of distance to crack
double NuTo::StructureIp::CalculatedDistanceToCrack2DdAlpha(double rCoordinates[2])const
{
    return cos(mCrackAngle)*(mlX/2-rCoordinates[0]) + sin(mCrackAngle)*(mlY/2-rCoordinates[1]);
}

//! @brief calculate the second derivative of the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return second derivative of distance to crack
double NuTo::StructureIp::Calculated2DistanceToCrack2Dd2Alpha(double rCoordinates[2])const
{
    return -sin(mCrackAngle)*(mlX/2-rCoordinates[0]) + cos(mCrackAngle)*(mlY/2-rCoordinates[1]);
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureIp::CalculateHomogeneousEngineeringStrain()
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double gamma(1./std::max(fabs(sinAlpha),fabs(cosAlpha)));
    assert(mlX==mlY);
    double factor=gamma/(mlX);

    mEpsilonHom.mEngineeringStrain[0] = mEpsilonTot.mEngineeringStrain[0] - factor * sinAlpha *(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[1] = mEpsilonTot.mEngineeringStrain[1] - factor * cosAlpha *(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[2] = mEpsilonTot.mEngineeringStrain[2] - factor * ((2.*cosAlpha*cosAlpha-1.)* mCrackOpening[0] - 2.*sinAlpha*cosAlpha * mCrackOpening[1]);
/*    mEpsilonHom.mEngineeringStrain[0] = 0;
    mEpsilonHom.mEngineeringStrain[1] = 0;
    mEpsilonHom.mEngineeringStrain[2] = 0;

    std::cout << "calculation of e hom  " << std::endl;
*/
    std::cout << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << std::endl;

    /*
//#ifdef MYDEBUG
    std::cout << "mEpsilonTot " << mEpsilonTot.mEngineeringStrain[0] << " " << mEpsilonTot.mEngineeringStrain[1] << " " << mEpsilonTot.mEngineeringStrain[2] << std::endl;
    std::cout << "mCrackOpening " << mCrackOpening[0] << " " << mCrackOpening[1] << " " << factor << std::endl;
    std::cout << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << std::endl;
//#endif
*/
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureIp::SetTotalEngineeringStrain(EngineeringStrain2D& rTotalEngineeringStrain)
{
    mEpsilonTot = rTotalEngineeringStrain;
#ifdef MYDEBUG
    std::cout << "set total strain to " << rTotalEngineeringStrain.mEngineeringStrain[0]
                                 << " " << rTotalEngineeringStrain.mEngineeringStrain[1]
                                 << " " << rTotalEngineeringStrain.mEngineeringStrain[2] <<std::endl;
#endif
   CalculateHomogeneousEngineeringStrain();
}

//! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
void NuTo::StructureIp::SetDeltaTotalEngineeringStrain(EngineeringStrain2D& rDeltaTotalEngineeringStrain)
{
    mDeltaEpsilonTot = rDeltaTotalEngineeringStrain;
}

//! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
void NuTo::StructureIp::SetPrevTotalEngineeringStrain(EngineeringStrain2D& rPrevTotalEngineeringStrain)
{
    mPrevEpsilonTot = rPrevTotalEngineeringStrain;
}

//! @brief returns the total strain
NuTo::EngineeringStrain2D NuTo::StructureIp::GetTotalEngineeringStrain()const
{
    return mEpsilonTot;
}

//! @brief returns the total strain
NuTo::EngineeringStrain2D NuTo::StructureIp::GetHomogeneousEngineeringStrain()const
{
    return mEpsilonHom;
}

//! @brief set the homgeneous strain (just for test purpose)
void NuTo::StructureIp::SetHomogeneousEngineeringStrain(NuTo::EngineeringStrain2D rStrain)
{
    mEpsilonHom = rStrain;
}

//! @brief return the previous crack angle
double NuTo::StructureIp::GetPrevCrackAngle()const
{
    return mPrevCrackAngle;
}

//! @brief sets the previous crack angle
void NuTo::StructureIp::SetPrevCrackAngle(double rPrevCrackAngle)
{
    mPrevCrackAngle = rPrevCrackAngle;
}

//! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack orientation and crack opening
//! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
//! @parameter rbHomAlpha dHom wrt alpha
//! @paramter rbHomU [0-2] wrt ux [3-5] wrt uy
//! @parameter bHessian depsilondalpha2[0-2], depsilondalphadux[3-5], depsilondalphadux[6-8]
void NuTo::StructureIp::GetdEpsilonHomdCrack(double rbHomAlpha[3], double rbHomU[6], double rbHessian[9])const
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    assert(mlX==mlY);

    if (fabs(sinAlpha)>fabs(cosAlpha))
    {
        double signdivl = sinAlpha<0 ? 1./mlX : -1./mlX;
        rbHomAlpha[0] = signdivl*(sinAlpha*mCrackOpening[0] + cosAlpha*mCrackOpening[1]);
        rbHomAlpha[1] = signdivl*(-sinAlpha*mCrackOpening[0] + cosAlpha/(sinAlpha*sinAlpha)*(cosAlpha*cosAlpha-2.)*mCrackOpening[1]);
        rbHomAlpha[2] = signdivl*(cosAlpha/(sinAlpha*sinAlpha)*(2.*cosAlpha*cosAlpha-3.)*mCrackOpening[0] +2.* sinAlpha*mCrackOpening[1]);

        rbHomU[0] = signdivl*(-cosAlpha);
        rbHomU[1] = signdivl*cosAlpha;
        rbHomU[2] = signdivl*(2.*cosAlpha*cosAlpha-1.)/sinAlpha;

        rbHomU[3] = signdivl*(sinAlpha);
        rbHomU[4] = signdivl*(cosAlpha*cosAlpha)/sinAlpha;
        rbHomU[5] = signdivl*(-2.*cosAlpha);

        if (rbHessian!=0)
        {
            rbHessian[0] = signdivl*(cosAlpha*mCrackOpening[0] - sinAlpha*mCrackOpening[1]);;
            rbHessian[1] = signdivl*(-cosAlpha*mCrackOpening[0] + (cosAlpha*cosAlpha*cosAlpha*cosAlpha-cosAlpha*cosAlpha+2.)/(sinAlpha*sinAlpha*sinAlpha)*mCrackOpening[1]);
            rbHessian[2] = signdivl*((2.*cosAlpha*cosAlpha*cosAlpha*cosAlpha-3.*cosAlpha*cosAlpha+3.)/(sinAlpha*sinAlpha*sinAlpha)*mCrackOpening[0] +2.* cosAlpha*mCrackOpening[1]);

            rbHessian[3] = signdivl*sinAlpha;
            rbHessian[4] = -signdivl*sinAlpha;
            rbHessian[5] = signdivl*(cosAlpha/(sinAlpha*sinAlpha)*(2.*cosAlpha*cosAlpha-3.));

            rbHessian[6] = signdivl*cosAlpha;
            rbHessian[7] = signdivl*cosAlpha/(sinAlpha*sinAlpha)*(cosAlpha*cosAlpha-2.);
            rbHessian[8] = signdivl*2.* sinAlpha;

        }
    }
    else
    {
        double signdivl = cosAlpha<0 ? 1./mlX : -1./mlX;
        rbHomAlpha[0] = signdivl*(-cosAlpha*mCrackOpening[0] + sinAlpha/(cosAlpha*cosAlpha)*(cosAlpha*cosAlpha+1)*mCrackOpening[1]);
        rbHomAlpha[1] = signdivl*( cosAlpha*mCrackOpening[0] - sinAlpha*mCrackOpening[1]);
        rbHomAlpha[2] = signdivl*(-sinAlpha/(cosAlpha*cosAlpha)*(2.*cosAlpha*cosAlpha+1)*mCrackOpening[0] -2.* cosAlpha*mCrackOpening[1]);

        rbHomU[0] = signdivl*(-sinAlpha);
        rbHomU[1] = signdivl*sinAlpha;
        rbHomU[2] = signdivl*(2*cosAlpha*cosAlpha-1.)/cosAlpha;

        rbHomU[3] = signdivl*(sinAlpha*sinAlpha)/cosAlpha;
        rbHomU[4] = signdivl*(cosAlpha);
        rbHomU[5] = signdivl*(-2.*sinAlpha);

        if (rbHessian!=0)
        {
            rbHessian[0] = signdivl*(sinAlpha*mCrackOpening[0] + (cosAlpha*cosAlpha*cosAlpha*cosAlpha-cosAlpha*cosAlpha+2)/(sinAlpha*sinAlpha*sinAlpha)*mCrackOpening[1]);
            rbHessian[1] = signdivl*(-sinAlpha*mCrackOpening[0] - cosAlpha*mCrackOpening[1]);
            rbHessian[2] = signdivl*(-(2.*cosAlpha*cosAlpha*cosAlpha*cosAlpha-cosAlpha*cosAlpha+2.)/(cosAlpha*cosAlpha*cosAlpha)*mCrackOpening[0] +2.* sinAlpha*mCrackOpening[1]);

            rbHessian[3] = -signdivl*cosAlpha;
            rbHessian[4] = signdivl*cosAlpha;
            rbHessian[5] = -signdivl*sinAlpha/(cosAlpha*cosAlpha)*(2.*cosAlpha*cosAlpha+1);

            rbHessian[6] = signdivl*(sinAlpha/(cosAlpha*cosAlpha)*(cosAlpha*cosAlpha+1));
            rbHessian[7] = -signdivl*sinAlpha;
            rbHessian[8] = signdivl*(-2.* cosAlpha);
        }
    }
/*
    std::cout << "crack opening "  << mCrackOpening[0] << " " << mCrackOpening[1]  << std::endl;
    std::cout << "dhom dalpha " << rbHomAlpha[0] << " " <<  rbHomAlpha[1] << " " <<  rbHomAlpha[2] << std::endl;
    std::cout << "dhom dut " << rbHomU[0] << " " <<  rbHomU[1] << " " <<  rbHomU[2] << std::endl;
    std::cout << "dhom dun " << rbHomU[3] << " " <<  rbHomU[4] << " " <<  rbHomU[5] << std::endl;
*/
}

//! @brief return the total strain
NuTo::EngineeringStrain2D NuTo::StructureIp::GetTotalStrain()const
{
    return mEpsilonTot;
}

//! @brief renumbers the global dofs in the structure after
void NuTo::StructureIp::ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOFCrackAngle = rMappingInitialToNewOrdering[mDOFCrackAngle];
    mDOFCrackOpening[0] = rMappingInitialToNewOrdering[mDOFCrackOpening[0]];
    mDOFCrackOpening[1] = rMappingInitialToNewOrdering[mDOFCrackOpening[1]];
    mDOFGlobalTotalStrain[0]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[0]];
    mDOFGlobalTotalStrain[1]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[1]];
    mDOFGlobalTotalStrain[2]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[2]];
}

//! @brief calculates the crack angle for elastic solutions (initial value, no scaling with previous crack angle)
//! @return crack angle in the range 0..Pi
double  NuTo::StructureIp::CalculateInitialCrackAngleElastic()const
{
    FullMatrix<double> strain(2,2);
    const double *dataPtr = mEpsilonTot.GetData();
    strain(0,0) = dataPtr[0];
    strain(0,1) = dataPtr[2]*0.5;
    strain(1,0) = strain(0,1);
    strain(1,1) = dataPtr[1];

//    strain.Info(12,5);

    NuTo::FullMatrix<double> eigenVectors;
    NuTo::FullMatrix<double> eigenValues;
    strain.EigenVectorsSymmetric(eigenVectors);
    strain.EigenValuesSymmetric(eigenValues);

    strain.EigenValuesSymmetric(eigenValues);
    /*    std::cout << "eigenvalues of epsilon tot" << std::endl;
    eigenValues.Trans().Info(12,3);
    std::cout << "eigenvectors of epsilon tot" << std::endl;
    eigenVectors.Info(12,3);
*/
    assert(eigenValues(0,0)<=eigenValues(1,0));
    double alpha = atan2(eigenVectors(1,0),eigenVectors(0,0));
    if (alpha<0)
        alpha+=M_PI;
    return alpha;
}

//! @brief calculates the crack angle for elastic solutions
//! @return crack angle in the range 0..Pi
double NuTo::StructureIp::CalculateCrackAngleElastic()const
{
    FullMatrix<double> strain(2,2);
    const double *dataPtr = mEpsilonTot.GetData();
    strain(0,0) = dataPtr[0];
    strain(0,1) = dataPtr[2]*0.5;
    strain(1,0) = strain(0,1);
    strain(1,1) = dataPtr[1];
//    strain.Info(12,5);

    NuTo::FullMatrix<double> eigenVectors;
    NuTo::FullMatrix<double> eigenValues;
    strain.EigenVectorsSymmetric(eigenVectors);
    strain.EigenValuesSymmetric(eigenValues);

    strain.EigenValuesSymmetric(eigenValues);
    /*    std::cout << "eigenvalues of epsilon tot" << std::endl;
    eigenValues.Trans().Info(12,3);
    std::cout << "eigenvectors of epsilon tot" << std::endl;
    eigenVectors.Info(12,3);
*/
    assert(eigenValues(0,0)<=eigenValues(1,0));
    double alpha_2(0);
    if (eigenValues(1,0)-eigenValues(0,0) > mToleranceElasticCrackAngleHigh)
    {
        alpha_2 = atan2(eigenVectors(1,0),eigenVectors(0,0));
        if (alpha_2<0)
            alpha_2+=M_PI;
    }
    else
    {
        if (eigenValues(1,0)-eigenValues(0,0)<mToleranceElasticCrackAngleLow)
        {
            alpha_2 = mPrevCrackAngle;
        }
        else
        {
            double s((eigenValues(1,0)-eigenValues(0,0)-mToleranceElasticCrackAngleLow)/(mToleranceElasticCrackAngleHigh - mToleranceElasticCrackAngleLow));
            alpha_2 = atan2(eigenVectors(1,0),eigenVectors(0,0));
            if (alpha_2<0)
                alpha_2+=M_PI;
            alpha_2 = alpha_2* (1-s) + mPrevCrackAngle*s;
        }
    }
    return alpha_2;
}

//! @brief calculates the difference between the crack angle of the elastic solution and the current angle
//! attention, the periodicity of the crack angle has to be taken into account
double NuTo::StructureIp::CalculateDeltaCrackAngleElastic()const
{
    //calculate angle orthogonal to second principal stress
     double alpha_2 = this->CalculateCrackAngleElastic();

     double delta_alpha = mCrackAngle-alpha_2;
     while (fabs(delta_alpha>M_PI))
     {
         if (delta_alpha>0)
         {
             delta_alpha-=M_PI;
             alpha_2+=M_PI;
         }
         else
         {
             delta_alpha+=M_PI;
             alpha_2-=M_PI;
         }
     }
     //now delta_alpha is in [0..PI]

     if (delta_alpha>0.5*M_PI)
     {
         delta_alpha = M_PI-delta_alpha;
         alpha_2 += M_PI;
     }
     //std::cout << "alpha " << mCrackAngle << " alpha2 " << alpha_2 << " delta " << delta_alpha << std::endl;

     return delta_alpha;
}

//! @brief add a constraint equation for alpha, which corresponds to an artificial spring
//! @parameter rPenaltyStiffness penalty stiffness
//! @parameter rScalingFactor scaling factor
//! @return id of the constraint
int NuTo::StructureIp::ConstraintNonlinearCrackAngle(double rPenaltyStiffness, double rScalingFactor)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintCrackAngle);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 1;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    ConstraintNonlinearGlobalCrackAngle2D *mConst = new NuTo::ConstraintNonlinearGlobalCrackAngle2D(this, rPenaltyStiffness, rScalingFactor);

    mConstraintMap.insert(id, mConst);
    mConstraintCrackAngle = id;
    return id;
}

//! @brief set the penalty stiffness for the nonlinear crack angle constraint
void NuTo::StructureIp::SetPenaltyStiffnessCrackAngle(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintCrackAngle);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackAngle2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackAngle2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessCrackAngle] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetPenaltyStiffness(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessCrackAngle] There is no constraint for the crack angle.");
    }
}

//! @brief set the scaling factor for the nonlinear crack angle constraint
void NuTo::StructureIp::SetPenaltyStiffnessScalingFactorCrackAngle(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintCrackAngle);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackAngle2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackAngle2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessScalingFactorCrackAngle] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetScalingFactor(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessScalingFactorCrackAngle] There is no constraint for the crack angle.");
    }

}

//! @brief add a constraint equation for the tangential crack opening, which corresponds to an artificial spring
//! @parameter rPenaltyStiffness penalty stiffness
//! @parameter rScalingFactor scaling factor
//! @return id of the constraint
int NuTo::StructureIp::ConstraintNonlinearTangentialCrackOpening(double rScalingFactor, double rPenaltyStiffness)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 1;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    ConstraintNonlinearGlobalCrackOpeningTangential2D *mConst = new NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D(this, rScalingFactor, rPenaltyStiffness);

    mConstraintMap.insert(id, mConst);
    mConstraintTangentialCrackOpening = id;
    return id;
}

//! @brief set the penalty stiffness for the nonlinear TangentialCrackOpening
void NuTo::StructureIp::SetPenaltyStiffnessTangentialCrackOpening(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackOpeningTangential2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackOpeningTangential2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessTangentialCrackOpening] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetPenaltyStiffness(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessTangentialCrackOpening] There is no constraint for the crack angle.");
    }
}

//! @brief set the scaling factor for the nonlinear TangentialCrackOpening
void NuTo::StructureIp::SetPenaltyStiffnessScalingFactorTangentialCrackOpening(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackOpeningTangential2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackOpeningTangential2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessScalingFactorTangentialCrackOpening] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetScalingFactor(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureIp::SetPenaltyStiffnessScalingFactorTangentialCrackOpening] There is no constraint for the crack angle.");
    }

}

//! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
void NuTo::StructureIp::SetToleranceElasticCrackAngleHigh(double rParameter)
{
    if (rParameter<mToleranceElasticCrackAngleLow || rParameter<0)
        throw MechanicsException("[NuTo::StructureIp::SetToleranceElasticCrackAngleHigh] Tolerance has to be positive and higher than the low tolerance.");
    mToleranceElasticCrackAngleHigh = rParameter;
}

//! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
void NuTo::StructureIp::SetToleranceElasticCrackAngleLow(double rParameter)
{
    if (rParameter<0 || rParameter>mToleranceElasticCrackAngleHigh)
        throw MechanicsException("[NuTo::StructureIp::SetToleranceElasticCrackAngleLow] Tolerance has to be positive and lower than the high tolerance.");
    mToleranceElasticCrackAngleLow = rParameter;
}

//! @brief add a constraint equation for the crack opening (normal crack opening non negativ)
//! @parameter rPenaltyStiffness penalty stiffness for augmented Lagrangian
//! @return id of the constraint
int NuTo::StructureIp::ConstraintLagrangeCrackOpening(double rPenaltyStiffness)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintNormalCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 1;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    ConstraintLagrangeGlobalCrackOpening2D *mConst = new NuTo::ConstraintLagrangeGlobalCrackOpening2D(this, rPenaltyStiffness);

    mConstraintMap.insert(id, mConst);
    mConstraintNormalCrackOpening = id;
    return id;
}

//! @brief add a constraint equation for the total strain
//! @parameter rStrain applied strain (rhs)
//! @return id of the constraint
int NuTo::StructureIp::ConstraintLinearGlobalTotalStrain(const EngineeringStrain2D& rStrain)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTotalStrain);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 1;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    NuTo::ConstraintLinearGlobalTotalStrain *mConst = new NuTo::ConstraintLinearGlobalTotalStrain(this, rStrain);

    mConstraintMap.insert(id, mConst);
    mConstraintTotalStrain = id;
    return id;
}

//! @brief initializes some variables etc. before the Newton-Raphson routine is executed
void NuTo::StructureIp::InitBeforeNewtonRaphson()
{
    mSavedToStringStream = false;
    NewtonRaphsonInfo(10);
}

//! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
//! @param load factor
void NuTo::StructureIp::SetLoadFactor(double rLoadFactor)
{
    //set the total strain and calculate from the existing crack opening the homogeneous strain
    EngineeringStrain2D curEngineeringStrain(mPrevEpsilonTot+mDeltaEpsilonTot*rLoadFactor);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTotalStrain);
    if (it!=mConstraintMap.end())
    {
        dynamic_cast<NuTo::ConstraintLinearGlobalTotalStrain*>(it->second)->SetRHS(curEngineeringStrain);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureIp::SetLoadFactor] Linear constraint for total strain not found.");
    }
}

//! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureIp::PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor)const
{
    std::cout << " Finescale convergence after " << rNumNewtonIterations << " Newton iterations, curLoadFactor " << rLoadFactor << ", deltaLoadFactor "<< rDeltaLoadFactor << std::endl<< std::endl;
    std::cout << " crack angle " << mCrackAngle*180./M_PI << " Crack opening " << mCrackOpening[0] << " " << mCrackOpening[1] <<  std::endl;
}

//! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureIp::PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor)const
{
    std::cout << " Finescale step, iteration " << rNewtonIteration <<
                 ", line search factor " << rLineSearchFactor <<
                 ", load factor " << rLoadFactor << std::endl;
    std::cout << " crack angle " << mCrackAngle*180./M_PI << " Crack opening " << mCrackOpening[0] << " " << mCrackOpening[1] <<  std::endl;
    std::stringstream ssLoadStep;
    ssLoadStep << rLoadStep;
    std::stringstream ssIteration;
    ssIteration << rNewtonIteration;
    this->ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FinescaleConcurrentMultiscale") + ssLoadStep.str()+"_" + ssIteration.str() + std::string(".vtk"));
}

//! @brief only for debugging, info at some stages of the Newton Raphson iteration
void NuTo::StructureIp::NewtonRaphsonInfo(int rVerboseLevel)const
{
    std::cout << "DOFs : alpha "   << this->GetCrackAngle()                << "(" << this->GetDofCrackAngle()<< ") "
              << " tangential " << this->GetGlobalCrackOpening2D()[0]   << "(" << this->GetDofGlobalCrackOpening2D()[0] << ") "
              << " normal " << this->GetGlobalCrackOpening2D()[1]   << "(" << this->GetDofGlobalCrackOpening2D()[1] << ") "
              << std::endl;

}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::StructureIp)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
