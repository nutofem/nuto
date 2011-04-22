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
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/groups/GroupEnum.h"

#include <ANN/ANN.h>

#define PRINTRESULT true
//#define MYDEBUG
#define tolCrackOpeningConstraint 1e-5
//! @brief constructor
//! @param mDimension number of nodes
NuTo::StructureMultiscale::StructureMultiscale ( int rDimension)  : Structure ( rDimension )
{
    if (rDimension!=2)
        throw MechanicsException("[NuTo::StructureMultiscale::StructureMultiscale] The concurrent multiscale model is only implemented for 2D.");
    mCrackAngle = M_PI*(0.7000);
    mInitCrackAngle = 0.5*M_PI;
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
    mSquareCoarseScaleModel = true;
    mlCoarseScale = 0.;
    mlFineScaleDamage = 0.;
    mlFineScaleHomogeneous = 0.;
    mCenterDamage[0] = 0.;
    mCenterDamage[1] = 0.;
    mCenterHomogeneous[0] = 0.;
    mCenterHomogeneous[1] = 0.;
    mFineScaleAreaDamage = 0.;
    mFineScaleAreaHomogeneous = 0.;
    mCrackTransitionRadius = 0;
    mIPName = std::string("fineScaleIp");
    mGroupBoundaryNodesDamage = -1;
    mGroupBoundaryNodesHomogeneous = -1;
    mGroupElementsDamage = -1;
    mGroupElementsHomogeneous = -1;
    mBoundaryNodesElementsAssigned = false;
    mThickness = 1.;

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
void NuTo::StructureMultiscale::Info()const
{
    Structure::Info();
    //add crack angle and crack orientation

}

#ifdef ENABLE_SERIALIZATION
template void NuTo::StructureMultiscale::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::StructureMultiscale::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::StructureMultiscale::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::StructureMultiscale::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::StructureMultiscale::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::StructureMultiscale::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::StructureMultiscale::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of StructureMultiscale" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Structure)
       & BOOST_SERIALIZATION_NVP(mCrackAngle)
       & BOOST_SERIALIZATION_NVP(mDOFCrackAngle)
       & BOOST_SERIALIZATION_NVP(mInitCrackAngle)
       & BOOST_SERIALIZATION_NVP(mCrackOpening)
       & BOOST_SERIALIZATION_NVP(mDOFCrackOpening)
       & BOOST_SERIALIZATION_NVP(mEpsilonTot)
       & BOOST_SERIALIZATION_NVP(mDOFGlobalTotalStrain)
       & BOOST_SERIALIZATION_NVP(mEpsilonHom)
       & BOOST_SERIALIZATION_NVP(mPrevEpsilonHom)
       & BOOST_SERIALIZATION_NVP(mPrevAverageStress)
       & BOOST_SERIALIZATION_NVP(mSquareCoarseScaleModel)
       & BOOST_SERIALIZATION_NVP(mlCoarseScale)
       & BOOST_SERIALIZATION_NVP(mlFineScaleDamage)
       & BOOST_SERIALIZATION_NVP(mlFineScaleHomogeneous)
       & BOOST_SERIALIZATION_NVP(mCrackTransitionRadius)
       & BOOST_SERIALIZATION_NVP(mCenterDamage)
       & BOOST_SERIALIZATION_NVP(mCenterHomogeneous)
       & BOOST_SERIALIZATION_NVP(mFineScaleAreaDamage)
       & BOOST_SERIALIZATION_NVP(mFineScaleAreaHomogeneous)
       & BOOST_SERIALIZATION_NVP(mThickness)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleX)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleY)
       & BOOST_SERIALIZATION_NVP(mConstraintNormalCrackOpening)
       & BOOST_SERIALIZATION_NVP(mConstraintTangentialCrackOpening)
       & BOOST_SERIALIZATION_NVP(mConstraintCrackAngle)
       & BOOST_SERIALIZATION_NVP(mConstraintTotalStrain)
       & BOOST_SERIALIZATION_NVP(mBoundaryNodesElementsAssigned)
       & BOOST_SERIALIZATION_NVP(mGroupBoundaryNodesDamage)
       & BOOST_SERIALIZATION_NVP(mGroupBoundaryNodesHomogeneous)
       & BOOST_SERIALIZATION_NVP(mGroupMultiscaleNodesDamage)
       & BOOST_SERIALIZATION_NVP(mGroupMultiscaleNodesHomogeneous)
       & BOOST_SERIALIZATION_NVP(mGroupElementsDamage)
       & BOOST_SERIALIZATION_NVP(mGroupElementsHomogeneous)
       & BOOST_SERIALIZATION_NVP(mIPName)
       & BOOST_SERIALIZATION_NVP(mResultDirectoryAfterLineSearch)
       & BOOST_SERIALIZATION_NVP(mResultFileAfterConvergence)
       & BOOST_SERIALIZATION_NVP(mToleranceElasticCrackAngleLow)
       & BOOST_SERIALIZATION_NVP(mToleranceElasticCrackAngleHigh);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of StructureMultiscale" << std::endl;
#endif
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureMultiscale::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException("[NuTo::StructureMultiscale::Save] Error opening file.");
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
            throw MechanicsException ( "[NuTo::StructureMultiscale::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureMultiscale::Save]File save exception in boost - " ) + std::string ( e.what() ) );
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
void NuTo::StructureMultiscale::Restore (const std::string &filename, std::string rType )
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
            throw MechanicsException("[NuTo::StructureMultiscale::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureMultiscale::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureMultiscale::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureMultiscale::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::StructureMultiscale::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureMultiscale::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
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
        std::cout<<"[NuTo::StructureMultiscale::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
#endif// ENABLE_SERIALIZATION


void NuTo::StructureMultiscale::SaveStructure(std::stringstream& rSaveStringStream)const
{
#ifdef ENABLE_SERIALIZATION
    boost::archive::binary_oarchive oba(rSaveStringStream);
    oba << (*this);
#else
    throw MechanicsException("[NuTo::StructureMultiscale::Save] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

void NuTo::StructureMultiscale::RestoreStructure(std::stringstream& rSaveStringStream)
{
#ifdef ENABLE_SERIALIZATION
        boost::archive::binary_iarchive iba(rSaveStringStream);
        iba >> (*this);
#else
    throw MechanicsException("[NuTo::StructureMultiscale::Restore] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

//! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
//! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
void NuTo::StructureMultiscale::TransformMultiscaleNodes()
{
	if (GroupGetNumMembers(mGroupElementsDamage)+GroupGetNumMembers(mGroupElementsHomogeneous)!=GetNumElements())
	{
		std::cout << "damage elements " << GroupGetNumMembers(mGroupElementsDamage) << " homogeneous elements " << GroupGetNumMembers(mGroupElementsHomogeneous) << " total nodes" << GetNumElements() << std::endl;
		throw MechanicsException("[NuTo::StructureMultiscale::TransformMultiscaleNodes] there is something wrong with your elements groups");
	}
	if (GroupGetNumMembers(mGroupElementsDamage)+GroupGetNumMembers(mGroupElementsHomogeneous)!=GetNumElements())
		throw MechanicsException("[NuTo::StructureMultiscale::TransformMultiscaleNodes] there is something wrong with your elements");
	TransformMultiscaleNodes(mGroupBoundaryNodesDamage, mGroupMultiscaleNodesDamage, mGroupElementsDamage, mCenterDamage, mlFineScaleDamage, mFineScaleAreaDamage, true);
	TransformMultiscaleNodes(mGroupBoundaryNodesHomogeneous, mGroupMultiscaleNodesHomogeneous, mGroupElementsHomogeneous, mCenterHomogeneous, mlFineScaleHomogeneous, mFineScaleAreaHomogeneous,false);
}

//! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
//! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
void NuTo::StructureMultiscale::TransformMultiscaleNodes(int rGroupBoundaryNodes, int rGroupNodes, int rGroupElements, boost::array<double,2>& rCenter, double& rLength, double& rArea, bool rCrackedDomain)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (mBoundaryNodesElementsAssigned==false)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Groups of boundary nodes and elements (damage,homogeneous) has not been assigned yet.");
    boost::ptr_map<int,GroupBase>::iterator itGroupBoundary = mGroupMap.find(rGroupBoundaryNodes);
    if (itGroupBoundary==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Group with the given identifier does not exist.");
    if (itGroupBoundary->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Group for boundary nodes is not a node group.");
    Group<NodeBase>* nodeGroupBoundary = itGroupBoundary->second->AsGroupNode();

    //first calculate bounding box and dimension of specimen (round, square, size etc.)
    double minX(0.), minY(0.), maxX(0.), maxY(0.);

    //determine from the first three nodes the center and radius in case the model is round
    double vec1[2], vec2[2], coord1[2], coord2[2], coord3[2], midPoint1[2], midPoint2[2];
    if (nodeGroupBoundary->GetNumMembers()<3)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Number of nodes in the fine scale model has to be at least 3.");
    Group<NodeBase>::iterator itBoundaryNode=nodeGroupBoundary->begin();
    itBoundaryNode->second->GetCoordinates2D(coord1);itBoundaryNode++;
    itBoundaryNode->second->GetCoordinates2D(coord2);itBoundaryNode++;
    itBoundaryNode->second->GetCoordinates2D(coord3);
    //check for being on a straight line
    bool onLine(true);
    while (onLine)
    {
        vec1[0] = coord2[0]-coord1[0];
        vec1[1] = coord2[1]-coord1[1];
        vec2[0] = coord3[0]-coord1[0];
        vec2[1] = coord3[1]-coord1[1];
        std::cout << (vec1[0]*vec2[0]+vec2[1]*vec1[1])/(sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1])*sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1])) << std::endl;
        if (fabs((vec1[0]*vec2[0]+vec2[1]*vec1[1])/(sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1])*sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1])))>1.-1e-3)
        {
        	itBoundaryNode++;
            if (itBoundaryNode==nodeGroupBoundary->end())
                throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] All nodes are on the same line.");
            itBoundaryNode->second->GetCoordinates2D(coord3);
        }
        else
            onLine = false;
    }

    vec1[0] = -(coord2[1]-coord1[1]);
    vec1[1] = coord2[0]-coord1[0];
    vec2[0] = -(coord3[1]-coord1[1]);
    vec2[1] = coord3[0]-coord1[0];
    midPoint1[0] = 0.5*(coord2[0]+coord1[0]);
    midPoint1[1] = 0.5*(coord2[1]+coord1[1]);
    midPoint2[0] = 0.5*(coord3[0]+coord1[0]);
    midPoint2[1] = 0.5*(coord3[1]+coord1[1]);

    //calculate the intersection
    double s = -(-(midPoint1[0] -midPoint2[0])* vec2[1]+(midPoint1[1] -midPoint2[1])*vec2[0])/
            (-vec1[0]*vec2[1]+vec1[1]*vec2[0]);
    double t = -(-(midPoint1[0] -midPoint2[0])* vec1[1]+(midPoint1[1] -midPoint2[1])*vec1[0])/
            (-vec1[0]*vec2[1]+vec1[1]*vec2[0]);

    rCenter[0] = midPoint1[0] + s*vec1[0];
    rCenter[1] = midPoint1[1] + s*vec1[1];

    if (fabs(rCenter[0]-(midPoint2[0] + t*vec2[0]))>1e-3)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] calculation of midpoint is wrong.");
    if (fabs(rCenter[1]-(midPoint2[1] + t*vec2[1]))>1e-3)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] calculation of midpoint is wrong.");
    double fineScaleRadius = sqrt((coord1[0]-rCenter[0])*(coord1[0]-rCenter[0])+(coord1[1]-rCenter[1])*(coord1[1]-rCenter[1]));
    //check for a round fine scale model
    bool squareFineScaleModel = false;
    for (Group<NodeBase>::iterator itNode=nodeGroupBoundary->begin(); itNode!=nodeGroupBoundary->end();itNode++)
    {
        double coordinates[2];
        itNode->second->GetCoordinates2D(coordinates);
        if (squareFineScaleModel==false)
        {
            if (fabs((rCenter[0]-coordinates[0])*(rCenter[0]-coordinates[0])+(rCenter[1]-coordinates[1])*(rCenter[1]-coordinates[1])-fineScaleRadius*fineScaleRadius)>0.001*fineScaleRadius*fineScaleRadius)
            {
                std::cout << (rCenter[0]-coordinates[0])*(rCenter[0]-coordinates[0])+(rCenter[1]-coordinates[1])*(rCenter[1]-coordinates[1]) << " " << fineScaleRadius*fineScaleRadius << std::endl;
                squareFineScaleModel = true;
            }
        }
    }

    // transform all the standard nodes into multiscale nodes
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupNodes);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] Group for boundary nodes is not a node group.");
    Group<NodeBase>* nodeGroup = itGroup->second->AsGroupNode();
    //copy the boundary nodes first, since the iterators are invalidated in NodeExchangePtr
    std::vector<std::pair<int, NodeBase*> >nodeVec(nodeGroup->size());
    unsigned int countNode=0;
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end(); itNode++, countNode++)
    {
        nodeVec[countNode].first  = itNode->first;
        nodeVec[countNode].second = itNode->second;
    }

    for (countNode=0; countNode<nodeVec.size(); countNode++)
    {
        Node::eNodeType nodeType = nodeVec[countNode].second->GetNodeType();
        NodeBase* newNode(0);

        double coordinates[2];
        nodeVec[countNode].second->GetCoordinates2D(coordinates);

        switch (nodeType)
        {
        case Node::NodeCoordinatesDisplacements2D:
            newNode = new NodeCoordinatesDisplacementsMultiscale2D(this, rCrackedDomain);
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
            //old node is deleted in Exchange routine when being removed from node table
            NodeExchangePtr(nodeVec[countNode].first,nodeVec[countNode].second,newNode);
            break;
        case Node::NodeCoordinatesDisplacementsMultiscale2D:
            break;
        break;
        default:
            throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] node type not implemented.");
        }

    }
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0) = 1.;
    direction(1,0) = 0.;
    mConstraintFineScaleX = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    direction(0,0) = 0.;
    direction(1,0) = 1.;
    mConstraintFineScaleY = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    ConstraintInfo(10);
    if (minX==maxX || minY==maxY)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] structure has zero width or height, either check your boundary group or the full structure.");
    if (squareFineScaleModel)
    {
        if (fabs((maxX-minX)-(maxY-minY))>1e-3)
            throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] domain is not square.");
        rCenter[0] = 0.5*(maxX+minX);
        rCenter[1] = 0.5*(maxY+minY);
        rArea = (maxX-minX)* (maxY-minY);
        if (mlCoarseScale<=0)
            throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] coarse length is not correct (<=0)");
        rLength = maxX-minX;
        std::cout<< "square fine scale model" << std::endl;
    }
    else
    {
        boost::ptr_map<int,ElementBase>::iterator elementIter;
        //this area calculation takes into account the missing area between the circle and the secants
        rArea = nodeGroupBoundary->GetNumMembers()*fineScaleRadius*fineScaleRadius*0.5*sin(2.*M_PI/nodeGroupBoundary->GetNumMembers());
        if (mlCoarseScale<=0)
            throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] coarse length is not correct (<=0)");
        rLength = 2.*fineScaleRadius;
        std::cout<< "round fine scale model" << std::endl;
    }
    std::cout<< "center of fine scale at (" << rCenter[0] << "," << rCenter[1] << "), area(" << rArea << "), fine scale length "<< rLength << std::endl;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureMultiscale::TransformBoundaryNodes] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief numbers non standard DOFs' e.g. in StructureMultiscale, for standard structures this routine is empty
void NuTo::StructureMultiscale::NumberAdditionalGlobalDofs()
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
        std::cout << "mDofs for crack angle "<< mDOFCrackAngle << std::endl;
    }
    else
        throw MechanicsException("[NuTo::StructureMultiscale::SetAdditionalGlobalDofs] Only implemented for 2D.");
}

// merge dof values
void NuTo::StructureMultiscale::NodeMergeAdditionalGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
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
        throw MechanicsException("[NuTo::StructureMultiscale::NodeMergeActiveDofValues] Only implemented for 2D");

    //update the homogeneous strain
    CalculateHomogeneousEngineeringStrain();
}

//! @brief extract dof values additional dof values
//! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
//! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
void NuTo::StructureMultiscale::NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const
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

void NuTo::StructureMultiscale::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
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
    std::vector<boost::array<double,6> > dDOF;
    std::vector<boost::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);
    double bAlphaRow, bURow[2], bMRow[3];

    // loop over all elements in the damaged region
	boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
    const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorDamage = this->GetScalingFactorDamage();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
    	itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
    	elementVector*=scalingFactorDamage;
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
                {
                    rActiveDofGradientVector(mDOFCrackOpening[1],0) += bURow[1] * elementVector(rowCount,0);
                }
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
    }
    // loop over all elements in the damaged region
	itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) is not an element group.");
    elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
    	itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
    	elementVector*=scalingFactorHomogeneous;
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
                bURow[0] = dDOF[theDofMapRow][1];;
                bURow[1] = dDOF[theDofMapRow][2];;

                //influence of epsilon_tot = epsilon_hom
                bMRow[0] = dDOF[theDofMapRow][3];;
                bMRow[1] = dDOF[theDofMapRow][4];;
                bMRow[2] = dDOF[theDofMapRow][5];;

                if (mDOFCrackAngle< this->mNumActiveDofs)
                    rActiveDofGradientVector(mDOFCrackAngle,0) += bAlphaRow * elementVector(rowCount,0);
                else
                    rDependentDofGradientVector(mDOFCrackAngle - this->mNumActiveDofs,0) += bAlphaRow * elementVector(rowCount,0);

                if (mDOFCrackOpening[0]< this->mNumActiveDofs)
                    rActiveDofGradientVector(mDOFCrackOpening[0],0) += bURow[0] * elementVector(rowCount,0);
                else
                    rDependentDofGradientVector(mDOFCrackOpening[0] - this->mNumActiveDofs,0) += bURow[0] * elementVector(rowCount,0);

                if (mDOFCrackOpening[1]< this->mNumActiveDofs)
                {
                    rActiveDofGradientVector(mDOFCrackOpening[1],0) += bURow[1] * elementVector(rowCount,0);
                }
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
    }

    //write contribution of Lagrange Multipliers
    ConstraintBuildGlobalGradientInternalPotentialSubVectors(rActiveDofGradientVector, rDependentDofGradientVector);
}


void NuTo::StructureMultiscale::CalculateCohesiveForce(NuTo::FullMatrix<double>& rCohesiveForce) const
{
    // initialize vectors
    assert(rCohesiveForce.GetNumRows() == 2);
    assert(rCohesiveForce.GetNumColumns() == 1);
    rCohesiveForce(0,0) = 0.0;
    rCohesiveForce(1,0) = 0.0;

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to ux, uy without assuming this to be balanced with epsilon_hom
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<boost::array<double,2> > dDOF;
    CalculatedDispdCrackOpening(mappingDofMultiscaleNode,dDOF);
    double bURow[2];

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

            //add contribution of the global degrees of freedom (crack opening)
            if (mappingDofMultiscaleNode[globalRowDof]!=-1)
            {
                //influence of alpha
                int theDofMapRow = mappingDofMultiscaleNode[globalRowDof];

                bURow[0] = dDOF[theDofMapRow][0];
                bURow[1] = dDOF[theDofMapRow][1];

                rCohesiveForce(0,0) += bURow[0] * elementVector(rowCount,0);
                rCohesiveForce(1,0) += bURow[1] * elementVector(rowCount,0);
            }
        }
        elementIter++;
    }
}



// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
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
    std::vector<boost::array<double,6> > dDOF;
    std::vector<boost::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);

    std::cout << "alpha "  << mCrackAngle << " crack opening "<< mCrackOpening[0] << " " << mCrackOpening[1] <<
    		     " epsilon hom" << mEpsilonHom.GetData()[0] << " "<< mEpsilonHom.GetData()[1]<< " "<< mEpsilonHom.GetData()[2]<<std::endl;

    // loop over all elements in the damaged region
	boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
    const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorDamage = this->GetScalingFactorDamage();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
        bool symmetryFlag = false;
        itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
        elementMatrix*=scalingFactorDamage;
        elementVector*=scalingFactorDamage;

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
    }

    // loop over all elements in the homogeneous region
	itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(homogeneous) is not an element group.");
    elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
        bool symmetryFlag = false;
        itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
        elementMatrix*=scalingFactorHomogeneous;
        elementVector*=scalingFactorHomogeneous;

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
/*    {
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
void NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK) const
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
    std::vector<boost::array<double,6> > dDOF;
    std::vector<boost::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);

    // loop over all elements in the damaged region
	boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
    const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorDamage = this->GetScalingFactorDamage();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
        bool symmetryFlag = false;
        itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
        elementMatrix*=scalingFactorDamage;
        elementVector*=scalingFactorDamage;

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
    }
    // loop over all elements in the homogeneous region
	itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) is not an element group.");
    elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        // calculate element contribution
        bool symmetryFlag = false;
        itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
        elementMatrix*=scalingFactorHomogeneous;
        elementVector*=scalingFactorHomogeneous;

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
    }

    /*
    {
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
/*
    {
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

void NuTo::StructureMultiscale::AddElementMatrixToGlobalSubMatricesGeneral(
        NuTo::FullMatrix<double>& elementMatrix,
        std::vector<int>& elementMatrixGlobalDofsRow,
        std::vector<int>& elementMatrixGlobalDofsColumn,
        NuTo::FullMatrix<double>& elementVector,
        std::vector<int>& mappingDofMultiscaleNode,
        std::vector<boost::array<double,6> >&dDOF,
        std::vector<boost::array<double,3> >&dDOF2,
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
    double bAlphaRow(0), bURow[2]={ 0, 0}, bMRow[3]={0,0,0}, bAlphaCol(0), bUCol[2]={0,0}, bMCol[3]={0,0,0};

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
                    //std::cout.precision(15);
                    //std::cout << "add F(" << elementVector(rowCount, 0) << ")*dDOF2(" << dDOF2[theDofMapRow][1] << ")=" << elementVector(rowCount, 0) * dDOF2[theDofMapRow][1] <<  std::endl;
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
                            //std::cout << "add at ("<< mDOFCrackAngle << "," << mDOFCrackOpening[0]<< ") " << bAlphaRow<< "(" << bAlphaRow << ")*elementMatrix(" << elementMatrix(rowCount, colCount) << ")*bUCol(" <<  bUCol[0] << ")=" << bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0] << std::endl;
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

//! @brief calculate the derivative of the displacements at the nodes with respect to crack opening without considering this to be balanced by epsilon_hom
//! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
//! @param rDOF return value, for each dof, the corresponding derivatives (ux, uy)
void NuTo::StructureMultiscale::CalculatedDispdCrackOpening(std::vector<int>& rMappingDofMultiscaleNode,std::vector<boost::array<double,2> >& rDOF)const
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
                throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] This multiscale node type has not been implemented.");
        }
    }

    rDOF.resize(numMultiscaleDofs);
    int countMultiscaleDofs(0);
    //second loop, calculate derivates
    for (boost::ptr_map<int,NodeBase>::const_iterator itNode=mNodeMap.begin(); itNode!=mNodeMap.end(); itNode++)
    {
        if (itNode->second->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
        {
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(0)] = countMultiscaleDofs;
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(1)] = countMultiscaleDofs+1;

            double coord[2];
            itNode->second->GetCoordinates2D(coord);

            //derivative of displacement with respect to discontinuity (crack orientation)
            GetdDisplacementdCrackOpening(coord,&(rDOF[countMultiscaleDofs][0]),&(rDOF[countMultiscaleDofs+1][0]));
            countMultiscaleDofs+=2;

            if (countMultiscaleDofs>numMultiscaleDofs)
            {
                throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdCrackOpening] countMultiscaleDofs>numMultiscaleDofs internal error.");
            }
        }
    }
    if (countMultiscaleDofs!=numMultiscaleDofs)
        throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdCrackOpening] countMultiscaleDofs!=numMultiscaleDofs internal error.");
}


//! @brief calculate the derivative of the displacements at the nodes with respect to crack opening and crack orientation
//! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
//! @param rDOF return value, for each dof, the corresponding derivatives alpha, ux, uy, exx, eyy, gamma_xy  [0..5]
//! @param r2DOF return value, for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy)
void NuTo::StructureMultiscale::CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<boost::array<double,6> >& rDOF, std::vector<boost::array<double,3> >& rDOF2)const
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
                throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] This multiscale node type has not been implemented.");
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
    const_cast<StructureMultiscale*>(&*this)->mCrackAngle += interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    EngineeringStrain2D strain2 = mEpsilonHom;
    std::cout << "cdf alpha " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
                        << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
                        << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< std::endl;
    const_cast<StructureMultiscale*>(&*this)->mCrackAngle -= interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();

    for (int count=0; count<2; count++)
    {
        std::cout << "algo u" << bHomU[0+3*count] << "  " << bHomU[1+3*count] << "  " <<  bHomU[2+3*count] << std::endl;
        strain1 = mEpsilonHom;
        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count]+=interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
        EngineeringStrain2D strain2 = mEpsilonHom;
        std::cout << "cdf " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
                            << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
                            << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< std::endl;
        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count] -= interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    }
    std::cout << "crack opening1 " << mCrackOpening[0] << " " << mCrackOpening[1]<< std::endl;
    //check hessian
    double hessianCDF[9],bHomAlpha2[3],bHomU2[6];
    std::cout << "algo hessian " << std::endl;
    std::cout << bHessian[0] << "  " << bHessian[3] << "  " <<  bHessian[6] << std::endl;
    std::cout << bHessian[1] << "  " << bHessian[4] << "  " <<  bHessian[7] << std::endl;
    std::cout << bHessian[2] << "  " << bHessian[5] << "  " <<  bHessian[8] << std::endl;
    const_cast<StructureMultiscale*>(&*this)->mCrackAngle += interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    GetdEpsilonHomdCrack(bHomAlpha2,bHomU2,0);
    hessianCDF[0] = (bHomAlpha2[0]-bHomAlpha[0])/interval;
    hessianCDF[1] = (bHomAlpha2[1]-bHomAlpha[1])/interval;
    hessianCDF[2] = (bHomAlpha2[2]-bHomAlpha[2])/interval;

    const_cast<StructureMultiscale*>(&*this)->mCrackAngle -= interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();

    for (int count=0; count<2; count++)
    {
        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count]+=interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
        GetdEpsilonHomdCrack(bHomAlpha2,bHomU2,0);
        hessianCDF[3+3*count] = (bHomAlpha2[0]-bHomAlpha[0])/interval;
        hessianCDF[4+3*count] = (bHomAlpha2[1]-bHomAlpha[1])/interval;
        hessianCDF[5+3*count] = (bHomAlpha2[2]-bHomAlpha[2])/interval;

        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count] -= interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    }
    std::cout << "cdf hessian " << std::endl;
    std::cout << hessianCDF[0] << "  " << hessianCDF[3] << "  " <<  hessianCDF[6] << std::endl;
    std::cout << hessianCDF[1] << "  " << hessianCDF[4] << "  " <<  hessianCDF[7] << std::endl;
    std::cout << hessianCDF[2] << "  " << hessianCDF[5] << "  " <<  hessianCDF[8] << std::endl;

    if (fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 ||
        fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 ||
        fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 )
        throw MechanicsException("CalculatedDispdGlobalDofs hessian is not correct.");

    std::cout << std::endl;
#endif

    //second loop, calculate derivates
    for (boost::ptr_map<int,NodeBase>::const_iterator itNode=mNodeMap.begin(); itNode!=mNodeMap.end(); itNode++)
    {
    	if (itNode->second->GetNodeType()==Node::NodeCoordinatesDisplacementsMultiscale2D)
        {
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(0)] = countMultiscaleDofs;
            rMappingDofMultiscaleNode[itNode->second->GetDofFineScaleDisplacement(1)] = countMultiscaleDofs+1;

            bool nodeInCrackedDomain = itNode->second->IsInCrackedDomain();

            double coord[2];
            itNode->second->GetCoordinates2D(coord);
            //derivative of displacement with respect to homogeneous strain
            double dDOFdEpsilonHomX[3];
            double dDOFdEpsilonHomY[3];
            if (nodeInCrackedDomain)
            {
                GetdDisplacementdEpsilonHom(coord,dDOFdEpsilonHomX,dDOFdEpsilonHomY,mCenterDamage);
            }
            else
            {
            	GetdDisplacementdEpsilonHom(coord,dDOFdEpsilonHomX,dDOFdEpsilonHomY,mCenterHomogeneous);
            }
#ifdef MYDEBUG
            {
            double interval(1e-3);
            double disp1[2],disp2[2];
            double dDOFdEpsilonHomXcdf[3];
            double dDOFdEpsilonHomYcdf[3];
            std::cout << "dDofxdEpsilonHom algo " << dDOFdEpsilonHomX[0] << "  " << dDOFdEpsilonHomX[1] << "  " << dDOFdEpsilonHomX[2] << std::endl;
            std::cout << "dDofydEpsilonHom algo " << dDOFdEpsilonHomY[0] << "  " << dDOFdEpsilonHomY[1] << "  " << dDOFdEpsilonHomY[2] << std::endl;

            if (nodeInCrackedDomain)
                GetDisplacementsEpsilonHom2D(coord, disp1, mCenterDamage);
            else
            	GetDisplacementsEpsilonHom2D(coord, disp1, mCenterHomogeneous);

            for (int count=0; count<3; count++)
            {
                const_cast<StructureMultiscale*>(&*this)->mEpsilonHom.mEngineeringStrain[count] += interval;
                if (nodeInCrackedDomain)
                    GetDisplacementsEpsilonHom2D(coord, disp2, mCenterDamage);
                else
                	GetDisplacementsEpsilonHom2D(coord, disp2, mCenterHomogeneous);
                dDOFdEpsilonHomXcdf[count] = (disp2[0]-disp1[0])/interval;
                dDOFdEpsilonHomYcdf[count] = (disp2[1]-disp1[1])/interval;
                const_cast<StructureMultiscale*>(&*this)->mEpsilonHom.mEngineeringStrain[count] -= interval;
            }
            std::cout << "dDofxdEpsilonHom cdf " << dDOFdEpsilonHomXcdf[0] << "  " << dDOFdEpsilonHomXcdf[1] << "  " << dDOFdEpsilonHomXcdf[2] << std::endl;
            std::cout << "dDofydEpsilonHom cdf " << dDOFdEpsilonHomYcdf[0] << "  " << dDOFdEpsilonHomYcdf[1] << "  " << dDOFdEpsilonHomYcdf[2] << std::endl;
            std::cout << std::endl;
            }
#endif

            //derivative of displacement with respect to discontinuity (crack orientation)
            if (nodeInCrackedDomain)
            {
                GetdDisplacementdCrackOrientation(coord,&(rDOF[countMultiscaleDofs][0]),&(rDOF[countMultiscaleDofs+1][0]));
            }
            else
            {
            	rDOF[countMultiscaleDofs][0] = 0.;
            	rDOF[countMultiscaleDofs+1][0] = 0.;
            }

            //std::cout << "node in cracked domain " << nodeInCrackedDomain << std::endl;
            //std::cout << " dDofEpsilonX" <<  dDOFdEpsilonHomX[0] << " "<<  dDOFdEpsilonHomX[1] << " " << dDOFdEpsilonHomX[2] << std::endl;
            //std::cout << " dDofEpsilonY" <<  dDOFdEpsilonHomY[0] << " "<<  dDOFdEpsilonHomY[1] << " " << dDOFdEpsilonHomY[2] << std::endl;
            //std::cout << " bHomAlpha" <<  bHomAlpha[0] << " "<<  bHomAlpha[1] << " " << bHomAlpha[2] << std::endl;
            rDOF[countMultiscaleDofs][0]  +=  dDOFdEpsilonHomX[0]*bHomAlpha[0] +
                    dDOFdEpsilonHomX[1]*bHomAlpha[1] +
                    dDOFdEpsilonHomX[2]*bHomAlpha[2];

            rDOF[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHomAlpha[0] +
                    dDOFdEpsilonHomY[1]*bHomAlpha[1] +
                    dDOFdEpsilonHomY[2]*bHomAlpha[2];
#ifdef MYDEBUG
            {
            double interval(1e-8);
            double disp1[2],disp2[2],disptmp[2];
            std::cout << "dDofdAlpha algo " << rDOF[countMultiscaleDofs][0] << "  " << rDOF[countMultiscaleDofs+1][0] << std::endl;
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp1);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp1[0] += disptmp[0];
                disp1[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp1,mCenterHomogeneous);
            }
            //std::cout << "coord " << coord[0] << " " << coord[1] << " disp " << disp1[0] << " " <<  disp1[1] << std::endl;
            const_cast<StructureMultiscale*>(&*this)->mCrackAngle += interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp2);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp2[0] += disptmp[0];
                disp2[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp2,mCenterHomogeneous);
            }
            //std::cout << "coord " << coord[0] << " " << coord[1] << " disp " << disp2[0] << " " <<  disp2[1] << std::endl;
            std::cout << "dDofdAlpha cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][0])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][0])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureMultiscale*>(&*this)->mCrackAngle -= interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            }
#endif
            //derivative of displacement with respect to discontinuity (crack opening)
            if (nodeInCrackedDomain)
            {
                GetdDisplacementdCrackOpening(coord,&(rDOF[countMultiscaleDofs][1]),&(rDOF[countMultiscaleDofs+1][1]));
            }
            else
            {
            	rDOF[countMultiscaleDofs][1] = 0.;
            	rDOF[countMultiscaleDofs+1][1] = 0.;
            	rDOF[countMultiscaleDofs][2] = 0.;
            	rDOF[countMultiscaleDofs+1][2] = 0.;
            }

            rDOF[countMultiscaleDofs][1]  += dDOFdEpsilonHomX[0]*bHomU[0] +
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
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp1);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp1[0] += disptmp[0];
                disp1[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp1,mCenterHomogeneous);
            }
            const_cast<StructureMultiscale*>(&*this)->mCrackOpening[0] += interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp2);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp2[0] += disptmp[0];
                disp2[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp2,mCenterHomogeneous);
            }
            std::cout << "dDofdU1 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][1])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][1])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureMultiscale*>(&*this)->mCrackOpening[0] -= interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            }
#endif
            rDOF[countMultiscaleDofs][2]  += dDOFdEpsilonHomX[0]*bHomU[3] +
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
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp1);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp1[0] += disptmp[0];
                disp1[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp1,mCenterHomogeneous);
            }
            const_cast<StructureMultiscale*>(&*this)->mCrackOpening[1] += interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            if (nodeInCrackedDomain)
            {
                GetDisplacementsCrack2D(coord, disp2);
                GetDisplacementsEpsilonHom2D(coord, disptmp,mCenterDamage);
                disp2[0] += disptmp[0];
                disp2[1] += disptmp[1];
            }
            else
            {
                GetDisplacementsEpsilonHom2D(coord, disp2,mCenterHomogeneous);
            }
            std::cout << "dDofdU2 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << std::endl;
            if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][2])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][2])>1e-2)
                throw MechanicsException("Hier ist ein Fehler.");
            const_cast<StructureMultiscale*>(&*this)->mCrackOpening[1] -= interval;
            const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
            std::cout << std::endl;
            std::cout << "crack opening2 " << mCrackOpening[0] << " " << mCrackOpening[1]<< std::endl;

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
            if (nodeInCrackedDomain)
            {
                Getd2Displacementd2CrackOrientation(coord,&(rDOF2[countMultiscaleDofs][0]),&(rDOF2[countMultiscaleDofs+1][0]));
            }
            else
            {
            	rDOF2[countMultiscaleDofs][0] = 0.;
            	rDOF2[countMultiscaleDofs+1][0] = 0.;
            }

            rDOF2[countMultiscaleDofs][0]   += dDOFdEpsilonHomX[0]*bHessian[0] +
                                               dDOFdEpsilonHomX[1]*bHessian[1] +
                                               dDOFdEpsilonHomX[2]*bHessian[2];
            rDOF2[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHessian[0] +
                                               dDOFdEpsilonHomY[1]*bHessian[1] +
                                               dDOFdEpsilonHomY[2]*bHessian[2];

            //derivative of displacement with respect to alpha and discontinuity (crack opening)
            if (nodeInCrackedDomain)
            {
                Getd2Displacementd2CrackOpening(coord, &(rDOF2[countMultiscaleDofs][1]), &(rDOF2[countMultiscaleDofs+1][1]));
            }
            else
            {
            	rDOF2[countMultiscaleDofs][1] = 0.;
            	rDOF2[countMultiscaleDofs][2] = 0.;
            	rDOF2[countMultiscaleDofs+1][1] = 0.;
            	rDOF2[countMultiscaleDofs+1][2] = 0.;
            }
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
                throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofs] countMultiscaleDofs>numMultiscaleDofs internal error.");
            }
        }
    }
#ifdef MYDEBUG
    std::cout << "End of routine" << std::endl << std::endl;
#endif
    if (countMultiscaleDofs!=numMultiscaleDofs)
        throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofs] countMultiscaleDofs!=numMultiscaleDofs internal error.");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
void NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const
{
    throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric] not implemented for StructureMultiscale");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
//! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
void NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const
{
    throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric] not implemented for StructureMultiscale");
}

//! @brief ... calculate the displacement based on the homogeneous strain
//! @param rCoordinates ... coordinates of the point
//! @param rCoarseDisplacements ... return array of displacements
void NuTo::StructureMultiscale::GetDisplacementsEpsilonHom2D(double rCoordinates[2], double rDisplacements[2], const boost::array<double,2>& rCenter)const
{
	rDisplacements[0] = mEpsilonHom.mEngineeringStrain[0] * (rCoordinates[0]-rCenter[0]) + 0.5 * mEpsilonHom.mEngineeringStrain[2] * (rCoordinates[1]-rCenter[1]);
    rDisplacements[1] = mEpsilonHom.mEngineeringStrain[1] * (rCoordinates[1]-rCenter[1]) + 0.5 * mEpsilonHom.mEngineeringStrain[2] * (rCoordinates[0]-rCenter[0]);
}

//! @brief ... calculate the displacement based on the crack opening
//! @param rCoordinates ... coordinates of the point
//! @param rCoarseDisplacements ... return array of displacements
void NuTo::StructureMultiscale::GetDisplacementsCrack2D(double rCoordinates[2], double rDisplacements[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    double factor;
    if (d<-mCrackTransitionRadius)
    {
        factor = -0.5;
    }
    else
    {
        if (d>mCrackTransitionRadius)
        {
            factor = 0.5;
        }
        else
        {
            //smooth transition from cracking to none cracking
            factor = 0.5*sin(0.5*M_PI*d/mCrackTransitionRadius);
        }
    }

    rDisplacements[0] = factor*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
    rDisplacements[1] = factor*(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1]);

}

//! @brief derivative of displacement with respect to homogeneous strain
//! @param rdX_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
//! @param rdY_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
void NuTo::StructureMultiscale::GetdDisplacementdEpsilonHom(double rCoordinates[2], double rdX_dEpsilonHom[3], double rdY_dEpsilonHom[3], const boost::array<double,2>& rCenter)const
{
    rdX_dEpsilonHom[0] = rCoordinates[0]-rCenter[0];
    rdX_dEpsilonHom[1] = 0.;
    rdX_dEpsilonHom[2] = 0.5 * (rCoordinates[1]-rCenter[1]);

    rdY_dEpsilonHom[0] = 0.;
    rdY_dEpsilonHom[1] = rCoordinates[1]-rCenter[1];
    rdY_dEpsilonHom[2] = 0.5 * (rCoordinates[0]-rCenter[0]);
}

//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
//! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to crack opening (ux, uy)
void NuTo::StructureMultiscale::GetdDisplacementdCrackOpening(double rCoordinates[2], double rdX_dCrackOpening[2], double rdY_dCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double factor;

    if (d<-mCrackTransitionRadius)
    {
        factor = -0.5;
    }
    else
    {
        if (d>mCrackTransitionRadius)
        {
            factor = 0.5;
        }
        else
        {
            //smooth transition from cracking to none cracking
            factor=0.5*sin(0.5*M_PI*d/mCrackTransitionRadius);
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
void NuTo::StructureMultiscale::Getd2Displacementd2CrackOpening(double rCoordinates[2], double rdX_dAlphaCrackOpening[2], double rdY_dAlphaCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double factor;
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionRadius)
    {
        factor=-0.5;
        rdX_dAlphaCrackOpening[0] = -factor*sinAlpha;
        rdX_dAlphaCrackOpening[1] = -factor*cosAlpha;
        rdY_dAlphaCrackOpening[0] =  factor*cosAlpha;
        rdY_dAlphaCrackOpening[1] = -factor*sinAlpha;
    }
    else
    {
        if (d>mCrackTransitionRadius)
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
            double dFactor(0.25*M_PI/mCrackTransitionRadius*cos(0.5*M_PI*d/mCrackTransitionRadius)*dDdAlpha);
            factor = 0.5*sin(0.5*M_PI*d/mCrackTransitionRadius);

            rdX_dAlphaCrackOpening[0] = -factor*sinAlpha+dFactor*cosAlpha;
            rdX_dAlphaCrackOpening[1] = -factor*cosAlpha-dFactor*sinAlpha;
            rdY_dAlphaCrackOpening[0] =  factor*cosAlpha+dFactor*sinAlpha;
            rdY_dAlphaCrackOpening[1] = -factor*sinAlpha+dFactor*cosAlpha;
        }
    }
#ifdef MYDEBUG
    {
        double delta(1e-8);
        double rdX_dAlphaCrackOpeningCDF[2],rdY_dAlphaCrackOpeningCDF[2];
        double rdX_dCrackOpening1[2],rdX_dCrackOpening2[2],rdY_dCrackOpening1[2],rdY_dCrackOpening2[2];
        GetdDisplacementdCrackOpening(rCoordinates, rdX_dCrackOpening1, rdY_dCrackOpening1);
        const_cast<StructureMultiscale*>(this)->mCrackAngle+=delta;
        GetdDisplacementdCrackOpening(rCoordinates, rdX_dCrackOpening2, rdY_dCrackOpening2);
        rdX_dAlphaCrackOpeningCDF[0] = (rdX_dCrackOpening2[0] - rdX_dCrackOpening1[0])/delta;
        rdX_dAlphaCrackOpeningCDF[1] = (rdX_dCrackOpening2[1] - rdX_dCrackOpening1[1])/delta;
        rdY_dAlphaCrackOpeningCDF[0] = (rdY_dCrackOpening2[0] - rdY_dCrackOpening1[0])/delta;
        rdY_dAlphaCrackOpeningCDF[1] = (rdY_dCrackOpening2[1] - rdY_dCrackOpening1[1])/delta;
        if (fabs(rdX_dAlphaCrackOpeningCDF[0]-rdX_dAlphaCrackOpening[0])>1e-3 || fabs(rdX_dAlphaCrackOpeningCDF[1]-rdX_dAlphaCrackOpening[1])>1e-3 ||
            fabs(rdY_dAlphaCrackOpeningCDF[0]-rdY_dAlphaCrackOpening[0])>1e-3 || fabs(rdY_dAlphaCrackOpeningCDF[1]-rdY_dAlphaCrackOpening[1])>1e-3)
            throw MechanicsException("[Getd2Displacementd2CrackOpening] something is wrong");
        const_cast<StructureMultiscale*>(this)->mCrackAngle-=delta;
    }
#endif
}


//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, derivative of y-displacement with respect to crack orientation (alpha)
void NuTo::StructureMultiscale::GetdDisplacementdCrackOrientation(double rCoordinates[2], double rdX_dAlpha[1], double rdY_dAlpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionRadius)
    {
        rdX_dAlpha[0] = -0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
        rdY_dAlpha[0] = -0.5*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
    }
    else
    {
        if (d>mCrackTransitionRadius)
        {
            rdX_dAlpha[0] = 0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
            rdY_dAlpha[0] = 0.5*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double factor(0.5*sin(0.5*M_PI*d/mCrackTransitionRadius));
            double dFactor(0.25*M_PI/mCrackTransitionRadius*cos(0.5*M_PI*d/mCrackTransitionRadius)*dDdAlpha);
#ifdef MYDEBUG
            {
                double delta(1e-8);
                const_cast<StructureMultiscale*>(this)->mCrackAngle+=delta;
                double d2 = CalculateDistanceToCrack2D(rCoordinates);
                double factor2(0.5*sin(0.5*M_PI*d2/mCrackTransitionRadius));
                double dFactorCDF = (factor2 - factor)/delta;
                double dDdAlphaCDF = (d2 - d)/delta;
                std::cout << "dFactor " << dFactor<< " " << dFactorCDF << std::endl;
                std::cout << "dDdAlpha " << dDdAlpha<< " " << dDdAlphaCDF << std::endl;
                if (fabs(dFactor-dFactorCDF)>1e-3)
                    throw MechanicsException("[GetdDisplacementdCrackOrientation] dFactor calculation is wrong");
                const_cast<StructureMultiscale*>(this)->mCrackAngle-=delta;
            }
#endif

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
    const_cast<NuTo::StructureMultiscale*> (this)->mCrackAngle+=interval;
    GetDisplacementsCrack2D(rCoordinates, displacements2);
    const_cast<NuTo::StructureMultiscale*> (this)->mCrackAngle-=interval;
    std::cout << "rdX_dAlpha "<< rdX_dAlpha[0] << " " << (displacements2[0]-displacements1[0])/interval << std::endl;
    std::cout << "rdY_dAlpha "<< rdY_dAlpha[0] << " " << (displacements2[1]-displacements1[1])/interval << std::endl;
    std::cout << d << std::endl;
    if (fabs(rdX_dAlpha[0]-(displacements2[0]-displacements1[0])/interval)>1e-3)
        throw MechanicsException("[NuTo::StructureMultiscale::GetdDisplacementdCrackOrientation] Here is something wrong.");
    if (fabs(rdY_dAlpha[0]-(displacements2[1]-displacements1[1])/interval)>1e-3)
        throw MechanicsException("[NuTo::StructureMultiscale::GetdDisplacementdCrackOrientation] Here is something wrong.");
#endif
}
//! @brief second derivative of displacement with respect to orientation of discontinuity (crack angle)
//! @param rdX_dAlpha[2] return value, second derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, second derivative of y-displacement with respect to crack orientation (alpha)
void NuTo::StructureMultiscale::Getd2Displacementd2CrackOrientation(double rCoordinates[2], double rd2X_d2Alpha[1], double rd2Y_d2Alpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (d<-mCrackTransitionRadius)
    {
        rd2X_d2Alpha[0] = -0.5*(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
        rd2Y_d2Alpha[0] = -0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
    }
    else
    {
        if (d>mCrackTransitionRadius)
        {
            rd2X_d2Alpha[0] = 0.5*(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
            rd2Y_d2Alpha[0] = 0.5*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double d2DdAlpha = Calculated2DistanceToCrack2Dd2Alpha(rCoordinates);
            double factor(0.5*sin(0.5*M_PI*d/mCrackTransitionRadius));
            double dFactor(0.25*M_PI/mCrackTransitionRadius*cos(0.5*M_PI*d/mCrackTransitionRadius)*dDdAlpha);
            double dFactor2(0.25*M_PI/mCrackTransitionRadius*(-0.5*M_PI/mCrackTransitionRadius*sin(0.5*M_PI*d/mCrackTransitionRadius)*dDdAlpha*dDdAlpha
                                      +cos(0.5*M_PI*d/mCrackTransitionRadius)*d2DdAlpha));

            rd2X_d2Alpha[0] = (dFactor2-factor)*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1])+2.*dFactor*(-sinAlpha*mCrackOpening[0]-cosAlpha*mCrackOpening[1]);
            rd2Y_d2Alpha[0] = (dFactor2-factor)*(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1])+2.*dFactor*(cosAlpha*mCrackOpening[0]-sinAlpha*mCrackOpening[1]);
        }
    }
#ifdef MYDEBUG

    //check dd2d2alpha
    double interval(1e-8);
    double dX_dAlpha1,dY_dAlpha1,dX_dAlpha2,dY_dAlpha2;
    GetdDisplacementdCrackOrientation(rCoordinates, &dX_dAlpha1, &dY_dAlpha1);
    const_cast<NuTo::StructureMultiscale*> (this)->mCrackAngle+=interval;
    GetdDisplacementdCrackOrientation(rCoordinates, &dX_dAlpha2, &dY_dAlpha2);
    const_cast<NuTo::StructureMultiscale*> (this)->mCrackAngle-=interval;

    std::cout << "dd2d2alpha algo " << rd2X_d2Alpha[0] << " " << rd2Y_d2Alpha[0] << std::endl;
    std::cout << "dd2d2alpha cdf " << (dX_dAlpha2-dX_dAlpha1)/interval << " " << (dY_dAlpha2-dY_dAlpha1)/interval << std::endl;

    if (fabs(rd2X_d2Alpha[0]-(dX_dAlpha2-dX_dAlpha1)/interval)>1e-2)
        throw MechanicsException("[Getd2Displacementd2CrackOrientation] there is something wrong.");
    if (fabs(rd2Y_d2Alpha[0]-(dY_dAlpha2-dY_dAlpha1)/interval)>1e-2)
        throw MechanicsException("[Getd2Displacementd2CrackOrientation] there is something wrong.");
#endif
}

//! @brief calculate the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return distance to crack
double NuTo::StructureMultiscale::CalculateDistanceToCrack2D(double rCoordinates[2])const
{
    return sin(mCrackAngle)*(mCenterDamage[0]-rCoordinates[0]) - cos(mCrackAngle)*(mCenterDamage[1]-rCoordinates[1]);
}

//! @brief calculate the derivative of the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return derivative of distance to crack
double NuTo::StructureMultiscale::CalculatedDistanceToCrack2DdAlpha(double rCoordinates[2])const
{
    return cos(mCrackAngle)*(mCenterDamage[0]-rCoordinates[0]) + sin(mCrackAngle)*(mCenterDamage[1]-rCoordinates[1]);
}

//! @brief calculate the second derivative of the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return second derivative of distance to crack
double NuTo::StructureMultiscale::Calculated2DistanceToCrack2Dd2Alpha(double rCoordinates[2])const
{
    return -sin(mCrackAngle)*(mCenterDamage[0]-rCoordinates[0]) + cos(mCrackAngle)*(mCenterDamage[1]-rCoordinates[1]);
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureMultiscale::CalculateHomogeneousEngineeringStrain()
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double gamma(0.);
    if (mSquareCoarseScaleModel)
    {
        gamma = 1./std::max(fabs(sinAlpha),fabs(cosAlpha));
    }
    else
    {
        gamma = 1.;
    }
    double factor=gamma/(mlCoarseScale);

    mEpsilonHom.mEngineeringStrain[0] = mEpsilonTot.mEngineeringStrain[0] - factor * sinAlpha *(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[1] = mEpsilonTot.mEngineeringStrain[1] - factor * cosAlpha *(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[2] = mEpsilonTot.mEngineeringStrain[2] - factor * ((2.*cosAlpha*cosAlpha-1.)* mCrackOpening[0] - 2.*sinAlpha*cosAlpha * mCrackOpening[1]);
    /*    mEpsilonHom.mEngineeringStrain[0] = 0;
    mEpsilonHom.mEngineeringStrain[1] = 0;
    mEpsilonHom.mEngineeringStrain[2] = 0;

    std::cout << "calculation of e hom  " << std::endl;

//    std::cout << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << std::endl;

//#ifdef MYDEBUG
    std::cout << "mEpsilonTot " << mEpsilonTot.mEngineeringStrain[0] << " " << mEpsilonTot.mEngineeringStrain[1] << " " << mEpsilonTot.mEngineeringStrain[2] << std::endl;
    std::cout << "mCrackOpening " << mCrackOpening[0] << " " << mCrackOpening[1] << " " << factor << std::endl;
    std::cout << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << std::endl;
//#endif
*/
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureMultiscale::SetTotalEngineeringStrain(EngineeringStrain2D& rTotalEngineeringStrain)
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
void NuTo::StructureMultiscale::SetDeltaTotalEngineeringStrain(EngineeringStrain2D& rDeltaTotalEngineeringStrain)
{
    mDeltaEpsilonTot = rDeltaTotalEngineeringStrain;
}

//! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
void NuTo::StructureMultiscale::SetPrevTotalEngineeringStrain(EngineeringStrain2D& rPrevTotalEngineeringStrain)
{
    mPrevEpsilonTot = rPrevTotalEngineeringStrain;
}

//! @brief returns the total strain
const NuTo::EngineeringStrain2D& NuTo::StructureMultiscale::GetTotalEngineeringStrain()const
{
    return mEpsilonTot;
}

//! @brief returns the total strain
const NuTo::EngineeringStrain2D& NuTo::StructureMultiscale::GetHomogeneousEngineeringStrain()const
{
    return mEpsilonHom;
}

//! @brief set the homgeneous strain (just for test purpose)
void NuTo::StructureMultiscale::SetHomogeneousEngineeringStrain(NuTo::EngineeringStrain2D rStrain)
{
    mEpsilonHom = rStrain;
}

//! @brief return the previous crack angle
double NuTo::StructureMultiscale::GetInitCrackAngle()const
{
    return mInitCrackAngle;
}

//! @brief sets the previous crack angle
void NuTo::StructureMultiscale::SetInitCrackAngle(double rInitCrackAngle)
{
    mInitCrackAngle = rInitCrackAngle;
}

//! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack orientation and crack opening
//! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
//! @parameter rbHomAlpha dHom wrt alpha
//! @paramter rbHomU [0-2] wrt ux [3-5] wrt uy
//! @parameter bHessian depsilondalpha2[0-2], depsilondalphadux[3-5], depsilondalphaduy[6-8]
void NuTo::StructureMultiscale::GetdEpsilonHomdCrack(double rbHomAlpha[3], double rbHomU[6], double rbHessian[9])const
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    if (mSquareCoarseScaleModel)
    {
        if (fabs(sinAlpha)>fabs(cosAlpha))
        {
            double signdivl = sinAlpha<0 ? 1./mlCoarseScale : -1./mlCoarseScale;
            rbHomAlpha[0] = signdivl*(sinAlpha*mCrackOpening[0] + cosAlpha*mCrackOpening[1]);
            rbHomAlpha[1] = signdivl*(-sinAlpha*mCrackOpening[0] +(cosAlpha/(sinAlpha*sinAlpha)*(cosAlpha*cosAlpha-2.))*mCrackOpening[1]);
            rbHomAlpha[2] = signdivl*((cosAlpha/(sinAlpha*sinAlpha)*(2.*cosAlpha*cosAlpha-3.))*mCrackOpening[0] +2.* sinAlpha*mCrackOpening[1]);

            rbHomU[0] = signdivl*(-cosAlpha);
            rbHomU[1] = signdivl*cosAlpha;
            rbHomU[2] = signdivl*(2.*cosAlpha*cosAlpha-1.)/sinAlpha;

            rbHomU[3] = signdivl*(sinAlpha);
            rbHomU[4] = signdivl*(cosAlpha*cosAlpha)/sinAlpha;
            rbHomU[5] = signdivl*(-2.*cosAlpha);

            if (rbHessian!=0)
            {
                rbHessian[0] = signdivl*(cosAlpha*mCrackOpening[0] - sinAlpha*mCrackOpening[1]);;
                rbHessian[1] = signdivl*(-cosAlpha*mCrackOpening[0] + (cosAlpha*cosAlpha*(cosAlpha*cosAlpha-1.)+2.)/(sinAlpha*sinAlpha*sinAlpha)*mCrackOpening[1]);
                rbHessian[2] = signdivl*((cosAlpha*cosAlpha*(2.*cosAlpha*cosAlpha-3.)+3.)/(sinAlpha*sinAlpha*sinAlpha)*mCrackOpening[0] +2.* cosAlpha*mCrackOpening[1]);

                rbHessian[3] = signdivl*sinAlpha;
                rbHessian[4] = -signdivl*sinAlpha;
                rbHessian[5] = signdivl*(cosAlpha/(sinAlpha*sinAlpha)*(2.*cosAlpha*cosAlpha-3.));

                rbHessian[6] = signdivl*cosAlpha;
                rbHessian[7] = signdivl*(cosAlpha*(cosAlpha*cosAlpha-2.)/(sinAlpha*sinAlpha));
                rbHessian[8] = signdivl*2.* sinAlpha;
            }
        }
        else
        {
            double signdivl = cosAlpha<0 ? 1./mlCoarseScale : -1./mlCoarseScale;
            rbHomAlpha[0] = signdivl*(-cosAlpha*mCrackOpening[0] + ((cosAlpha*cosAlpha+1.)*sinAlpha/(cosAlpha*cosAlpha))*mCrackOpening[1]);
            rbHomAlpha[1] = signdivl*( cosAlpha*mCrackOpening[0] - sinAlpha*mCrackOpening[1]);
            rbHomAlpha[2] = signdivl*(-(sinAlpha/(cosAlpha*cosAlpha)*(2.*cosAlpha*cosAlpha+1.))*mCrackOpening[0] -2.* cosAlpha*mCrackOpening[1]);

            rbHomU[0] = signdivl*(-sinAlpha);
            rbHomU[1] = signdivl*sinAlpha;
            rbHomU[2] = signdivl*(2.*cosAlpha*cosAlpha-1.)/cosAlpha;

            rbHomU[3] = signdivl*(sinAlpha*sinAlpha)/cosAlpha;
            rbHomU[4] = signdivl*(cosAlpha);
            rbHomU[5] = signdivl*(-2.*sinAlpha);

            if (rbHessian!=0)
            {
                rbHessian[0] = signdivl*(sinAlpha*mCrackOpening[0] +(cosAlpha*cosAlpha*(cosAlpha*cosAlpha-1.)+2.)/(cosAlpha*cosAlpha*cosAlpha)*mCrackOpening[1]);
                rbHessian[1] = signdivl*(-sinAlpha*mCrackOpening[0] - cosAlpha*mCrackOpening[1]);
                rbHessian[2] = signdivl*(-(cosAlpha*cosAlpha*(2.*cosAlpha*cosAlpha-1.)+2)/(cosAlpha*cosAlpha*cosAlpha)*mCrackOpening[0] +2.* sinAlpha*mCrackOpening[1]);

                rbHessian[3] = -signdivl*cosAlpha;
                rbHessian[4] = signdivl*cosAlpha;
                rbHessian[5] = signdivl*(-sinAlpha/(cosAlpha*cosAlpha)*(2.*cosAlpha*cosAlpha+1.));

                rbHessian[6] = signdivl*(sinAlpha*(cosAlpha*cosAlpha+1.)/(cosAlpha*cosAlpha));
                rbHessian[7] = -signdivl*sinAlpha;
                rbHessian[8] = -signdivl*(2.*cosAlpha);
            }
        }
    }
    else
    {
        // coarse scale model is not a square or not aligned with the principal axes
        double factor = -1./(mlCoarseScale);
        rbHomAlpha[0] = factor*((1.-2.*cosAlpha*cosAlpha)*mCrackOpening[0] + 2.*sinAlpha*cosAlpha*mCrackOpening[1]);
        rbHomAlpha[1] = -rbHomAlpha[0];
        rbHomAlpha[2] = factor*((-4.*sinAlpha*cosAlpha)*mCrackOpening[0] + (2.-4.*cosAlpha*cosAlpha)*mCrackOpening[1]);

        rbHomU[0] = factor*(-sinAlpha*cosAlpha);
        rbHomU[1] = -rbHomU[0];
        rbHomU[2] = factor*(2.*cosAlpha*cosAlpha-1.);

        rbHomU[3] = factor*(sinAlpha*sinAlpha);
        rbHomU[4] = factor*(cosAlpha*cosAlpha);
        rbHomU[5] = factor*(-2.*sinAlpha*cosAlpha);

        if (rbHessian!=0)
        {
            rbHessian[0] = factor*(4.*sinAlpha*cosAlpha*mCrackOpening[0] + (4.*cosAlpha*cosAlpha-2.)*mCrackOpening[1]);
            rbHessian[1] = -rbHessian[0];
            rbHessian[2] = factor*((4.-8.*cosAlpha*cosAlpha)*mCrackOpening[0] + 8.*sinAlpha*cosAlpha*mCrackOpening[1]);

            rbHessian[3] = factor*(1.-2.*cosAlpha*cosAlpha);
            rbHessian[4] = -rbHessian[3];
            rbHessian[5] = factor*(-4.*sinAlpha*cosAlpha);

            rbHessian[6] = factor*2.*sinAlpha*cosAlpha;
            rbHessian[7] = -rbHessian[6];
            rbHessian[8] = factor*(2.-4.*cosAlpha*cosAlpha);
        }
    }
/*
    std::cout << "crack opening "  << mCrackOpening[0] << " " << mCrackOpening[1]  << std::endl;
    std::cout << "dhom dalpha " << rbHomAlpha[0] << " " <<  rbHomAlpha[1] << " " <<  rbHomAlpha[2] << std::endl;
    std::cout << "dhom dut " << rbHomU[0] << " " <<  rbHomU[1] << " " <<  rbHomU[2] << std::endl;
    std::cout << "dhom dun " << rbHomU[3] << " " <<  rbHomU[4] << " " <<  rbHomU[5] << std::endl;
*/
}

//! @brief renumbers the global dofs in the structure after
void NuTo::StructureMultiscale::ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
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
double  NuTo::StructureMultiscale::CalculateInitialCrackAngleElastic()const
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

//! @brief calculates the difference between the crack angle of the initial elastic solution and the current angle
//! attention, the periodicity of the crack angle has to be taken into account
double NuTo::StructureMultiscale::CalculateDeltaCrackAngleElastic()const
{
     double alpha_2 = mInitCrackAngle;

     double delta_alpha = mCrackAngle-alpha_2;
     while (fabs(delta_alpha>M_PI))
     {
         if (delta_alpha>0)
         {
             delta_alpha-=2.*M_PI;
         }
         else
         {
             delta_alpha+=2.*M_PI;
         }
     }
     //now delta_alpha is in [0..PI]
     //std::cout << "alpha " << mCrackAngle << " alpha2 " << alpha_2 << " delta " << delta_alpha << std::endl;
     return delta_alpha;
}

//! @brief add a constraint equation for alpha, which corresponds to an artificial spring
//! @parameter rPenaltyStiffness penalty stiffness
//! @parameter rScalingFactor scaling factor
//! @return id of the constraint
int NuTo::StructureMultiscale::ConstraintNonlinearCrackAngle(double rPenaltyStiffness, double rScalingFactor)
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
void NuTo::StructureMultiscale::SetPenaltyStiffnessCrackAngle(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintCrackAngle);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackAngle2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackAngle2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessCrackAngle] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetPenaltyStiffness(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessCrackAngle] There is no constraint for the crack angle.");
    }
}

//! @brief set the scaling factor for the nonlinear crack angle constraint
void NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorCrackAngle(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintCrackAngle);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackAngle2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackAngle2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorCrackAngle] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetScalingFactor(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorCrackAngle] There is no constraint for the crack angle.");
    }

}

//! @brief add a constraint equation for the tangential crack opening, which corresponds to an artificial spring
//! @parameter rPenaltyStiffness penalty stiffness
//! @parameter rScalingFactor scaling factor
//! @return id of the constraint
int NuTo::StructureMultiscale::ConstraintNonlinearTangentialCrackOpening(double rScalingFactor, double rPenaltyStiffness)
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
void NuTo::StructureMultiscale::SetPenaltyStiffnessTangentialCrackOpening(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackOpeningTangential2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackOpeningTangential2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessTangentialCrackOpening] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetPenaltyStiffness(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessTangentialCrackOpening] There is no constraint for the crack angle.");
    }
}

//! @brief set the scaling factor for the nonlinear TangentialCrackOpening
void NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorTangentialCrackOpening(double rParameter)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        ConstraintNonlinearGlobalCrackOpeningTangential2D *constraintPtr = dynamic_cast<ConstraintNonlinearGlobalCrackOpeningTangential2D*>(it->second);
        if (constraintPtr==0)
            throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorTangentialCrackOpening] Constraint for the crack angle is not of the type ConstraintNonlinearGlobalCrackAngle2D.");
        constraintPtr->SetScalingFactor(rParameter);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureMultiscale::SetPenaltyStiffnessScalingFactorTangentialCrackOpening] There is no constraint for the crack angle.");
    }

}

//! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
void NuTo::StructureMultiscale::SetToleranceElasticCrackAngleHigh(double rParameter)
{
    if (rParameter<mToleranceElasticCrackAngleLow || rParameter<0)
        throw MechanicsException("[NuTo::StructureMultiscale::SetToleranceElasticCrackAngleHigh] Tolerance has to be positive and higher than the low tolerance.");
    mToleranceElasticCrackAngleHigh = rParameter;
}

//! @brief set the tolerance for the transition between crack angle from principal strain and previous strain
void NuTo::StructureMultiscale::SetToleranceElasticCrackAngleLow(double rParameter)
{
    if (rParameter<0 || rParameter>mToleranceElasticCrackAngleHigh)
        throw MechanicsException("[NuTo::StructureMultiscale::SetToleranceElasticCrackAngleLow] Tolerance has to be positive and lower than the high tolerance.");
    mToleranceElasticCrackAngleLow = rParameter;
}

//! @brief add a constraint equation for the crack opening (normal crack opening non negativ)
//! @parameter rPenaltyStiffness penalty stiffness for augmented Lagrangian
//! @return id of the constraint
int NuTo::StructureMultiscale::ConstraintLagrangeCrackOpening(double rPenaltyStiffness)
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
int NuTo::StructureMultiscale::ConstraintLinearGlobalTotalStrain(const EngineeringStrain2D& rStrain)
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
void NuTo::StructureMultiscale::InitBeforeNewtonRaphson()
{
    NewtonRaphsonInfo(10);
}

//! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
//! @param load factor
void NuTo::StructureMultiscale::SetLoadFactor(double rLoadFactor)
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
        throw MechanicsException("[NuTo::StructureMultiscale::SetLoadFactor] Linear constraint for total strain not found.");
    }
    //std::cout << "set total strain " << curEngineeringStrain.mEngineeringStrain[0]<< " " << curEngineeringStrain.mEngineeringStrain[1]<<" " << curEngineeringStrain.mEngineeringStrain[2] << std::endl;
}

//! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureMultiscale::PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor)const
{
    std::cout << " Finescale convergence after " << rNumNewtonIterations << " Newton iterations, curLoadFactor " << rLoadFactor << ", deltaLoadFactor "<< rDeltaLoadFactor << std::endl;
    std::cout << " crack angle " << mCrackAngle*180./M_PI << "[0..360], Crack opening " << mCrackOpening[0] << "[t] " << mCrackOpening[1] << "[n]"<< std::endl;
    std::cout << " total strain " << mEpsilonTot.mEngineeringStrain[0] << " " << mEpsilonTot.mEngineeringStrain[1] << " " << mEpsilonTot.mEngineeringStrain[2]  << std::endl;
    std::cout << " hom strain " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2]  << std::endl;
    NuTo::FullMatrix<double> multiplier;
    if (mConstraintNormalCrackOpening!=-1)
    {
        ConstraintLagrangeGetMultiplier(mConstraintNormalCrackOpening,multiplier);
        std::cout << " lambda crack opening " << multiplier(0,0)  << std::endl;    //! @brief calculates the average stress

    }
/*    NuTo::FullMatrix<double> cohesiveForce(2,1);
    CalculateCohesiveForce(cohesiveForce);
    std::cout << " cohesive force " << cohesiveForce(0,0) << "(t) " << cohesiveForce(1,0) << "(n) " << std::endl;

    NuTo::FullMatrix<double> engineeringStress(3,1);
    ElementTotalGetAverageStress(mFineScaleAreaDamage, engineeringStress);
    std::cout << " average stress vector " << engineeringStress(0,0) << " " << engineeringStress(1,0)<< " " <<  engineeringStress(2,0) << std::endl;
    std::cout << " stress vector on the crack " << engineeringStress(0,0)*(-sin(mCrackAngle))+ engineeringStress(2,0)*cos(mCrackAngle)<< "(t) " ;
                                      std::cout << engineeringStress(2,0)*(-sin(mCrackAngle))+ engineeringStress(1,0)*cos(mCrackAngle)<< "(n) " << std::endl;
*/
    std::stringstream ssLoadStep;
    ssLoadStep << rLoadStep;
    if (rLoadFactor==1)
        this->ExportVtkDataFile(mResultDirectoryAfterLineSearch+std::string("/../FinescaleConcurrentMultiscale") + ssLoadStep.str()+ std::string(".vtk"));
}

//! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureMultiscale::PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor)const
{
    std::cout << " Finescale step, iteration " << rNewtonIteration <<
                 ", line search factor " << rLineSearchFactor <<
                 ", load factor " << rLoadFactor << std::endl;
    std::cout << " crack angle " << mCrackAngle*180./M_PI << "[0..360], Crack opening " << mCrackOpening[0] << "[t] " << mCrackOpening[1] << "[n]"<< std::endl;
    std::stringstream ssLoadStep;
    ssLoadStep << rLoadStep;
    std::stringstream ssIteration;
    ssIteration << rNewtonIteration;
    this->ExportVtkDataFile(mResultFileAfterConvergence);
}

//! @brief only for debugging, info at some stages of the Newton Raphson iteration
void NuTo::StructureMultiscale::NewtonRaphsonInfo(int rVerboseLevel)const
{
    std::cout << "DOFs : alpha "   << this->GetCrackAngle()                << "(" << this->GetDofCrackAngle()<< ") "
              << " tangential " << this->GetGlobalCrackOpening2D()[0]   << "(" << this->GetDofGlobalCrackOpening2D()[0] << ") "
              << " normal " << this->GetGlobalCrackOpening2D()[1]   << "(" << this->GetDofGlobalCrackOpening2D()[1] << ") "
              << std::endl;

}

//! @brief calculates the total energy of the system
//! @return total energy
double NuTo::StructureMultiscale::ElementGroupGetTotalEnergy(int rGroupId)const
{
	throw MechanicsException("[NuTo::StructureMultiscale::ElementGroupGetTotalEnergy] not implemented if you do, consider the different scaling of each element.");
	// e.g. make an intersection with the damage and homogeneous part into new groups, call the Structure routine for these tmp groups scale the energies.
    return 0;
}


//! @brief calculates the total energy of the system
//! @return total energy
double NuTo::StructureMultiscale::ElementTotalGetTotalEnergy()const
{
	double energyDamage(Structure::ElementGroupGetTotalEnergy(mGroupElementsDamage));
	double energyHomogeneous(Structure::ElementGroupGetTotalEnergy(mGroupElementsHomogeneous));
	std::cout << "energy damage " << energyDamage << " * " << GetScalingFactorDamage() << " = " << energyDamage*GetScalingFactorDamage() << std::endl;
	std::cout << "energy homogeneous " << energyHomogeneous << " * " << GetScalingFactorHomogeneous() << " = " << energyHomogeneous*GetScalingFactorHomogeneous() << std::endl;
	return energyDamage*GetScalingFactorDamage()+energyHomogeneous*GetScalingFactorHomogeneous();
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::StructureMultiscale)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
