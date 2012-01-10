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

# ifdef _OPENMP
#include <omp.h>
# endif

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearFineScaleDisplacementsPeriodic2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalCrackOpening.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalTotalStrain.h"
#include "nuto/mechanics/constraints/ConstraintLinearPeriodicBoundaryShapeFunctions.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupFineScaleDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackOpening2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsMultiscale2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
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
    mShiftCenterDamage[0] = 0.0;
    mShiftCenterDamage[1] = 0.;
    mCrackOpening[0] = 0.0;
    mCrackOpening[1] = 0.0;
    mDOFCrackOpening[0] = -1;
    mDOFCrackOpening[1] = -1;
    mEpsilonHom.mEngineeringStrain[0] = 0.;
    mEpsilonHom.mEngineeringStrain[1] = 0.;
    mEpsilonHom.mEngineeringStrain[2] = 0.;
    mScalingFactorCrackOpening = 1e-1;
    mScalingFactorEpsilon = 1e-5;
    mDOFGlobalTotalStrain[0] = -1;
    mDOFGlobalTotalStrain[1] = -1;
    mDOFGlobalTotalStrain[2] = -1;
    mDOFPeriodicBoundaryDisplacements[0] = -1.;
    mDOFPeriodicBoundaryDisplacements[1] = -1.;
    mDOFPeriodicBoundaryDisplacements[2] = -1.;
    mPeriodicBoundaryDisplacements[0] = 0.;
    mPeriodicBoundaryDisplacements[1] = 0.;
    mPeriodicBoundaryDisplacements[2] = 0.;
    mConstraintFineScaleDamageX = -1;
    mConstraintFineScaleDamageY = -1;
    mConstraintFineScaleHomogeneousX = -1;
    mConstraintFineScaleHomogeneousY = -1;
    mConstraintFineScalePeriodicDamage = -1;
    mConstraintFineScalePeriodicHomogeneous = -1;
    mConstraintNormalCrackOpening = -1;
    mConstraintTangentialCrackOpening = -1;
    mConstraintPeriodicBoundaryShapeFunctions[0] = -1;
    mConstraintPeriodicBoundaryShapeFunctions[1] = -1;
    mConstraintPeriodicBoundaryShapeFunctions[2] = -1;
    mlCoarseScaleCrack = 0.;
    mlFineScaleCrack = 0.;
    mlFineScale = 0.;
    mCenterDamage[0] = 0.;
    mCenterDamage[1] = 0.;
    mCenterHomogeneous[0] = 0.;
    mCenterHomogeneous[1] = 0.;
    mCenterMacro[0] = 0.;
    mCenterMacro[1] = 0.;
    mScaleWRTDamageCenter = true;
    mFineScaleArea = 0.;
    mCoarseScaleArea = 0.;
    mSquareFineScaleModel = true;
    mCrackTransitionRadius = 0;
    mIPName = std::string("fineScaleIp");
    mResultDirectory = std::string(".");
    mLoadStepMacro = std::string("loadstep_1");
    mGroupBoundaryNodesDamage = -1;
    mGroupBoundaryNodesHomogeneous = -1;
    mGroupElementsDamage = -1;
    mGroupElementsHomogeneous = -1;
    mBoundaryNodesElementsAssigned = false;
    mThickness = 1.;

    // generate the constrain equation
    mConstraintTotalStrain = CreateConstraintLinearGlobalTotalStrain();

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
}

NuTo::StructureMultiscale::StructureMultiscale()
{

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
    mLogger << "start serialization of StructureMultiscale" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Structure)
    & BOOST_SERIALIZATION_NVP(mCrackAngle)
    & BOOST_SERIALIZATION_NVP(mShiftCenterDamage)
    & BOOST_SERIALIZATION_NVP(mCrackOpening)
    & BOOST_SERIALIZATION_NVP(mDOFCrackOpening)
    & BOOST_SERIALIZATION_NVP(mEpsilonTot)
    & BOOST_SERIALIZATION_NVP(mEpsilonTotConstraint)
    & BOOST_SERIALIZATION_NVP(mDOFGlobalTotalStrain)
    & BOOST_SERIALIZATION_NVP(mDOFPeriodicBoundaryDisplacements)
    & BOOST_SERIALIZATION_NVP(mPeriodicBoundaryDisplacements)
    & BOOST_SERIALIZATION_NVP(mEpsilonHom)
    & BOOST_SERIALIZATION_NVP(mScalingFactorCrackOpening)
    & BOOST_SERIALIZATION_NVP(mScalingFactorEpsilon)
    & BOOST_SERIALIZATION_NVP(mlCoarseScaleCrack)
    & BOOST_SERIALIZATION_NVP(mlFineScaleCrack)
    & BOOST_SERIALIZATION_NVP(mlFineScale)
    & BOOST_SERIALIZATION_NVP(mCrackTransitionRadius)
    & BOOST_SERIALIZATION_NVP(mCenterDamage)
    & BOOST_SERIALIZATION_NVP(mCenterHomogeneous)
    & BOOST_SERIALIZATION_NVP(mCenterMacro)
    & BOOST_SERIALIZATION_NVP(mScaleWRTDamageCenter)
    & BOOST_SERIALIZATION_NVP(mFineScaleArea)
    & BOOST_SERIALIZATION_NVP(mCoarseScaleArea)
    & BOOST_SERIALIZATION_NVP(mSquareFineScaleModel)
    & BOOST_SERIALIZATION_NVP(mThickness)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScaleDamageX)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScaleDamageY)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScaleHomogeneousX)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScaleHomogeneousY)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScalePeriodicDamage)
    & BOOST_SERIALIZATION_NVP(mConstraintFineScalePeriodicHomogeneous)
    & BOOST_SERIALIZATION_NVP(mConstraintNormalCrackOpening)
    & BOOST_SERIALIZATION_NVP(mConstraintTangentialCrackOpening)
    & BOOST_SERIALIZATION_NVP(mConstraintTotalStrain)
    & BOOST_SERIALIZATION_NVP(mConstraintPeriodicBoundaryShapeFunctions)
    & BOOST_SERIALIZATION_NVP(mBoundaryNodesElementsAssigned)
    & BOOST_SERIALIZATION_NVP(mGroupBoundaryNodesDamage)
    & BOOST_SERIALIZATION_NVP(mGroupBoundaryNodesHomogeneous)
    & BOOST_SERIALIZATION_NVP(mGroupNodesDamage)
    & BOOST_SERIALIZATION_NVP(mGroupNodesHomogeneous)
    & BOOST_SERIALIZATION_NVP(mGroupElementsDamage)
    & BOOST_SERIALIZATION_NVP(mGroupElementsHomogeneous)
    & BOOST_SERIALIZATION_NVP(mIPName)
    & BOOST_SERIALIZATION_NVP(mResultDirectory)
    & BOOST_SERIALIZATION_NVP(mLoadStepMacro)
    & BOOST_SERIALIZATION_NVP(mPrevEpsilonTot)
    & BOOST_SERIALIZATION_NVP(mDeltaEpsilonTot)
    ;
#ifdef DEBUG_SERIALIZATION
    mLogger << "finish serialization of StructureMultiscale" << "\n";
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
        mLogger<<"[NuTo::Structure::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
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
        mLogger<<"[NuTo::StructureMultiscale::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
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
    if (GroupGetNumMembers(mGroupElementsHomogeneous)!=GetNumElements())
    {
        mLogger << " homogeneous elements " << GroupGetNumMembers(mGroupElementsHomogeneous) << " total nodes" << GetNumElements() << "\n";
        throw MechanicsException("[NuTo::StructureMultiscale::TransformMultiscaleNodes] there is something wrong with your elements groups");
    }
    if (GroupGetNumMembers(mGroupElementsHomogeneous)!=GetNumElements())
        throw MechanicsException("[NuTo::StructureMultiscale::TransformMultiscaleNodes] there is something wrong with your elements");
    TransformMultiscaleNodes(mGroupBoundaryNodesHomogeneous, mGroupBoundaryNodesHomogeneous, mGroupElementsHomogeneous, mCenterHomogeneous, mFineScaleArea,false);
}

//! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
//! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
void NuTo::StructureMultiscale::TransformMultiscaleNodes(int rGroupBoundaryNodes, int rGroupNodes, int rGroupElements, boost::array<double,2>& rCenter, double& rArea, bool rCrackedDomain)
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
    itBoundaryNode->second->GetCoordinates2D(coord1);
    itBoundaryNode++;
    itBoundaryNode->second->GetCoordinates2D(coord2);
    itBoundaryNode++;
    itBoundaryNode->second->GetCoordinates2D(coord3);
    //check for being on a straight line
    bool onLine(true);
    while (onLine)
    {
        vec1[0] = coord2[0]-coord1[0];
        vec1[1] = coord2[1]-coord1[1];
        vec2[0] = coord3[0]-coord1[0];
        vec2[1] = coord3[1]-coord1[1];
        mLogger << (vec1[0]*vec2[0]+vec2[1]*vec1[1])/(sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1])*sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1])) << "\n";
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
    if (nodeGroupBoundary->GetNumMembers()==4)
        mSquareFineScaleModel = true;
    else
        mSquareFineScaleModel = false;
    for (Group<NodeBase>::iterator itNode=nodeGroupBoundary->begin(); itNode!=nodeGroupBoundary->end(); itNode++)
    {
        double coordinates[2];
        itNode->second->GetCoordinates2D(coordinates);
        if (mSquareFineScaleModel==false)
        {
            if (fabs((rCenter[0]-coordinates[0])*(rCenter[0]-coordinates[0])+(rCenter[1]-coordinates[1])*(rCenter[1]-coordinates[1])-fineScaleRadius*fineScaleRadius)>0.001*fineScaleRadius*fineScaleRadius)
            {
                mLogger << (rCenter[0]-coordinates[0])*(rCenter[0]-coordinates[0])+(rCenter[1]-coordinates[1])*(rCenter[1]-coordinates[1]) << " " << fineScaleRadius*fineScaleRadius << "\n";
                mSquareFineScaleModel = true;
                break;
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

        switch (nodeType)
        {
        case Node::NodeCoordinatesDisplacements2D:
            newNode = new NodeCoordinatesDisplacementsMultiscale2D(this, rCrackedDomain);
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
    //NuTo::FullMatrix<double> direction(2,1);
    //direction(0,0) = 1.;
    //direction(1,0) = 0.;
    //mConstraintFineScaleX = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    //direction(0,0) = 0.;
    //direction(1,0) = 1.;
    //mConstraintFineScaleY = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup, direction, 0);

    if (minX==maxX || minY==maxY)
        throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] structure has zero width or height, either check your boundary group or the full structure.");
    mlFineScale = maxX-minX;
    if (mSquareFineScaleModel)
    {
        if (fabs((maxX-minX)-(maxY-minY))>1e-3)
            throw MechanicsException("[NuTo::StructureMultiscale::TransformBoundaryNodes] domain is not square.");
        rCenter[0] = 0.5*(maxX+minX);
        rCenter[1] = 0.5*(maxY+minY);
        rArea = (maxX-minX)* (maxY-minY);

        mLogger<< "square fine scale model" << "\n";
    }
    else
    {
        /* this only works for problems without holes
        rArea = 0;
        std::cout << "number of elements " << mElementMap.size() << std::endl;
        boost::ptr_map<int,ElementBase>::iterator elementIter;
        for (elementIter=mElementMap.begin(); elementIter!=mElementMap.end(); elementIter++)
        {
        	std::vector<double> volume;
        	elementIter->second->GetIntegrationPointVolume(volume);
        	for (unsigned int count=0; count<volume.size(); count++)
        	{
        		rArea+=volume[count];
        	}
        }*/

        //this area calculation takes into account the missing area between the circle and the secants
        if (mElementMap.begin()->second->GetNumNodes()==6)
            rArea = fineScaleRadius*fineScaleRadius*M_PI;
        else
        {
            if (mElementMap.begin()->second->GetNumNodes()==3)
                rArea = nodeGroupBoundary->GetNumMembers()*fineScaleRadius*fineScaleRadius*0.5*sin(2.*M_PI/nodeGroupBoundary->GetNumMembers());
            else
                throw MechanicsException("[NuTo::StructureMultiscale::TransformMultiscaleNodes] the area of the surrounding boundary is only calculated for linear and quadratric tringular meshes.");
        }
        mLogger << "round fine scale model" << "\n";
    }
    mLogger<< "center of fine scale at (" << rCenter[0] << "," << rCenter[1] << "), area(" << rArea << ")" << "\n";
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureMultiscale::TransformBoundaryNodes] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief numbers non standard DOFs' e.g. in StructureMultiscale, for standard structures this routine is empty
void NuTo::StructureMultiscale::NumberAdditionalGlobalDofs()
{
    // DOFs related to the crack angle and crack opening
    if (mDimension==2)
    {
        mDOFCrackOpening[0] = mNumDofs++;
        mDOFCrackOpening[1] = mNumDofs++;
        mDOFGlobalTotalStrain[0] = mNumDofs++;
        mDOFGlobalTotalStrain[1] = mNumDofs++;
        mDOFGlobalTotalStrain[2] = mNumDofs++;
        mDOFPeriodicBoundaryDisplacements[0] = mNumDofs++;
        mDOFPeriodicBoundaryDisplacements[1] = mNumDofs++;
        mDOFPeriodicBoundaryDisplacements[2] = mNumDofs++;
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
            value*=mScalingFactorCrackOpening;
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
            value*=mScalingFactorEpsilon;
            this->mEpsilonTot.mEngineeringStrain[count] = value;
        }

        //merge global displacements for periodic bc
        for (int count=0; count<3; count++)
        {
            int dof = this->mDOFPeriodicBoundaryDisplacements[count];
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
            this->mPeriodicBoundaryDisplacements[count] = value;
        }
    }
    else
        throw MechanicsException("[NuTo::StructureMultiscale::NodeMergeActiveDofValues] Only implemented for 2D");

    //update the homogeneous strain and the elastic crack angle
    CalculateHomogeneousEngineeringStrain();
}

void NuTo::StructureMultiscale::SetCrackOpening(NuTo::FullMatrix<double>& crackOpening)
{
    if (crackOpening.GetNumRows()!=2 || crackOpening.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::StructureMultiscale::SetCrackOpening] the matrix has to be of dimension (2,1).");
    this->mCrackOpening[0] = crackOpening(0,0);
    this->mCrackOpening[1] = crackOpening(1,0);

    //update the homogeneous strain and the elastic crack angle
    CalculateHomogeneousEngineeringStrain();
}

//! @brief extract dof values additional dof values
//! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
//! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
void NuTo::StructureMultiscale::NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const
{
    if (mDimension==2)
    {
        for (int count=0; count<2; count++)
        {
            int dof = this->mDOFCrackOpening[count];
            double value = this->mCrackOpening[count]/mScalingFactorCrackOpening;
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
            double value = this->mEpsilonTot.mEngineeringStrain[count]/mScalingFactorEpsilon;
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
            int dof = this->mDOFPeriodicBoundaryDisplacements[count];
            double value = this->mPeriodicBoundaryDisplacements[count];
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

NuTo::Error::eError NuTo::StructureMultiscale::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
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
    std::vector<boost::array<double,5> > dDOF;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF);
    double bURow[2], bMRow[3];
    Error::eError errorGlobal(Error::SUCCESSFUL);

    if (mGroupElementsDamage!=-1)
    {
        // loop over all elements in the damaged region
        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Elements)
            throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
        const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
        assert(elementGroup!=0);
        double scalingFactorDamage = this->GetScalingFactorDamage();
#ifdef _OPENMP
    	if (mNumProcessors!=0)
    		omp_set_num_threads(mNumProcessors);
		#pragma omp parallel default(shared) private(elementVector,elementVectorGlobalDofs)
#endif
        for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
        {
#ifdef _OPENMP
			#pragma omp single nowait
#endif
			{
				// calculate element contribution
				Error::eError error = itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
				if (error!=Error::SUCCESSFUL)
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
				}
				else
				{
					elementVector*=scalingFactorDamage;
					assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());
					assert(static_cast<unsigned int>(elementVector.GetNumColumns()) == 1);

					// write element contribution to global vectors
					#pragma omp critical (StructureMultiscaleBuildGlobalGradientInternalPotentialSubVectorsDamage)
					{
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
								//bAlphaRow =dDOF[theDofMapRow][0];

								//influence of u
								bURow[0] = dDOF[theDofMapRow][0];
								bURow[1] = dDOF[theDofMapRow][1];

								//influence of epsilon_tot = epsilon_hom
								bMRow[0] = dDOF[theDofMapRow][2];
								bMRow[1] = dDOF[theDofMapRow][3];
								bMRow[2] = dDOF[theDofMapRow][4];

								/*                if (mDOFCrackAngle< this->mNumActiveDofs)
												{
													rActiveDofGradientVector(mDOFCrackAngle,0) += bAlphaRow * elementVector(rowCount,0);
													//std::cout << "add " << bAlphaRow << "*" << elementVector(rowCount,0) << "=" << rActiveDofGradientVector(mDOFCrackAngle,0) << "(DOF "<< mDOFCrackAngle<<")" <<std::endl;
												}
												else
													rDependentDofGradientVector(mDOFCrackAngle - this->mNumActiveDofs,0) += bAlphaRow * elementVector(rowCount,0);
								*/
								if (mDOFCrackOpening[0]< this->mNumActiveDofs)
								{
									rActiveDofGradientVector(mDOFCrackOpening[0],0) += bURow[0] * elementVector(rowCount,0);
									//std::cout << "add " << bURow[0] << "*" << elementVector(rowCount,0) << "=" << rActiveDofGradientVector(mDOFCrackOpening[0],0) << "(DOF "<< globalRowDof<<","<<mDOFCrackOpening[0]<<")" <<std::endl;
								}
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
                }
            }
        }
    }
    // loop over all elements in the homogeneous region
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) is not an element group.");
    const Group<ElementBase>* elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
#endif
	#pragma omp parallel default(shared) private(elementVector,elementVectorGlobalDofs)
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
    {
		#pragma omp single nowait
		{
//			std::cout << "thread " << omp_get_thread_num( ) << "\n";
			// calculate element contribution
			Error::eError error = itElement->second->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
			if (error!=Error::SUCCESSFUL)
			{
				if (errorGlobal==Error::SUCCESSFUL)
					errorGlobal = error;
				else if (errorGlobal!=error)
					throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
			}
			else
			{
				elementVector*=scalingFactorHomogeneous;
				assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());
				assert(static_cast<unsigned int>(elementVector.GetNumColumns()) == 1);

				// write element contribution to global vectors
				#pragma omp critical (StructureMultiscaleBuildGlobalGradientInternalPotentialSubVectorsHomogeneous)
				{
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
							//bAlphaRow =dDOF[theDofMapRow][0];

							//influence of u
							bURow[0] = dDOF[theDofMapRow][0];;
							bURow[1] = dDOF[theDofMapRow][1];;

							//influence of epsilon_tot = epsilon_hom
							bMRow[0] = dDOF[theDofMapRow][2];;
							bMRow[1] = dDOF[theDofMapRow][3];;
							bMRow[2] = dDOF[theDofMapRow][4];;

							/*                if (mDOFCrackAngle< this->mNumActiveDofs)
											{
												rActiveDofGradientVector(mDOFCrackAngle,0) += bAlphaRow * elementVector(rowCount,0);
												//std::cout << "add " << bAlphaRow << "*" << elementVector(rowCount,0) << "=" << rActiveDofGradientVector(mDOFCrackAngle,0) << "(DOF "<< mDOFCrackAngle<<")" <<std::endl;
											}
											else
												rDependentDofGradientVector(mDOFCrackAngle - this->mNumActiveDofs,0) += bAlphaRow * elementVector(rowCount,0);
							*/
							if (mDOFCrackOpening[0]< this->mNumActiveDofs)
							{
								rActiveDofGradientVector(mDOFCrackOpening[0],0) += bURow[0] * elementVector(rowCount,0);
								//std::cout << "add " << bURow[0] << "*" << elementVector(rowCount,0) << "=" << rActiveDofGradientVector(mDOFCrackOpening[0],0) << "(DOF "<< globalRowDof<<","<<mDOFCrackOpening[0]<<")" <<std::endl;
							}
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
			}
        }
    }
    /*std::cout << "without constraint rActiveDofGradientVector "<< "\n";
    rActiveDofGradientVector.Trans().Info(12,10);
    std::cout << "without constraint rDependentDofGradientVector "<< "\n";
    rDependentDofGradientVector.Trans().Info(12,10);
    */
    //write contribution of Lagrange Multipliers

    ConstraintBuildGlobalGradientInternalPotentialSubVectors(rActiveDofGradientVector, rDependentDofGradientVector);
    /*    std::cout << "with constraint rActiveDofGradientVector "<< "\n";
        rActiveDofGradientVector.Trans().Info(12,10);
        std::cout << "with constraint rDependentDofGradientVector "<< "\n";
        rDependentDofGradientVector.Trans().Info(12,10);
    */
    return errorGlobal;
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
NuTo::Error::eError NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
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
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<boost::array<double,5> > dDOF;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF);

    //mLogger << "alpha "  << mCrackAngle << " crack opening "<< mCrackOpening[0] << " " << mCrackOpening[1] <<
    //		     " epsilon hom" << mEpsilonHom.GetData()[0] << " "<< mEpsilonHom.GetData()[1]<< " "<< mEpsilonHom.GetData()[2]<<"\n";
    Error::eError errorGlobal(Error::SUCCESSFUL);
    if (mGroupElementsDamage!=-1)
    {
        // loop over all elements in the damaged region
        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Elements)
            throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
        const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
        assert(elementGroup!=0);
        double scalingFactorDamage = this->GetScalingFactorDamage();
#ifdef _OPENMP
    	if (mNumProcessors!=0)
    		omp_set_num_threads(mNumProcessors);
        #pragma omp parallel default(shared) private(elementMatrix,elementMatrixGlobalDofsRow,elementMatrixGlobalDofsColumn)
#endif
        for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
        {
			#ifdef _OPENMP
			#pragma omp single nowait
			#endif
			{
				// calculate element contribution
				bool symmetryFlag = false;
				Error::eError error = itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
				if (error!=Error::SUCCESSFUL)
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
				}
				else
				{
					elementMatrix*=scalingFactorDamage;

					assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
					assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

					#ifdef _OPENMP
					#pragma omp critical (StructureMultiscaleBuildGlobalCoefficientSubMatrices0GeneralDamage)
					#endif
					this->AddElementMatrixToGlobalSubMatricesGeneral(
						elementMatrix,
						elementMatrixGlobalDofsRow,
						elementMatrixGlobalDofsColumn,
						mappingDofMultiscaleNode,
						dDOF,
						&rMatrixJJ,
						&rMatrixJK,
						0,
						0,
						false);
				}
			}
        }
    }

    // loop over all elements in the homogeneous region
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(homogeneous) is not an element group.");
    const Group<ElementBase>* elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
	#pragma omp parallel default(shared) private(elementMatrix,elementMatrixGlobalDofsRow,elementMatrixGlobalDofsColumn)
#endif
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
    {
#ifdef _OPENMP
		#pragma omp single nowait
#endif
	    {
			// calculate element contribution
			bool symmetryFlag = false;
			Error::eError error = itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
			if (error!=Error::SUCCESSFUL)
			{
				if (errorGlobal==Error::SUCCESSFUL)
					errorGlobal = error;
				else if (errorGlobal!=error)
					throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
			}
			else
			{
				elementMatrix*=scalingFactorHomogeneous;

				assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
				assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

#ifdef _OPENMP
				#pragma omp critical (StructureMultiscaleBuildGlobalCoefficientSubMatrices0GeneralHomogeneous)
#endif
				this->AddElementMatrixToGlobalSubMatricesGeneral(
					elementMatrix,
					elementMatrixGlobalDofsRow,
					elementMatrixGlobalDofsColumn,
					mappingDofMultiscaleNode,
					dDOF,
					&rMatrixJJ,
					&rMatrixJK,
					0,
					0,
					false);
			}
	    }
    }


    /*
        {
            mLogger << "JJ before "<< "\n";
            NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
            rMatrixJJFull.Info(12,5);
            mLogger << "JK before "<< "\n";
            NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
            rMatrixJKFull.Info(12,5);

        }
    */
    //add contribution of Lagrange Multipliers
    ConstraintsBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK);
    /*    {
            mLogger << "JJ after "<< "\n";
            NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
            rMatrixJJFull.Info(12,5);
            mLogger << "JK after "<< "\n";
            NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
            rMatrixJKFull.Info(12,5);

        }
    */
    return errorGlobal;
}

// based on the global dofs build submatrices of the global coefficent matrix0
NuTo::Error::eError NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK) const
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
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;
    std::vector<int> elementVectorGlobalDofs;

    // calculate for all multiscale dofs the derivatives
    // with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<boost::array<double,5> > dDOF;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF);

    Error::eError errorGlobal(Error::SUCCESSFUL);

    if (mGroupElementsDamage!=-1)
    {
        // loop over all elements in the damaged region
        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Elements)
            throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(damage) is not an element group.");
        const Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
        assert(elementGroup!=0);
        double scalingFactorDamage = this->GetScalingFactorDamage();
#ifdef _OPENMP
    	if (mNumProcessors!=0)
    		omp_set_num_threads(mNumProcessors);
        #pragma omp parallel default(shared) private(elementMatrix,elementMatrixGlobalDofsRow,elementMatrixGlobalDofsColumn)
#endif
        for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
        {
#ifdef _OPENMP
			#pragma omp single nowait
#endif
			{
				// calculate element contribution
				bool symmetryFlag = false;
				Error::eError error = itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
				if (error!=Error::SUCCESSFUL)
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
				}
				else
				{
					elementMatrix*=scalingFactorDamage;

					assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
					assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

#ifdef _OPENMP
					#pragma omp critical (StructureMultiscaleBuildGlobalCoefficientSubMatrices0GeneralDamage)
#endif
					this->AddElementMatrixToGlobalSubMatricesGeneral(
						elementMatrix,
						elementMatrixGlobalDofsRow,
						elementMatrixGlobalDofsColumn,
						mappingDofMultiscaleNode,
						dDOF,
						&rMatrixJJ,
						&rMatrixJK,
						&rMatrixKJ,
						&rMatrixKK,
						true);
				}
			}
        }
    }

    // loop over all elements in the homogeneous region
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupElementsHomogeneous);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0General] Group(Homogeneous) is not an element group.");
    const Group<ElementBase>* elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);
    double scalingFactorHomogeneous = this->GetScalingFactorHomogeneous();
#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared) private(elementMatrix,elementMatrixGlobalDofsRow,elementMatrixGlobalDofsColumn)
#endif
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end(); itElement++)
    {
#ifdef _OPENMP
		#pragma omp single nowait
#endif
		{
			// calculate element contribution
			bool symmetryFlag = false;
			Error::eError error = itElement->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
			if (error!=Error::SUCCESSFUL)
			{
				if (errorGlobal==Error::SUCCESSFUL)
					errorGlobal = error;
				else if (errorGlobal!=error)
					throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
			}
			else
			{
				elementMatrix*=scalingFactorHomogeneous;

				assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
				assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());

#ifdef _OPENMP
				#pragma omp critical (StructureMultiscaleBuildGlobalCoefficientSubMatrices0GeneralHomogeneous)
#endif
				this->AddElementMatrixToGlobalSubMatricesGeneral(
					elementMatrix,
					elementMatrixGlobalDofsRow,
					elementMatrixGlobalDofsColumn,
					mappingDofMultiscaleNode,
					dDOF,
					&rMatrixJJ,
					&rMatrixJK,
					&rMatrixKJ,
					&rMatrixKK,
					true);
			}
		}
    }

    /*
    {
        mLogger << "JJ before "<< "\n";
        NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
        rMatrixJJFull.Info(12,5);
        mLogger << "JK before "<< "\n";
        NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
        rMatrixJKFull.Info(12,5);
        mLogger << "JJ before "<< "\n";
        NuTo::FullMatrix<double> rMatrixKJFull(rMatrixKJ);
        rMatrixKJFull.Info(12,5);
        mLogger << "KK before "<< "\n";
        NuTo::FullMatrix<double> rMatrixKKFull(rMatrixKK);
        rMatrixKKFull.Info(12,5);

    }
    */
    //add contribution of Lagrange Multipliers
    ConstraintBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK, rMatrixKJ, rMatrixKK);
    /*
        {
            mLogger << "JJ after "<< "\n";
            NuTo::FullMatrix<double> rMatrixJJFull(rMatrixJJ);
            rMatrixJJFull.Info(12,5);
            mLogger << "JK after "<< "\n";
            NuTo::FullMatrix<double> rMatrixJKFull(rMatrixJK);
            rMatrixJKFull.Info(12,5);
            mLogger << "JK after "<< "\n";
            NuTo::FullMatrix<double> rMatrixKJFull(rMatrixKJ);
            rMatrixKJFull.Info(12,5);
            mLogger << "KK after "<< "\n";
            NuTo::FullMatrix<double> rMatrixKKFull(rMatrixKK);
            rMatrixKKFull.Info(12,5);
        }
    */
    return errorGlobal;
}

void NuTo::StructureMultiscale::AddElementMatrixToGlobalSubMatricesGeneral(
    NuTo::FullMatrix<double>& elementMatrix,
    std::vector<int>& elementMatrixGlobalDofsRow,
    std::vector<int>& elementMatrixGlobalDofsColumn,
    std::vector<int>& mappingDofMultiscaleNode,
    std::vector<boost::array<double,5> >&dDOF,
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
    double bURow[2]= { 0, 0}, bMRow[3]= {0,0,0}, bUCol[2]= {0,0}, bMCol[3]= {0,0,0};

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
            //bAlphaRow =dDOF[theDofMapRow][0];

            //influence of u
            bURow[0] = dDOF[theDofMapRow][0];
            bURow[1] = dDOF[theDofMapRow][1];

            //influence of e_hom = e_tot
            bMRow[0] = dDOF[theDofMapRow][2];
            bMRow[1] = dDOF[theDofMapRow][3];
            bMRow[2] = dDOF[theDofMapRow][4];
        }
        for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
        {
            int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
            if (mappingDofMultiscaleNode[globalColumnDof]!=-1)
            {
                //include influence of global dofs like crack orientation (alpha) and crackopening (ux,uy)
                //influence of alpha
                int theDofMapCol = mappingDofMultiscaleNode[globalColumnDof];
                //bAlphaCol = dDOF[theDofMapCol][0];

                //influence of u
                bUCol[0] = dDOF[theDofMapCol][0];
                bUCol[1] = dDOF[theDofMapCol][1];

                //influence of e_hom and e_tot
                bMCol[0] = dDOF[theDofMapCol][2];
                bMCol[1] = dDOF[theDofMapCol][3];
                bMCol[2] = dDOF[theDofMapCol][4];

                if (globalRowDof < this->mNumActiveDofs)
                {
                    /*                    if (mDOFCrackAngle < this->mNumActiveDofs)
                                        {
                                            rMatrixJJ->AddEntry(globalRowDof, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                                        }
                                        else
                                        {
                                            rMatrixJK->AddEntry(globalRowDof, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                                        }
                    */
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
                        /*                        if (mDOFCrackAngle < this->mNumActiveDofs)
                                                {
                                                    rMatrixKJ->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                                                }
                                                else
                                                {
                                                    rMatrixKK->AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                                                }
                        */
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
                    /*                    if (mDOFCrackAngle < this->mNumActiveDofs)
                                        {
                                            rMatrixJJ->AddEntry(mDOFCrackAngle,      mDOFCrackAngle,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                                            {
                                                rMatrixJJ->AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                                                //mLogger << "add at ("<< mDOFCrackAngle << "," << mDOFCrackOpening[0]<< ") " << bAlphaRow<< "(" << bAlphaRow << ")*elementMatrix(" << elementMatrix(rowCount, colCount) << ")*bUCol(" <<  bUCol[0] << ")=" << bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0] << "\n";
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
                    */
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        /*                        if (mDOFCrackAngle < this->mNumActiveDofs)
                                                {
                                                    rMatrixJJ->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                }
                                                else
                                                {
                                                    rMatrixJK->AddEntry(mDOFCrackOpening[0], mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                }
                        */
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
                            /*                            if (mDOFCrackAngle < this->mNumActiveDofs)
                                                        {
                                                            rMatrixKJ->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                        }
                                                        else
                                                        {
                                                            rMatrixKK->AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                        }
                            */
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
                        /*                        if (mDOFCrackAngle < this->mNumActiveDofs)
                                                {
                                                    rMatrixJJ->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                }
                                                else
                                                {
                                                    rMatrixJK->AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                                                }
                        */
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
                            /*                            if (mDOFCrackAngle < this->mNumActiveDofs)
                                                        {
                                                            rMatrixKJ->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                                                        }
                                                        else
                                                        {
                                                            rMatrixKK->AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                                                        }
                            */
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
                        /*                       if (mDOFCrackAngle < this->mNumActiveDofs)
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
                        */
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
                            /*                           if (mDOFCrackAngle < this->mNumActiveDofs)
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
                            */
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
                    /*                    if (mDOFCrackAngle < this->mNumActiveDofs)
                                        {
                                            rMatrixJJ->AddEntry(mDOFCrackAngle,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                                        }
                                        else
                                        {
                                            if (rCalcMatrixKJ_KK)
                                                rMatrixKJ->AddEntry(mDOFCrackAngle - this->mNumActiveDofs,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                                        }
                    */
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
                    /*                    if (mDOFCrackAngle < this->mNumActiveDofs)
                                        {
                                            rMatrixJK->AddEntry(mDOFCrackAngle, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                                        }
                                        else
                                        {
                                            if (rCalcMatrixKJ_KK)
                                                rMatrixKK->AddEntry(mDOFCrackAngle - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                                        }
                    */
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
//! @param rDOF return value, for each dof, the corresponding derivatives ux, uy, exx, eyy, gamma_xy  [0..5]
//! @param r2DOF return value, for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy)
void NuTo::StructureMultiscale::CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<boost::array<double,5> >& rDOF)const
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
    double bHomU[6]; //bHomU[0-2] for ux [3-5] for uy
    GetdEpsilonHomdCrack(bHomU);

    //check the derivatives
#ifdef MYDEBUG
    double interval(1e-8);
    EngineeringStrain2D strain1 = mEpsilonHom;
    mLogger << "algo alpha " << bHomAlpha[0] << "  " << bHomAlpha[1] << "  " <<  bHomAlpha[2] << "\n";
    const_cast<StructureMultiscale*>(&*this)->mCrackAngle += interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    EngineeringStrain2D strain2 = mEpsilonHom;
    mLogger << "cdf alpha " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
            << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
            << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< "\n";
    const_cast<StructureMultiscale*>(&*this)->mCrackAngle -= interval;
    const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();

    for (int count=0; count<2; count++)
    {
        mLogger << "algo u" << bHomU[0+3*count] << "  " << bHomU[1+3*count] << "  " <<  bHomU[2+3*count] << "\n";
        strain1 = mEpsilonHom;
        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count]+=interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
        EngineeringStrain2D strain2 = mEpsilonHom;
        mLogger << "cdf " << (strain2.mEngineeringStrain[0]-strain1.mEngineeringStrain[0])/interval << "  "
                << (strain2.mEngineeringStrain[1]-strain1.mEngineeringStrain[1])/interval << "  "
                << (strain2.mEngineeringStrain[2]-strain1.mEngineeringStrain[2])/interval<< "\n";
        const_cast<StructureMultiscale*>(&*this)->mCrackOpening[count] -= interval;
        const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
    }
    mLogger << "crack opening1 " << mCrackOpening[0] << " " << mCrackOpening[1]<< "\n";
    //check hessian
    double hessianCDF[9],bHomAlpha2[3],bHomU2[6];
    mLogger << "algo hessian " << "\n";
    mLogger << bHessian[0] << "  " << bHessian[3] << "  " <<  bHessian[6] << "\n";
    mLogger << bHessian[1] << "  " << bHessian[4] << "  " <<  bHessian[7] << "\n";
    mLogger << bHessian[2] << "  " << bHessian[5] << "  " <<  bHessian[8] << "\n";
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
    mLogger << "cdf hessian " << "\n";
    mLogger << hessianCDF[0] << "  " << hessianCDF[3] << "  " <<  hessianCDF[6] << "\n";
    mLogger << hessianCDF[1] << "  " << hessianCDF[4] << "  " <<  hessianCDF[7] << "\n";
    mLogger << hessianCDF[2] << "  " << hessianCDF[5] << "  " <<  hessianCDF[8] << "\n";

    if (fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 ||
            fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 ||
            fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 || fabs(hessianCDF[0]-bHessian[0])>1e-3 )
    {
        mLogger << "[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] hessian is not correct" << "\n";
        throw MechanicsException("CalculatedDispdGlobalDofs hessian is not correct.");
    }

    mLogger << "\n";
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
                mLogger << "dDofxdEpsilonHom algo " << dDOFdEpsilonHomX[0] << "  " << dDOFdEpsilonHomX[1] << "  " << dDOFdEpsilonHomX[2] << "\n";
                mLogger << "dDofydEpsilonHom algo " << dDOFdEpsilonHomY[0] << "  " << dDOFdEpsilonHomY[1] << "  " << dDOFdEpsilonHomY[2] << "\n";

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
                    if (fabs(dDOFdEpsilonHomX[count]-dDOFdEpsilonHomXcdf[count])>1e-2 || fabs(dDOFdEpsilonHomY[count]-dDOFdEpsilonHomYcdf[count])>1e-2)
                    {
                        mLogger << "[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDOFdEpsilonHomX/Y is not correct" << "\n";
                        throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDOFdEpsilonHomX/Y is not correct.");
                    }
                }
                mLogger << "dDofxdEpsilonHom cdf " << dDOFdEpsilonHomXcdf[0] << "  " << dDOFdEpsilonHomXcdf[1] << "  " << dDOFdEpsilonHomXcdf[2] << "\n";
                mLogger << "dDofydEpsilonHom cdf " << dDOFdEpsilonHomYcdf[0] << "  " << dDOFdEpsilonHomYcdf[1] << "  " << dDOFdEpsilonHomYcdf[2] << "\n";
                mLogger << "\n";
            }
#endif

            //derivative of displacement with respect to discontinuity (crack opening)
            if (nodeInCrackedDomain)
            {
                GetdDisplacementdCrackOpening(coord,&(rDOF[countMultiscaleDofs][0]),&(rDOF[countMultiscaleDofs+1][0]));
            }
            else
            {
                rDOF[countMultiscaleDofs][0] = 0.;
                rDOF[countMultiscaleDofs+1][0] = 0.;
                rDOF[countMultiscaleDofs][1] = 0.;
                rDOF[countMultiscaleDofs+1][1] = 0.;
            }

            rDOF[countMultiscaleDofs][0]  += dDOFdEpsilonHomX[0]*bHomU[0] +
                                             dDOFdEpsilonHomX[1]*bHomU[1] +
                                             dDOFdEpsilonHomX[2]*bHomU[2];

            rDOF[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHomU[0] +
                                              dDOFdEpsilonHomY[1]*bHomU[1] +
                                              dDOFdEpsilonHomY[2]*bHomU[2];

#ifdef MYDEBUG
            {
                double interval(1e-5);
                double disp1[2],disp2[2],disptmp[2];
                mLogger << "dDofdU1 algo " << rDOF[countMultiscaleDofs][1] << "  " << rDOF[countMultiscaleDofs+1][1] << "\n";
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
                mLogger << "dDofdU1 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << "\n";
                if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][1])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][1])>1e-2)
                {
                    mLogger << "[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDofdU1 is not correct" << "\n";
                    throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDofdU1 is not correct");
                }
                const_cast<StructureMultiscale*>(&*this)->mCrackOpening[0] -= interval;
                const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
                mLogger << "\n";
            }
#endif
            //scale the derivative
            rDOF[countMultiscaleDofs][0]   *= mScalingFactorCrackOpening;
            rDOF[countMultiscaleDofs+1][0] *= mScalingFactorCrackOpening;

            rDOF[countMultiscaleDofs][1]  += dDOFdEpsilonHomX[0]*bHomU[3] +
                                             dDOFdEpsilonHomX[1]*bHomU[4] +
                                             dDOFdEpsilonHomX[2]*bHomU[5];

            rDOF[countMultiscaleDofs+1][1] += dDOFdEpsilonHomY[0]*bHomU[3] +
                                              dDOFdEpsilonHomY[1]*bHomU[4] +
                                              dDOFdEpsilonHomY[2]*bHomU[5];
#ifdef MYDEBUG
            {
                double interval(1e-5);
                double disp1[2],disp2[2],disptmp[2];
                mLogger << "dDofdU2 algo " << rDOF[countMultiscaleDofs][2] << "  " << rDOF[countMultiscaleDofs+1][2] << "\n";
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
                mLogger << "dDofdU2 cdf " << (disp2[0]-disp1[0])/interval << "  " << (disp2[1]-disp1[1])/interval << "\n";
                if (fabs((disp2[0]-disp1[0])/interval-rDOF[countMultiscaleDofs][2])>1e-2 || fabs((disp2[1]-disp1[1])/interval-rDOF[countMultiscaleDofs+1][2])>1e-2)
                {
                    mLogger << "[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDofdU2 is not correct" << "\n";
                    throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofsCalculatedDispdGlobalDofs] dDofdU2 is not correct");
                }
                const_cast<StructureMultiscale*>(&*this)->mCrackOpening[1] -= interval;
                const_cast<StructureMultiscale*>(&*this)->CalculateHomogeneousEngineeringStrain();
                mLogger << "\n";
                mLogger << "crack opening2 " << mCrackOpening[0] << " " << mCrackOpening[1]<< "\n";

            }
#endif
            //scale the derivative
            rDOF[countMultiscaleDofs][1]   *= mScalingFactorCrackOpening;
            rDOF[countMultiscaleDofs+1][1] *= mScalingFactorCrackOpening;

            //derivative of displacement with respect to homogeneous strain (is equal to derivative w.r.t. total strain)
            rDOF[countMultiscaleDofs][2]   = dDOFdEpsilonHomX[0]*mScalingFactorEpsilon;
            rDOF[countMultiscaleDofs][3]   = dDOFdEpsilonHomX[1]*mScalingFactorEpsilon;
            rDOF[countMultiscaleDofs][4]   = dDOFdEpsilonHomX[2]*mScalingFactorEpsilon;

            //std::cout << "rDOF[countMultiscaleDofs][3] " << rDOF[countMultiscaleDofs][3] << " " << rDOF[countMultiscaleDofs][4] << " " << rDOF[countMultiscaleDofs][5] << std::endl;

            rDOF[countMultiscaleDofs+1][2]   = dDOFdEpsilonHomY[0]*mScalingFactorEpsilon;
            rDOF[countMultiscaleDofs+1][3]   = dDOFdEpsilonHomY[1]*mScalingFactorEpsilon;
            rDOF[countMultiscaleDofs+1][4]   = dDOFdEpsilonHomY[2]*mScalingFactorEpsilon;

            countMultiscaleDofs+=2;

            if (countMultiscaleDofs>numMultiscaleDofs)
            {
                throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofs] countMultiscaleDofs>numMultiscaleDofs internal error.");
            }
        }
    }
#ifdef MYDEBUG
    mLogger << "End of routine" << "\n" << "\n";
#endif
    /*    std::cout << "rDOF" << std::endl;
        for (int count=0; count<rDOF.size(); count++)
        {
        	for (int count2=0; count2<6; count2++)
        		std::cout << rDOF[count][count2] << " " ;
        	std::cout << std::endl;
        }
    */
    if (countMultiscaleDofs!=numMultiscaleDofs)
        throw MechanicsException("[NuTo::StructureMultiscale::CalculatedDispdGlobalDofs] countMultiscaleDofs!=numMultiscaleDofs internal error.");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
NuTo::Error::eError NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const
{
    throw MechanicsException("[NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric] not implemented for StructureMultiscale");
}

//! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
//! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
//! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
//! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
NuTo::Error::eError NuTo::StructureMultiscale::BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const
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

//! @brief calculate the distance of a point to the crack
//! @param rCoordinates[2] coordinates of the point
//! @return distance to crack
double NuTo::StructureMultiscale::CalculateDistanceToCrack2D(double rCoordinates[2])const
{
    return sin(mCrackAngle)*(mCenterDamage[0]+mShiftCenterDamage[0]-rCoordinates[0]) - cos(mCrackAngle)*(mCenterDamage[1]+mShiftCenterDamage[1]-rCoordinates[1]);
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureMultiscale::CalculateHomogeneousEngineeringStrain()
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double factor=mlCoarseScaleCrack/mCoarseScaleArea;

    mEpsilonHom.mEngineeringStrain[0] = mEpsilonTot.mEngineeringStrain[0] - factor * sinAlpha *(-cosAlpha*mCrackOpening[0]+sinAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[1] = mEpsilonTot.mEngineeringStrain[1] - factor * cosAlpha *(sinAlpha*mCrackOpening[0]+cosAlpha*mCrackOpening[1]);
    mEpsilonHom.mEngineeringStrain[2] = mEpsilonTot.mEngineeringStrain[2] - factor * ((2.*cosAlpha*cosAlpha-1.)* mCrackOpening[0] - 2.*sinAlpha*cosAlpha * mCrackOpening[1]);
    /*    mEpsilonHom.mEngineeringStrain[0] = 0;
    mEpsilonHom.mEngineeringStrain[1] = 0;
    mEpsilonHom.mEngineeringStrain[2] = 0;

    mLogger << "calculation of e hom  " << "\n";

    //    mLogger << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << "\n";

    //#ifdef MYDEBUG
    mLogger << "mEpsilonTot " << mEpsilonTot.mEngineeringStrain[0] << " " << mEpsilonTot.mEngineeringStrain[1] << " " << mEpsilonTot.mEngineeringStrain[2] << "\n";
    mLogger << "mCrackOpening " << mCrackOpening[0] << " " << mCrackOpening[1] << " " << factor << "\n";
    mLogger << "mEpsilonHom " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2] << "\n";
    //#endif
    */
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureMultiscale::SetTotalEngineeringStrain(const EngineeringStrain2D& rTotalEngineeringStrain)
{
    mEpsilonTot = rTotalEngineeringStrain;
    CalculateHomogeneousEngineeringStrain();
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureMultiscale::SetTotalEngineeringStrainConstraint(const EngineeringStrain2D& rTotalEngineeringStrain)
{
    mEpsilonTotConstraint = rTotalEngineeringStrain;
}

//! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
void NuTo::StructureMultiscale::SetDeltaTotalEngineeringStrain(const EngineeringStrain2D& rDeltaTotalEngineeringStrain)
{
    mDeltaEpsilonTot = rDeltaTotalEngineeringStrain;
}

//! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
void NuTo::StructureMultiscale::SetPrevTotalEngineeringStrain(const EngineeringStrain2D& rPrevTotalEngineeringStrain)
{
    mPrevEpsilonTot = rPrevTotalEngineeringStrain;
}

//! @brief returns the total strain
const NuTo::EngineeringStrain2D& NuTo::StructureMultiscale::GetTotalEngineeringStrain()const
{
    return mEpsilonTot;
}

//! @brief returns the total strain
const NuTo::EngineeringStrain2D& NuTo::StructureMultiscale::GetTotalEngineeringStrainConstraint()const
{
    return mEpsilonTotConstraint;
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

//! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack orientation and crack opening
//! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
//! @paramter rbHomU [0-2] wrt ux [3-5] wrt uy
void NuTo::StructureMultiscale::GetdEpsilonHomdCrack(double rbHomU[6])const
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));

    // coarse scale model is not a square or not aligned with the principal axes
    double factor = -mlCoarseScaleCrack/mCoarseScaleArea;
    rbHomU[0] = factor*(-sinAlpha*cosAlpha);
    rbHomU[1] = -rbHomU[0];
    rbHomU[2] = factor*(2.*cosAlpha*cosAlpha-1.);

    rbHomU[3] = factor*(sinAlpha*sinAlpha);
    rbHomU[4] = factor*(cosAlpha*cosAlpha);
    rbHomU[5] = factor*(-2.*sinAlpha*cosAlpha);
}

//! @brief renumbers the global dofs in the structure after
void NuTo::StructureMultiscale::ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering)
{
    mDOFCrackOpening[0] = rMappingInitialToNewOrdering[mDOFCrackOpening[0]];
    mDOFCrackOpening[1] = rMappingInitialToNewOrdering[mDOFCrackOpening[1]];
    mDOFGlobalTotalStrain[0]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[0]];
    mDOFGlobalTotalStrain[1]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[1]];
    mDOFGlobalTotalStrain[2]  = rMappingInitialToNewOrdering[mDOFGlobalTotalStrain[2]];
    mDOFPeriodicBoundaryDisplacements[0] = rMappingInitialToNewOrdering[mDOFPeriodicBoundaryDisplacements[0]];
    mDOFPeriodicBoundaryDisplacements[1] = rMappingInitialToNewOrdering[mDOFPeriodicBoundaryDisplacements[1]];
    mDOFPeriodicBoundaryDisplacements[2] = rMappingInitialToNewOrdering[mDOFPeriodicBoundaryDisplacements[2]];
}

//! @brief calculates the crack angle for elastic solutions (initial value, no scaling with previous crack angle)
//! @return crack angle in the range 0..Pi
double  NuTo::StructureMultiscale::CalculateCrackAnglePrincipalStrain(const EngineeringStrain2D& rStrain)const
{
    FullMatrix<double> strain(2,2);

    const double *dataPtr;
    dataPtr = rStrain.GetData();

    strain(0,0) = dataPtr[0];
    strain(0,1) = dataPtr[2]*0.5;
    strain(1,0) = strain(0,1);
    strain(1,1) = dataPtr[1];

    //std::cout << "strain " << "\n";
    //strain.Info(12,5);

    NuTo::FullMatrix<double> eigenVectors;
    NuTo::FullMatrix<double> eigenValues;
    strain.EigenVectorsSymmetric(eigenVectors);
    strain.EigenValuesSymmetric(eigenValues);

    //align the eigenvector w.r.t. the previous crack angle
    /*    if (eigenVectors(0,0)*cos(mPrevCrackAngleElastic)+eigenVectors(1,0)*sin(mPrevCrackAngleElastic)<0)
        {
        	eigenVectors(0,0)*=-1.;
        	eigenVectors(1,0)*=-1.;
        }
    */
    //mLogger << "eigenvalues of epsilon tot" << "\n";
    //eigenValues.Trans().Info(12,10);
    //mLogger << "eigenvectors of epsilon tot" << "\n";
    //eigenVectors.Info(12,3);


    assert(eigenValues(0,0)<=eigenValues(1,0));
    double alpha = atan2(eigenVectors(1,0),eigenVectors(0,0));
    //std::cout << "alpha " << alpha*180./M_PI;
    /*    if (alpha-mPrevCrackAngleElastic<0)
        	while (fabs(alpha-mPrevCrackAngleElastic)>0.5*M_PI)
        	   alpha+=M_PI;
        else
        	while (fabs(alpha-mPrevCrackAngleElastic)>0.5*M_PI)
        	   alpha-=M_PI;
        //std::cout << " prevAngle " << mPrevCrackAngleElastic*180./M_PI << " alpha_mod " << alpha*180./M_PI << std::endl;
        //std::cout << "alpha " << alpha*180./M_PI << "\n";
    */
    return alpha;
}

//! @brief add a constraint equation for the crack opening (normal crack opening non negativ)
//! @parameter rPenaltyStiffness penalty stiffness for augmented Lagrangian
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLagrangeGlobalCrackOpeningNormal(double rPenaltyStiffness)
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
    mNodeNumberingRequired = true;
    return id;
}

//! @brief add a constraint equation for the total strain
//! @parameter rStrain applied strain (rhs)
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLinearGlobalTotalStrain()
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

    mConstraintMap.insert(id, new NuTo::ConstraintLinearGlobalTotalStrain(this));
    mConstraintTotalStrain = id;
    mNodeNumberingRequired = true;
    return id;
}

//! @brief add a linear constraint equation for the crack opening
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLinearGlobalCrackOpening(double rRHS, const NuTo::FullMatrix<double>& rDirection)
{
    int id = 0;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    ConstraintLinearGlobalCrackOpening *mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this, rDirection, rRHS);
    mNodeNumberingRequired = true;
    mConstraintMap.insert(id, mConst);
    return id;
}

//! @brief add a linear constraint equation for the crack opening
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLinearGlobalCrackOpeningTangential(double rRHS)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 0;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0)=1;
    direction(1,0)=0;
    ConstraintLinearGlobalCrackOpening *mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, rRHS);
    mConstraintMap.insert(id, mConst);
    mNodeNumberingRequired = true;
    mConstraintTangentialCrackOpening = id;
    return id;
}

//! @brief add a linear constraint equation for the crack opening
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLinearGlobalCrackOpeningNormal(double rRHS)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintNormalCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 0;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0)=0;
    direction(1,0)=1;
    ConstraintLinearGlobalCrackOpening *mConst = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, rRHS);
    mConstraintMap.insert(id, mConst);
    mNodeNumberingRequired = true;
    mConstraintNormalCrackOpening = id;
    return id;
}

//! @brief add a linear constraint equation for the additional shape functions depscribing the fluctuation boundary displacements
//! @return id of the constraint
int NuTo::StructureMultiscale::CreateConstraintLinearPeriodicBoundaryShapeFunctions(int rShapeFunction, double rRHS)
{
    if (rShapeFunction!=0 && rShapeFunction!=1 &&  rShapeFunction!=2)
        throw MechanicsException("[NuTo::StructureMultiscale::CreateConstraintLinearPeriodicBoundaryShapeFunctions] in 2D, there are only 3 additional boundary shape functions.");
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintPeriodicBoundaryShapeFunctions[rShapeFunction]);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }

    int id = 0;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    ConstraintLinearPeriodicBoundaryShapeFunctions *mConst = new NuTo::ConstraintLinearPeriodicBoundaryShapeFunctions(this, rShapeFunction, rRHS);
    mConstraintMap.insert(id, mConst);
    mConstraintPeriodicBoundaryShapeFunctions[rShapeFunction] = id;
    mNodeNumberingRequired = true;
    return id;
}

//! @brief delete the constraint for the tangential crack opening
void NuTo::StructureMultiscale::ConstraintDeleteTangentialCrackOpening()
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintTangentialCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }
    mConstraintTangentialCrackOpening = -1;
    mNodeNumberingRequired = true;
}

//! @brief delete the constraint for the normal crack opening
void NuTo::StructureMultiscale::ConstraintDeleteNormalCrackOpening()
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintNormalCrackOpening);
    if (it!=mConstraintMap.end())
    {
        mConstraintMap.erase(it);
    }
    mNodeNumberingRequired = true;
    mConstraintNormalCrackOpening = -1;
}

//! @brief set constraint for fine scale fluctuations on the boundary
void NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions(bool rHomogeneousDomain)
{
    this->mNodeNumberingRequired = true;

    if (rHomogeneousDomain)
    {
        boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintFineScalePeriodicHomogeneous);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
        mConstraintFineScalePeriodicHomogeneous = -1;

        //insert constraints for the homogeneous domain
        //find unused integer id
        int id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupBoundaryNodesHomogeneous);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions] Group(homogeneous) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureBase::CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions] Group(homogeneous) is not an node group.");
        const Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(this, nodeGroup));
        mConstraintFineScalePeriodicHomogeneous = id;

    }
    else
    {
        boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintFineScalePeriodicDamage);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
        mConstraintFineScalePeriodicDamage = -1;

        //insert constraints for the damage domain
        //find unused integer id
        int id(0);
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupBoundaryNodesDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureBase::CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions] Group(damage) is not an node group.");
        const Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(this, nodeGroup));
        mConstraintFineScalePeriodicDamage = id;
    }

    return ;

}

//! @brief add constraints for the finescale displacement of the boundary nodes
//! @parameter rDispX .. displacement in x-direction
//! @parameter rDispY .. displacement in y-direction
//! @parameter rHomogeneousDomain .. true for the homogeneous domain, false for the damage domain
void NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacements(double rDispX, double rDispY, bool rHomogeneousDomain)
{
    if (rHomogeneousDomain)
    {
        //insert constraints for the homogeneous domain
        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupBoundaryNodesHomogeneous);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacements] Group(homogeneous) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureMultiscale::CreateConstraintLinearFineScaleDisplacements] Group(homogeneous) is not an node group.");
        const Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

        //insert constraints for the homogeneous domain
        boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintFineScaleHomogeneousX);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }

        it = mConstraintMap.find(mConstraintFineScaleHomogeneousY);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }

        //insert constraints for the x direction
        //find unused integer id
        int id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }
        NuTo::FullMatrix<double> direction(2,1);
        direction(0,0) = 1;
        direction(1,0) = 0;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,rDispX));
        mConstraintFineScaleHomogeneousX = id;

        //insert constraints for the y direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }
        direction(0,0) = 0;
        direction(1,0) = 1;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,rDispY));
        mConstraintFineScaleHomogeneousY = id;
    }
    else
    {
        boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupBoundaryNodesDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureBase::CalculatePeriodicBoundaryShapeFunctions] Group(damage) is not an node group.");
        const Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

        boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(mConstraintFineScaleDamageX);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
        it = mConstraintMap.find(mConstraintFineScaleDamageY);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
        //insert constraints for the x direction
        //find unused integer id
        int id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        NuTo::FullMatrix<double> direction(2,1);
        direction(0,0) = 1;
        direction(1,0) = 0;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,rDispX));
        mConstraintFineScaleDamageX = id;

        //insert constraints for the y direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        direction(0,0) = 0;
        direction(1,0) = 1;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,rDispY));
        mConstraintFineScaleDamageY = id;
    }
}


//! @brief set periodic boundary conditions for a 2D structure
//! @parameter rGroupBoundaryNodes ... boundary nodes
//! @parameter rStrain ... strain
int NuTo::StructureMultiscale::ConstraintLinearSetFineScaleDisplacementsPeriodicNodeGroup(int rGroupBoundaryNodes, const EngineeringStrain2D& rStrain)
{
    if (mDimension!=2)
        throw MechanicsException("[NuTo::StructureMultiscale::ConstraintLagrangeSetDisplacementNodeGroup] only implemented for 2D structures.");
    this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupBoundaryNodes);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureMultiscale::ConstraintLinearSetFineScaleDisplacementsPeriodicNodeGroup] Group with the given identifier does not exist.");
    Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

    mConstraintMap.insert(id, new NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D(nodeGroup,rStrain));
    mNodeNumberingRequired = true;

    return id;
}

//! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
//! @param load factor
void NuTo::StructureMultiscale::SetLoadFactor(double rLoadFactor)
{
    //set the total strain and calculate from the existing crack opening the homogeneous strain
    //and the elastic crack angle
    EngineeringStrain2D curStrain(mPrevEpsilonTot+mDeltaEpsilonTot*rLoadFactor);
    SetTotalEngineeringStrainConstraint(curStrain);

    //mLogger << "set total strain " << curEngineeringStrain.mEngineeringStrain[0]<< " " << curEngineeringStrain.mEngineeringStrain[1]<<" " << curEngineeringStrain.mEngineeringStrain[2] << "\n";
}

//! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureMultiscale::PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const
{
    mLogger << "Convergence after " << rNumNewtonIterations << " Newton iterations, curLoadFactor " << rLoadFactor << ", deltaLoadFactor "<< rDeltaLoadFactor <<  "\n";
    mLogger << "   Crack opening " << mCrackOpening[0] << "[t] " << mCrackOpening[1] << "[n]"<<  "\n";
    mLogger << "   total strain " << mEpsilonTot.mEngineeringStrain[0] << " " << mEpsilonTot.mEngineeringStrain[1] << " " << mEpsilonTot.mEngineeringStrain[2]  <<  "\n";
    mLogger << "   hom strain " << mEpsilonHom.mEngineeringStrain[0] << " " << mEpsilonHom.mEngineeringStrain[1] << " " << mEpsilonHom.mEngineeringStrain[2]  <<  "\n";
    NuTo::FullMatrix<double> multiplier;
    if (mConstraintNormalCrackOpening!=-1)
    {
        boost::ptr_map<int,ConstraintBase>::const_iterator it = mConstraintMap.find(mConstraintNormalCrackOpening);
        if (it!=mConstraintMap.end())
        {
            if ((*it)->second->GetNumLagrangeMultipliers()>0)
            {
                ConstraintLagrangeGetMultiplier(mConstraintNormalCrackOpening,multiplier);
                mLogger << " lambda crack opening " << multiplier(0,0)  <<  "\n";    //! @brief calculates the average stress
            }
        }
    }


    /*    NuTo::FullMatrix<double> cohesiveForce(2,1);
        CalculateCohesiveForce(cohesiveForce);
        mLogger << " cohesive force " << cohesiveForce(0,0) << "(t) " << cohesiveForce(1,0) << "(n) " << "\n";

        NuTo::FullMatrix<double> engineeringStress(3,1);
        ElementTotalGetAverageStress(mFineScaleAreaDamage, engineeringStress);
        mLogger << " average stress vector " << engineeringStress(0,0) << " " << engineeringStress(1,0)<< " " <<  engineeringStress(2,0) << "\n";
        mLogger << " stress vector on the crack " << engineeringStress(0,0)*(-sin(mCrackAngle))+ engineeringStress(2,0)*cos(mCrackAngle)<< "(t) " ;
                                          mLogger << engineeringStress(2,0)*(-sin(mCrackAngle))+ engineeringStress(1,0)*cos(mCrackAngle)<< "(n) " << "\n";
    */
//    std::stringstream ssLoadStep;
//    ssLoadStep << rLoadStep;
//    if (rLoadFactor==1)
//        this->ExportVtkDataFile(mResultDirectoryAfterLineSearch+std::string("/../FinescaleConcurrentMultiscale") + ssLoadStep.str()+ std::string(".vtk"));
    /*    if (rNumNewtonIterations>5)
        {
    		mLogger << "please press return after update Newton iteration on fine scale of ip " << mIPName <<  "\n"<<  "\n";
    		std::string str;
    		getline (std::cin,str);
        }
    */
}

//! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureMultiscale::PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, const NuTo::FullMatrix<double>& rResidualVector)const
{
    mLogger << "After Linesearch, iteration " << rNewtonIteration <<
            ", line search factor " << rLineSearchFactor <<
            ", load factor " << rLoadFactor <<  "\n";
    mLogger << " crack angle " << mCrackAngle*180./M_PI << "[0..360]";

    if (mDOFCrackOpening[0]<rResidualVector.GetNumRows())
    {
        mLogger << ", crack opening " << mCrackOpening[0] << "[t](F="<< rResidualVector(mDOFCrackOpening[0],0) << "), ";
    }
    else
    {
        mLogger << ", crack opening " << mCrackOpening[0] << ", ";
    }
    if (mDOFCrackOpening[1]<rResidualVector.GetNumRows())
    {
        mLogger << mCrackOpening[1] << "[t](F=" << rResidualVector(mDOFCrackOpening[1],0)<< "), ";
    }
    else
    {
        mLogger <<  mCrackOpening[1] << "[t], ";
    }
    mLogger  << "\n";
    for (int count =0; count<3; count++)
    {
        if (mDOFPeriodicBoundaryDisplacements[count]<rResidualVector.GetNumRows())
        {
            mLogger << ", periodic boundary displacement " << count << ": "<< mPeriodicBoundaryDisplacements[count] << "(F="<< rResidualVector(mDOFPeriodicBoundaryDisplacements[count],0) << "\n";
        }
        else
        {
            mLogger << ", periodic boundary displacement " << count << ": "<< mPeriodicBoundaryDisplacements[count] << "\n";
        }

    }

    boost::ptr_map<int,ConstraintBase>::const_iterator it = mConstraintMap.find(mConstraintNormalCrackOpening);
    if (it!=mConstraintMap.end())
    {
        if ((*it)->second->GetNumLagrangeMultipliers()>0)
        {
            NuTo::FullMatrix<double> multiplier;
            (*it)->second->AsConstraintLagrange()->GetLagrangeMultiplier(multiplier);
            NuTo::FullMatrix<int> multiplierDOF;
            (*it)->second->AsConstraintLagrange()->GetDofsLagrangeMultiplier(multiplierDOF);
            double value = multiplier(0,0);
            int residual = rResidualVector(multiplierDOF(0,0),0);
            mLogger << ", Lagrange multiplier " << value << "(F=" << residual << ")";
        }
    }
    mLogger << "\n";
    int row, col;
    double maxError = rResidualVector.Max(row,col);
    mLogger << " max residual " << maxError << " at " << row << "\n";
    if (mGroupElementsDamage!=-1)
        mLogger << " scaling factor damage " << GetScalingFactorDamage() << " scaling factor homogeneous " << GetScalingFactorHomogeneous() << "\n";
    else
        mLogger << " scaling factor homogeneous " << GetScalingFactorHomogeneous() << "\n";
}

//! @brief do a postprocessing step in each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureMultiscale::PostProcessDataInLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, double rPrevResidual)const
{
    mLogger << "  in linesearch " << rNewtonIteration <<
            ", line search factor " << rLineSearchFactor <<
            ", residual " << rResidual <<
            ", previous residual " << rPrevResidual << "\n";
}

//! @brief initialize some stuff before a new load step (e.g. output directories for visualization, if required)
void NuTo::StructureMultiscale::InitBeforeNewLoadStep(int rLoadStep)
{
    mLogger <<" mDOF crack opening " << mDOFCrackOpening[0] << " " << mDOFCrackOpening[1] ;
    mLogger << " mDOF epsilontot " << mDOFGlobalTotalStrain[0] << " " << mDOFGlobalTotalStrain[1] << " " << mDOFGlobalTotalStrain[2] << "\n";
    mLogger << " mDOFPeriodicBoundaryDisplacements " << mDOFPeriodicBoundaryDisplacements[0] << " " << mDOFPeriodicBoundaryDisplacements[1] << " " << mDOFPeriodicBoundaryDisplacements[2];
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
    //mLogger << "energy damage " << energyDamage << " * " << GetScalingFactorDamage() << " = " << energyDamage*GetScalingFactorDamage() << "\n";
    //mLogger << "energy homogeneous " << energyHomogeneous << " * " << GetScalingFactorHomogeneous() << " = " << energyHomogeneous*GetScalingFactorHomogeneous() << "\n";
    return energyDamage*GetScalingFactorDamage()+energyHomogeneous*GetScalingFactorHomogeneous();
}

//! @brief is only true for structure used as multiscale (structure in a structure)
void NuTo::StructureMultiscale::ScaleCoordinates(double rCoordinates[3])const
{
    //rCoordinates[2] is unchanged
    if (mScaleWRTDamageCenter)
    {
        double ratio(0.6*sqrt(mCoarseScaleArea/mFineScaleArea));
        rCoordinates[0] = (rCoordinates[0]-mCenterDamage[0])*ratio + mCenterMacro[0];
        rCoordinates[1] = (rCoordinates[1]-mCenterDamage[1])*ratio + mCenterMacro[1];
    }
    else
    {
        double ratio(0.6*sqrt(mCoarseScaleArea/mFineScaleArea));
        rCoordinates[0] = (rCoordinates[0]-mCenterHomogeneous[0])*ratio + mCenterMacro[0];
        rCoordinates[1] = (rCoordinates[1]-mCenterHomogeneous[1])*ratio + mCenterMacro[1];
    }
}

//! @brief sets the result directory where the results are written to
void NuTo::StructureMultiscale::SetResultDirectory(std::string rResultDirectory)
{
    mResultDirectory = rResultDirectory;
    if (boost::filesystem::exists(mResultDirectory))    // does p actually exist?
    {
        if (!boost::filesystem::is_directory(mResultDirectory))      // is p a directory?
            throw MechanicsException("[NuTo::StructureMultiscale::SetResultDirectory] result directory is existing, but not a directory.");
    }
    else if (!boost::filesystem::create_directory(mResultDirectory))
        throw MechanicsException("[NuTo::StructureMultiscale::SetResultDirectory] result directory could not be created.");
    boost::filesystem::path resultDir(mResultDirectory);
    resultDir /= mIPName;
    if (boost::filesystem::exists(resultDir))    // does ip directory actually exist?
    {
        if (!boost::filesystem::is_directory(resultDir))      // is p a directory?
            throw MechanicsException("[NuTo::StructureMultiscale::SetResultDirectory] result directory for ip is existing, but not a directory.");
    }
    else
    {
        if (!boost::filesystem::create_directory(resultDir))
            throw MechanicsException("[NuTo::StructureMultiscale::SetResultDirectory] directory for ip could not be created in result directory.");
    }
}

void NuTo::StructureMultiscale::CalculateStiffness(NuTo::FullMatrix<double>& rStiffness, bool rPeriodic)
{
#ifdef ENABLE_MUMPS
    mLogger.OpenFile();

    //release the constraint for the total strain
    ConstraintBase* linearTotalStrainConstraintPtr = ConstraintRelease(mConstraintTotalStrain); // now this is not deleted

    //and the crack opening in tangential
    int id = 0;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    int constraintCrackOpeningTangential = id;
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0) = 1.;
    direction(1,0) = 0.;
    ConstraintLinearGlobalCrackOpening *mConst2 = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, 0);
    mConstraintMap.insert(id, mConst2);

    //and the crack opening in normal direction
    id = 0;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    int constraintCrackOpeningNormal = id;
    direction(0,0) = 0.;
    direction(1,0) = 1.;
    ConstraintLinearGlobalCrackOpening *mConst3 = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, 0);
    mConstraintMap.insert(id, mConst3);

    //and periodic boundary conditions
    if (rPeriodic)
    {
        EngineeringStrain2D strain;
        mNodeNumberingRequired = true;
        mConstraintFineScalePeriodicHomogeneous = ConstraintLinearSetFineScaleDisplacementsPeriodicNodeGroup(mGroupBoundaryNodesHomogeneous,strain);
    }
    else
    {
        //set fine scale disp of all boundary node groups (hom) to zero (x & y)
        bool homogeneousDomain(true);
        CreateConstraintLinearFineScaleDisplacements(0,0,homogeneousDomain);
    }

    //fix the additional dofs related to the periodic boundary displacements
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(0, 0);
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(1, 0);
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(2, 0);

    // use Schur complement to calculate the stiffness
    this->NodeBuildGlobalDofs();
    NuTo::FullMatrix<double> activeDOF, dependentDOF;
    this->NodeExtractDofValues(activeDOF,dependentDOF);
    this->NodeMergeActiveDofValues(activeDOF);
    this->ElementTotalUpdateTmpStaticData();
    SparseMatrixCSRVector2General<double>  matrixJJ(this->GetNumActiveDofs(), this->GetNumActiveDofs());
    FullMatrix<double> rhsVector(mNumDofs - mNumActiveDofs,1);
    this->BuildGlobalCoefficientMatrix0(matrixJJ, rhsVector);
    if (rhsVector.Abs().Max()>1e-6)
    {
        throw MechanicsException("[NuTo::Multiscale::CalculateStiffness] RHS vector not zero, previous equilibrium iteration was not successfull.");
    }
    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    matrixJJ.Max(maxValue);
    matrixJJ.Min(minValue);
    //mLogger << "min and max " << minValue << " , " << maxValue << "\n";

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    int numRemoved = matrixJJ.RemoveZeroEntries(ToleranceZeroStiffness,0);
    std::cout << "number of entries removed in stiffness matrix (int %) for calculation of initial macroscopic elastic stiffness " << 100*numRemoved/(matrixJJ.GetNumEntries()+numRemoved) << "% with tolerance " << ToleranceZeroStiffness << "\n";

    //calculate schur complement
    NuTo::FullMatrix<int> schurIndicesMatrix(3,1);
    NuTo::FullMatrix<double> stiffness(3,3);
    //attention, the index is in the zero based indexing system
    schurIndicesMatrix(0,0) = mDOFGlobalTotalStrain[0];
    schurIndicesMatrix(1,0) = mDOFGlobalTotalStrain[1];
    schurIndicesMatrix(2,0) = mDOFGlobalTotalStrain[2];

    NuTo::SparseDirectSolverMUMPS mumps;
    NuTo::SparseMatrixCSRGeneral<double> stiffnessFineScale(matrixJJ);

    //std::cout << "stiffness " << matrixJJ.GetNumRows() << " " << matrixJJ.GetNumColumns() << std::endl;
    //NuTo::FullMatrix<double> stiffnessFull(matrixJJ);
    //stiffnessFull.Info(12,3);

    //std::cout << "schurIndicesMatrix " <<  std::endl;
    //schurIndicesMatrix.Trans().Info(12,3);

    stiffnessFineScale.SetOneBasedIndexing();
    mumps.SchurComplement(stiffnessFineScale,schurIndicesMatrix,rStiffness);
    //scale with the dimension of the structure (area)
    rStiffness*=1./(mCoarseScaleArea*mScalingFactorEpsilon*mScalingFactorEpsilon);

    std::cout << "Schur stiffness algo" << std::endl;
    rStiffness.Info(12,3);

    //delete the constraints again
    this->ConstraintDelete(constraintCrackOpeningTangential);
    this->ConstraintDelete(constraintCrackOpeningNormal);
    if (rPeriodic)
    {
        //this->ConstraintDelete(mConstraintFineScalePeriodicDamage);
        //mConstraintFineScalePeriodicDamage = -1;
        this->ConstraintDelete(mConstraintFineScalePeriodicHomogeneous);
        mConstraintFineScalePeriodicHomogeneous = -1;
    }
    else
    {
        //this->ConstraintDelete(mConstraintFineScaleDamageX);
        //mConstraintFineScaleDamageX = -1;
        //this->ConstraintDelete(mConstraintFineScaleDamageY);
        //mConstraintFineScaleDamageY = -1;
        this->ConstraintDelete(mConstraintFineScaleHomogeneousX);
        mConstraintFineScaleHomogeneousX = -1;
        this->ConstraintDelete(mConstraintFineScaleHomogeneousY);
        mConstraintFineScaleHomogeneousY = -1;
    }

    for (int count=0; count<3; count++)
    {
        this->ConstraintDelete(mConstraintPeriodicBoundaryShapeFunctions[count]);
        this->mConstraintPeriodicBoundaryShapeFunctions[count] = -1;
    }

    //reinsert the constraint for the total strain
    this->ConstraintAdd(mConstraintTotalStrain,linearTotalStrainConstraintPtr);
    mLogger.CloseFile();
#else // ENABLE_MUMPS
	throw MechanicsException("[NuTo::StructureMultiscale::CalculateStiffness] mumps is not enabled.");
#endif //ENABLE_MUMPS
}

void NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions(double rDeltaStrain)
{
    mLogger.OpenFile();

    //and the crack opening in tangential
    int id = 0;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    int constraintCrackOpeningTangential = id;
    NuTo::FullMatrix<double> direction(2,1);
    direction(0,0) = 1.;
    direction(1,0) = 0.;
    ConstraintLinearGlobalCrackOpening *mConst2 = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, 0);
    mConstraintMap.insert(id, mConst2);

    //and the crack opening in normal direction
    id = 0;
    it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }
    int constraintCrackOpeningNormal = id;
    direction(0,0) = 0.;
    direction(1,0) = 1.;
    ConstraintLinearGlobalCrackOpening *mConst3 = new NuTo::ConstraintLinearGlobalCrackOpening(this, direction, 0);
    mConstraintMap.insert(id, mConst3);

    //and periodic boundary conditions
    EngineeringStrain2D strain;
    mConstraintFineScalePeriodicHomogeneous = ConstraintLinearSetFineScaleDisplacementsPeriodicNodeGroup(mGroupBoundaryNodesHomogeneous,strain);

    /*    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(mGroupBoundaryNodesDamage);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(damage) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        	throw MechanicsException("[NuTo::StructureBase::CalculatePeriodicBoundaryShapeFunctions] Group(damage) is not an node group.");
        const Group<NodeBase> *nodeGroup = itGroup->second->AsGroupNode();

        it = mConstraintMap.find(mConstraintFineScaleDamageX);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
    	it = mConstraintMap.find(mConstraintFineScaleDamageY);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }
        //insert constraints for the x direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        direction(0,0) = 1;
        direction(1,0) = 0;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,0));
        mConstraintFineScaleDamageX = id;

        //insert constraints for the y direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }

        direction(0,0) = 0;
        direction(1,0) = 1;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,0));
        mConstraintFineScaleDamageX = id;

        //insert constraints for the homogeneous domain
        itGroup = mGroupMap.find(mGroupBoundaryNodesHomogeneous);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(homogeneous) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        	throw MechanicsException("[NuTo::StructureBase::CalculatePeriodicBoundaryShapeFunctions] Group(homogeneous) is not an node group.");
        nodeGroup = itGroup->second->AsGroupNode();

        //insert constraints for the homogeneous domain
        it = mConstraintMap.find(mConstraintFineScaleHomogeneousX);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }

        it = mConstraintMap.find(mConstraintFineScaleHomogeneousY);
        if (it!=mConstraintMap.end())
        {
            mConstraintMap.erase(it);
        }

        //insert constraints for the x direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }
        direction(0,0) = 1;
        direction(1,0) = 0;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,0));
        mConstraintFineScaleHomogeneousX = id;

        //insert constraints for the y direction
        //find unused integer id
        id = 0;
        it = mConstraintMap.find(id);
        while (it!=mConstraintMap.end())
        {
            id++;
            it = mConstraintMap.find(id);
        }
        direction(0,0) = 0;
        direction(1,0) = 1;
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(nodeGroup,direction,0));
        mConstraintFineScaleHomogeneousY = id;
    */

    //fix the additional dofs related to the periodic boundary displacements
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(0, 0);
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(1, 0);
    CreateConstraintLinearPeriodicBoundaryShapeFunctions(2, 0);

    //ConstraintInfo(10);

    for (int countShapeFunction=0; countShapeFunction<3; countShapeFunction++)
    {
        EngineeringStrain2D prevStrain;
        this->SetPrevTotalEngineeringStrain(prevStrain);

        SetLoadFactor(0);
        NodeBuildGlobalDofs();
        NuTo::FullMatrix<double> activeDOF, dependentDOF;
        NodeExtractDofValues(activeDOF,dependentDOF);
        NodeMergeActiveDofValues(activeDOF);

        // calculate engineering strain
        EngineeringStrain2D engineeringStrain;
        engineeringStrain.mEngineeringStrain[countShapeFunction] = rDeltaStrain;
        EngineeringStrain2D deltaStrain(engineeringStrain-prevStrain);
        SetDeltaTotalEngineeringStrain(deltaStrain);

        std::stringstream saveStream;
        bool hasBeenSaved(false);
        GetLogger() << "\n" << "****************************************************" << "\n";
        GetLogger() << " Calculate Stress for the fine scale solution to get the periodic boundary shape functions "  << "\n";
        GetLogger() << " engineering strain " <<  engineeringStrain.mEngineeringStrain[0] << " "
                    <<  engineeringStrain.mEngineeringStrain[1] << " "
                    <<  engineeringStrain.mEngineeringStrain[2]  << "\n";
        try
        {
            NewtonRaphson(true, saveStream, hasBeenSaved);
        }
        catch(MechanicsException& e)
        {
            e.AddMessage(std::string("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Error in performing Newton-iteration on fine scale for ip ") + GetIPName());
            throw e;
        }
        catch(...)
        {
            throw MechanicsException(std::string("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Error in performing Newton-iteration on fine scale for ip ") + GetIPName());
        }
        std::stringstream ss;
        ss << countShapeFunction;
        std::cout << "total number of elements " << mElementMap.size();
        std::cout << ", active Dofs " << this->mNumActiveDofs << " dependent dofs " << this->mNumDofs-this->mNumActiveDofs << "\n";
#ifdef ENABLE_VISUALIZE
        ExportVtkDataFile(mResultDirectory + "/periodicBC" + ss.str() +".vtk");
#endif //ENABLE_VISUALIZE
        double max(0);
        // loop over all nodes in the damaged region
        /*		boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(mGroupBoundaryNodesDamage);
        	    if (itGroup==mGroupMap.end())
        	        throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(damage) with the given identifier does not exist.");
        	    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        	    	throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(damage) is not an node group.");
        	    Group<NodeBase> *nodeGroupDamage = itGroup->second->AsGroupNode();
        	    assert(nodeGroupDamage!=0);
        	    for (Group<NodeBase>::iterator itNode=nodeGroupDamage->begin(); itNode!=nodeGroupDamage->end();itNode++)
        	    {
        	    	itNode->second->SetShapeFunctionMultiscalePeriodic(countShapeFunction);
        	    	double X = itNode->second->GetShapeFunctionMultiscalePeriodicX()[countShapeFunction];
        	    	if (fabs(X)>max)
        	    		max=fabs(X);
        	    	double Y = itNode->second->GetShapeFunctionMultiscalePeriodicY()[countShapeFunction];
        	    	if (fabs(Y)>max)
        	    		max=fabs(Y);
        	    }
        	    */
        // loop over all nodes in the homogeneous region
        boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(mGroupBoundaryNodesHomogeneous);
        if (itGroup==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(homogeneous) with the given identifier does not exist.");
        if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] Group(homogeneous) is not an node group.");
        Group<NodeBase> *nodeGroupHomogeneous = itGroup->second->AsGroupNode();
        assert(nodeGroupHomogeneous!=0);
        for (Group<NodeBase>::iterator itNode=nodeGroupHomogeneous->begin(); itNode!=nodeGroupHomogeneous->end(); itNode++)
        {
            itNode->second->SetShapeFunctionMultiscalePeriodic(countShapeFunction);
            double X = itNode->second->GetShapeFunctionMultiscalePeriodicX()[countShapeFunction];
            if (fabs(X)>max)
                max=fabs(X);
            double Y = itNode->second->GetShapeFunctionMultiscalePeriodicY()[countShapeFunction];
            if (fabs(Y)>max)
                max=fabs(Y);
        }
        double scaleFactor(1./max);
        /*for (Group<NodeBase>::iterator itNode=nodeGroupDamage->begin(); itNode!=nodeGroupDamage->end();itNode++)
        {
        	itNode->second->ScaleShapeFunctionMultiscalePeriodic(countShapeFunction,scaleFactor);
        }
        */
        for (Group<NodeBase>::iterator itNode=nodeGroupHomogeneous->begin(); itNode!=nodeGroupHomogeneous->end(); itNode++)
        {
            itNode->second->ScaleShapeFunctionMultiscalePeriodic(countShapeFunction,scaleFactor);
        }

        //restore structure
        if (hasBeenSaved)
        {
            throw MechanicsException("[NuTo::StructureMultiscale::CalculatePeriodicBoundaryShapeFunctions] there should not be any nonlinearity here, check your input.");
        }
        else
        {
            //set load factor to zero in order to get the same ordering of the displacements as before the routine
            SetLoadFactor(0);
            NodeBuildGlobalDofs();
            NodeMergeActiveDofValues(activeDOF);
        }
    }

    //delete the constraints again
    this->ConstraintDelete(constraintCrackOpeningTangential);
    this->ConstraintDelete(constraintCrackOpeningNormal);
    this->ConstraintDelete(mConstraintFineScalePeriodicHomogeneous);
    mConstraintFineScalePeriodicHomogeneous = -1;
    for (int count=0; count<3; count++)
    {
        this->ConstraintDelete(mConstraintPeriodicBoundaryShapeFunctions[count]);
        this->mConstraintPeriodicBoundaryShapeFunctions[count] = -1;
    }

    //reinsert the constraint for the total strain
    mLogger.CloseFile();
}


//! @brief sets the result file for the converged solution where the results are written to
void NuTo::StructureMultiscale::SetResultLoadStepMacro(std::string rLoadStepMacro)
{
    mLoadStepMacro = rLoadStepMacro;
    mLogger.OpenFile(GetFileForLog());
    //close the file again, since there is only a limited number of open file handlers allowed
    mLogger.CloseFile();
}

#ifdef ENABLE_VISUALIZE
void NuTo::StructureMultiscale::VisualizeCrack(VisualizeUnstructuredGrid& rVisualize)const
{
    double GlobalPointCoor1[3];
    double ratio(0.3*sqrt(mCoarseScaleArea/mFineScaleArea));
    GlobalPointCoor1[0] = mCenterMacro[0]+(mShiftCenterDamage[0]+1.0*mlFineScaleCrack*cos(mCrackAngle))*ratio;
    GlobalPointCoor1[1] = mCenterMacro[1]+(mShiftCenterDamage[1]+1.0*mlFineScaleCrack*sin(mCrackAngle))*ratio;
    GlobalPointCoor1[2] = 0.;

    double GlobalPointCoor2[3];
    GlobalPointCoor2[0] = mCenterMacro[0]+(mShiftCenterDamage[0]-1.0*mlFineScaleCrack*cos(mCrackAngle))*ratio;
    GlobalPointCoor2[1] = mCenterMacro[1]+(mShiftCenterDamage[1]-1.0*mlFineScaleCrack*sin(mCrackAngle))*ratio;
    GlobalPointCoor2[2] = 0.;
    unsigned int Points[2];
    Points[0] = rVisualize.AddPoint(GlobalPointCoor1);
    Points[1] = rVisualize.AddPoint(GlobalPointCoor2);
    rVisualize.AddLineCell(Points);
}
#endif //ENABLE_VISUALIZE

/*
//! @brief performs a Newton Raphson iteration (displacement and/or load control)
//! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
//!             be careful, store it only once
NuTo::StructureEnum::eErrorCode NuTo::StructureMultiscale::NewtonRaphson(bool rSaveStructureBeforeUpdate,
        std::stringstream& rSaveStringStream,
        bool& rIsSaved, bool rInitialStateInEquilibrium)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
try
{

    // start analysis
    double deltaLoadFactor(mMaxDeltaLoadFactor);
    double curLoadFactor;

	this->SetLoadFactor(0);
	this->NodeBuildGlobalDofs();
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    this->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
    //std::cout << "last converged dofs " << std::endl;
    //displacementsActiveDOFsLastConverged.Trans().Info(12,5);

    //init some auxiliary variables
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> intForceVector;
    NuTo::FullMatrix<double> extForceVector;
    NuTo::FullMatrix<double> rhsVector;
    NuTo::FullMatrix<double> intForceVectorInit;

    //calculate the initial out of balance force
    if (rInitialStateInEquilibrium==false)
    {
		try
		{
			this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
		    //std::cout << "norm displacementsDependentDOFsLastConverged " << displacementsDependentDOFsLastConverged.Norm() << "\n";
    	    //std::cout << "displacements (norm) for zero load is " << displacementsActiveDOFsLastConverged.Norm() << "\n";
			//displacementsActiveDOFsLastConverged.Trans().Info(8,3);
			this->ElementTotalUpdateTmpStaticData();
    	    this->BuildGlobalGradientInternalPotentialVector(intForceVectorInit);
    	    //std::cout << "out of balance force (norm) for zero load is " << intForceVectorInit.Norm() << "\n";
    	    //intForceVectorInit.Trans().Info(8,3);
		}
		catch(MechanicsException& e)
		{
 			mLogger << "[NuTo::StructureBase::NewtonRaphson] no convergence for constitutive model in unloaded state. " << "\n";
			e.AddMessage("[NuTo::StructureBase::NewtonRaphson] [NuTo::StructureBase::NewtonRaphson] no convergence for constitutive model in unloaded state.");
			throw e;
		}
		catch(...)
		{
			throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] no convergence for constitutive model in unloaded state.");
		}
    }

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    //NuTo::SparseDirectSolverMKLPardiso mySolver;
    #ifdef SHOW_TIME
    if (mShowTime==true)
        mySolver.SetShowTime(true);
    else
        mySolver.SetShowTime(false);
    #endif
    bool convergenceInitialLoadStep(false);
    int loadStep(1);
    InitBeforeNewLoadStep(loadStep);
    while (convergenceInitialLoadStep==false)
    {
		try
		{
			//calculate stiffness
			curLoadFactor=deltaLoadFactor;
			this->SetLoadFactor(curLoadFactor);
			this->NodeBuildGlobalDofs();
			if (mNumActiveDofs==0)
			{
				this->SetLoadFactor(1);
				this->NodeBuildGlobalDofs();
				NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
				NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
				this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
				this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
				this->ElementTotalUpdateTmpStaticData();
				return;
			}
			this->ElementTotalUpdateTmpStaticData();

			//mLogger<<" calculate stiffness 3891" << "\n";
			this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
		//    NuTo::FullMatrix<double>(stiffnessMatrixCSRVector2).Info(12,3);
		//    mLogger << "disp force vector "<< "\n";
		//    dispForceVector.Trans().Info(12,10);
		//Check the stiffness matrix
		//CheckStiffness();
		//mLogger << "stiffness is calculated in Newton Raphson " << "\n";
			//mLogger << "total energy of system " << ElementTotalGetTotalEnergy() << "\n";

			//update displacements of all nodes according to the new conre mat
			{
				NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
				NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
				this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
				this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
				this->ElementTotalUpdateTmpStaticData();
			}

			// build global external load vector and RHS vector
			this->BuildGlobalExternalLoadVector(extForceVector);
			//mLogger<<" calculate gradient 3913" << "\n";
			this->BuildGlobalGradientInternalPotentialVector(intForceVector);
			convergenceInitialLoadStep = true;
		}
		catch(MechanicsException& e)
		{
            if (e.GetError()==MechanicsException::NOCONVERGENCE)
            {
				//decrease load step
				deltaLoadFactor*=mDecreaseFactor;

				//restore initial state
				this->SetLoadFactor(0);
				this->NodeBuildGlobalDofs();
				this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);

				//check for minimum delta (this mostly indicates an error in the software
				if (deltaLoadFactor<mMinDeltaLoadFactor)
				{
					mLogger << "No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum " << "\n";
					e.AddMessage("[NuTo::StructureMultiscale::NewtonRaphson] No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum.");
					throw e;
				}
            }
            else
            {
            	e.AddMessage("[NuTo::StructureMultiscale::NewtonRaphson] Error in Newton-Raphson iteration.");
            	throw e;
            }
		}
    }
    rhsVector = extForceVector - intForceVector;
    //attention this is only different for the first iteration step
    //since the internal force due to the applied constraints is not considered for the first iteration
    //in order to balance it (no localization in the boundary region)
    //for the linesearch this internal force has to be considered in order to obtain for a linesearch
    //factor of zero the normRHS
    double normRHS = rhsVector.Norm();
    //rhsVector.Trans().Info(12,10);
    rhsVector = extForceVector + dispForceVector;
    if (rInitialStateInEquilibrium==false)
    {
        rhsVector -= intForceVectorInit;
    }
    //rhsVector.Trans().Info(12,10);
    //std::cout << "norm rhsVector " << rhsVector.Norm() << "\n";

//    {
//        double energy;
//        energy = this->ElementTotalGetTotalEnergy();
//        energy += this->ConstraintTotalGetTotalEnergy();
//        mLogger << "energy " << energy << "\n";
//    }

    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    if (stiffnessMatrixCSRVector2.GetNumColumns()==0)
    {
        maxValue = 1.;
        minValue = 1.;
    }
    else
    {
        stiffnessMatrixCSRVector2.Max(maxValue);
        stiffnessMatrixCSRVector2.Min(minValue);
    }
    //mLogger << "min and max " << minValue << " , " << maxValue << "\n";

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

    //mySolver.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale") + std::string("0") + std::string(".vtk"));

    //store the structure only once in order to be able to restore the situation before entering the routine
    rIsSaved = false;

    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    while (!convergenceStatusLoadSteps)
    {
        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
        double alpha(1.);
        NuTo::FullMatrix<double> displacementsActiveDOFs;
        NuTo::FullMatrix<double> displacementsDependentDOFs;
        int convergenceStatus(0);
        //0 - not converged, continue Newton iteration
        //1 - converged
        //2 - stop iteration, decrease load step
        while(convergenceStatus==0)
        {
            numNewtonIterations++;

            if (numNewtonIterations>mMaxNumNewtonIterations && alpha<0.25)
            {
                if (mVerboseLevel>5)
                {
                    mLogger << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << mMaxNumNewtonIterations << ")" << "\n";
                }
                convergenceStatus = 2; //decrease load step
                break;
            }

            // solve
            NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
            this->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);
            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            try
            {
                mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
            }
            catch(...)
            {
            	mLogger << "Error solving system of equations using mumps." << "\n";
            	if (mNumActiveDofs<1000)
                {
            	    NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                    if (mNumActiveDofs<30)
                    {
                        mLogger << "stiffness full" << "\n";
                        mLogger.Out(stiffnessMatrixFull,12,3);
                    }
                    NuTo::FullMatrix<double> eigenValues;
                    stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                    mLogger << "eigenvalues" << "\n";
                    mLogger.Out(eigenValues.Trans(),12,3);
                    NuTo::FullMatrix<double> eigenVectors;
                    stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                    mLogger << "eigenvector 1" << "\n";
                    mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
                }
                throw MechanicsException("[NuTo::StructureMultiscale::NewtonRaphson] Error solving system of equations using mumps.");
            }

            mLogger << " rhs norm " << rhsVector.Norm() << "delta disp norm " << deltaDisplacementsActiveDOFs.Norm() << "\n";
            //perform a linesearch
            alpha = 1.;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;

                //mLogger << " displacementsActiveDOFs" << "\n";
                //displacementsActiveDOFs.Trans().Info(10,3);
                this->NodeMergeActiveDofValues(displacementsActiveDOFs);
                this->ElementTotalUpdateTmpStaticData();

                // calculate residual
                try
                {
        			//mLogger<<" calculate gradient 4086" << "\n";
                    this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                    //mLogger << "intForceVector "  << "\n";
                    //mLogger.Out(intForceVector.Trans(),10,3,false);
                    rhsVector = extForceVector - intForceVector;
                    normResidual = rhsVector.Norm();
                    maxResidual = rhsVector.Abs().Max();

                    //mLogger << "total energy of system " << ElementTotalGetTotalEnergy() << "\n";
                    PostProcessDataInLineSearch(loadStep, numNewtonIterations, alpha, curLoadFactor, normResidual, normRHS);
                }
                catch(MechanicsException& e)
                {
                	if (e.GetError()==MechanicsException::NOCONVERGENCE)
                	{
						convergenceStatus=2;
						mLogger << "Constitutive model is not converging, try with smaller load step" << "\n";
                	}
                	else
                	{
                		e.AddMessage("[NuTo::StructureMultiscale::NewtonRaphson] Error calling the gradient (resforce) routine.");
                		throw e;
                	}
                }

                alpha*=0.5;
            }
            while(convergenceStatus!=2 && alpha>mMinLineSearchFactor && normResidual>normRHS*(1-0.5*alpha) && normResidual>mToleranceResidualForce && maxResidual>mToleranceResidualForce);

            if (convergenceStatus==2)
            {
            	break;
            }

            this->PostProcessDataAfterLineSearch(loadStep, numNewtonIterations, 2.*alpha, curLoadFactor, normResidual, rhsVector);
    		//std::string str;
    		//getline (std::cin,str);

            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>mToleranceResidualForce && maxResidual>mToleranceResidualForce)
            {
                convergenceStatus=2;
                {
                    if (mNumActiveDofs<1000)
                    {
                    	mLogger << "System is not converging." << "\n";
                        NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                        if (mNumActiveDofs<30)
                        {
                            mLogger << "stiffness full" << "\n";
                            mLogger.Out(stiffnessMatrixFull,12,3);
                        }
                        NuTo::FullMatrix<double> eigenValues;
                        stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                        mLogger << "eigenvalues" << "\n";
                        mLogger.Out(eigenValues.Trans(),12,10);
                        NuTo::FullMatrix<double> eigenVectors;
                        stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                        mLogger << "eigenvector 1" << "\n";
                        mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
                        //mLogger << "alpha " << (eigenVectors(mDOFCrackAngle,0)) << "\n";
                    }

                }
                break;
            }

            //mLogger << "\n" << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<"\n";

            //check convergence
            if (normResidual<mToleranceResidualForce || maxResidual<mToleranceResidualForce)
            {
                this->PostProcessDataAfterConvergence(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                convergenceStatus=1;
                //CheckStiffness();
                //NodeInfo(12);
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            normRHS = rhsVector.Norm();
            //build new stiffness matrix
			//mLogger<<" calculate stiffness 4166" << "\n";
            this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            //mLogger << dispForceVector.Norm() << "\n";
//check stiffness
//CheckStiffness();
            //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";
        }

        if (convergenceStatus==1)
        {
            // the update is only required to allow for a stepwise solution procedure in the fine scale model
            // a final update is only required for an update on the macroscale, otherwise,the original state has
            // to be reconstructed.

            if (curLoadFactor>1-1e-8)
            {
            	if (rSaveStructureBeforeUpdate==false)
            	{
					this->ElementTotalUpdateStaticData();
	                this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
            	}
				convergenceStatusLoadSteps=true;
            }
            else
            {
            	this->ElementTotalUpdateStaticData();
                this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                //store the last converged step in order to be able to go back to that state
                displacementsActiveDOFsLastConverged  = displacementsActiveDOFs;

                //eventually increase load step
                if (mAutomaticLoadstepControl)
                {
                    if (numNewtonIterations<mMinNumNewtonIterations)
                    {
                        deltaLoadFactor*=mIncreaseFactor;
                    }
                    if (deltaLoadFactor>mMaxDeltaLoadFactor)
                        deltaLoadFactor = mMaxDeltaLoadFactor;
                }

                //increase displacement
                curLoadFactor+=deltaLoadFactor;
                if (curLoadFactor>1)
                {
                    deltaLoadFactor -= curLoadFactor -1.;
                    curLoadFactor=1;
                }
            }
            loadStep++;
        }
        else
        {

        	assert(convergenceStatus==2);
            if (mAutomaticLoadstepControl==false)
                throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence with the prescribed number of Newton iterations.",MechanicsException::NOCONVERGENCE);

            mLogger << "no convergence with current step size (" << deltaLoadFactor << "), current not converging load factor " << curLoadFactor << "\n";
            if (mNumActiveDofs<1000)
            {
            	InitBeforeNewLoadStep(0);
            	NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                NuTo::FullMatrix<double> eigenValues;
                stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                mLogger << "eigenvalues" << "\n";
                mLogger.Out(eigenValues.Trans(),12,10);
                NuTo::FullMatrix<double> eigenVectors;
                stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                int count = 0;
                while (eigenValues(count,0)<0.1)
                {
					mLogger << "eigenvector " << count+1 << "\n";
					mLogger.Out(eigenVectors.GetColumn(count).Trans(),12,3);
					//mLogger << "alpha " << (eigenVectors(mDOFCrackAngle,0)) << "\n";
					count++;
                }

                //CheckStiffness();
            }
            mLogger << "no convergence with current step size (" << deltaLoadFactor << "), current not converging load factor " << curLoadFactor << "\n";
            mLogger << "check stiffness " << "\n";
            //CheckStiffness();
            mLogger << "and continue with smaller load step " << "\n";

            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
            curLoadFactor-=deltaLoadFactor;

            //set the previous displacement state
            this->SetLoadFactor(curLoadFactor);

            // build global dof numbering
            this->NodeBuildGlobalDofs();

            //set previous converged displacements
            this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
            this->ElementTotalUpdateTmpStaticData();
            //this first part of the routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
            //for an update on the coarse scale
            //as a consequence, in an iterative solution with updates in between, the initial state has to be restored after leaving the routine
            if (rSaveStructureBeforeUpdate==true && rIsSaved==false)
            {
                assert(curLoadFactor==0);
                //store the structure only once in order to be able to restore the situation before entering the routine
            	this->SaveStructure(rSaveStringStream);
            	rIsSaved = true;
            }

            //decrease load step
            deltaLoadFactor*=mDecreaseFactor;
            curLoadFactor+=deltaLoadFactor;

            //check for minimum delta (this mostly indicates an error in the software
            if (deltaLoadFactor<mMinDeltaLoadFactor)
            {
            	mLogger << "return with a MechanicsNoConvergenceException " << "\n";
            	mLogger << "coordinates of the macropoint " << mCenterMacro[0] << " "<< mCenterMacro[1] << "\n";
                throw NuTo::MechanicsException("[NuTo::StructureMultiscale::NewtonRaphson]: No convergence, delta strain factor smaller than minimum.",MechanicsException::NOCONVERGENCE);
            }
        }

        if (!convergenceStatusLoadSteps)
        {
            //update new displacement of RHS

            // build global dof numbering
            bool convergenceConstitutive(false);
            while (convergenceConstitutive==false)
            {
				try
				{
		            this->SetLoadFactor(curLoadFactor);
					this->NodeBuildGlobalDofs();

					//update stiffness in order to calculate new dispForceVector, but still with previous displacement state
					//mLogger<<" calculate stiffness 4301" << "\n";
					this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
					//int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
					//int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
					//mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

					//update displacements of all nodes according to the new conre mat
					NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
					NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
					this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
					this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
					this->ElementTotalUpdateTmpStaticData();

					// calculate initial residual for next load step
					//mLogger<<" calculate gradient 4315" << "\n";
					this->BuildGlobalGradientInternalPotentialVector(intForceVector);
					convergenceConstitutive = true;

				}
				catch(MechanicsException& e)
				{
		            if (e.GetError()==MechanicsException::NOCONVERGENCE)
		            {
						if (mAutomaticLoadstepControl==false)
							throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence with the prescibed number of Newton iterations.",MechanicsException::NOCONVERGENCE);
						mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

						//calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
						//this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
						//otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
						curLoadFactor-=deltaLoadFactor;

						//set the previous displacement state
						this->SetLoadFactor(curLoadFactor);

						// build global dof numbering
						this->NodeBuildGlobalDofs();

						//set previous converged displacements
						this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
						this->ElementTotalUpdateTmpStaticData();

						//decrease load step
						deltaLoadFactor*=mDecreaseFactor;
						curLoadFactor+=deltaLoadFactor;

						//check for minimum delta (this mostly indicates an error in the software
						if (deltaLoadFactor<mMinDeltaLoadFactor)
						{
							mLogger << "return with a MechanicsNoConvergenceException " << "\n";
							throw NuTo::MechanicsException("[NuTo::StructureMultiscale::NewtonRaphson]: No convergence, delta strain factor smaller than minimum.",MechanicsException::NOCONVERGENCE);
						}
		            }
		            else
		            {
		            	e.AddMessage("[NuTo::StructureMultiscale::NewtonRaphson]: error calculating new displacement increment.");
		            	throw e;
		            }
				}
            }

            //update rhs vector for next Newton iteration
            rhsVector = extForceVector - intForceVector;
            normRHS = rhsVector.Norm();
            //attention this is only different for the first iteration step (load application)
            //since the internal force due to the applied constraints is not considered for the first iteration
            //in order to balance it (no localization in the boundary region)
            //for the linesearch this internal force has to be considered in order to obtain for a linesearch
            //factor of zero the normRHS
            rhsVector = dispForceVector  + extForceVector;
            if (rInitialStateInEquilibrium==false && curLoadFactor==deltaLoadFactor)
            {
                rhsVector -= intForceVectorInit;
            }
        }
    }
}
catch (MechanicsException& e)
{
    e.AddMessage("[NuTo::StructureMultiscale::NewtonRaphson] error performing Newton-Raphson iteration.");
    throw e;
}
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureMultiscale::NewtonRaphson] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureMultiscale::NewtonRaphson] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
}
*/

bool NuTo::StructureMultiscale::CheckStiffness()
{
    //be carefull this routine performs a node merge, which modifies the displacements
    //as a result, the stiffness calculated here might be different from the one if you just call the stiffness routine
    //this is especially true for the first step of the Newton iteration in displacement control situation
    //where the stiffness of the old state is calulated on purpose and the multiplied by the difference between prescribed dependent dofs and actual dependent dofs

    mLogger << "test of stiffness still included, node merge is called!!! " << "\n";
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;

    //recalculate stiffness
    this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    mLogger << "active dof values " << "\n";
    displacementsActiveDOFsCheck.Trans().Info(12,4);
    mLogger << "dependent dof values " << "\n";
    displacementsDependentDOFsCheck.Trans().Info(12,4);
    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();

    //check element stiffness
    //NuTo::FullMatrix<double> difference;
    //int maxErrorElement = ElementTotalCoefficientMatrix_0_Check(1e-8, difference);
    //std::cout << "maximum stiffness error in " << maxErrorElement;

    //calculate global stiffness
    this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    //this->ConstraintInfo(10);

    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
    //std::cout<<"stiffness matrix" << "\n";
    //stiffnessMatrixCSRVector2Full.Info(10,3);
    double interval(1e-11);
    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2_CDF(stiffnessMatrixCSRVector2.GetNumRows(), stiffnessMatrixCSRVector2.GetNumColumns());
    NuTo::FullMatrix<double> intForceVector1, intForceVector2, intForceVectorCDF(stiffnessMatrixCSRVector2.GetNumRows(),1);
    double energy1,energy2;
    this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();
    this->BuildGlobalGradientInternalPotentialVector(intForceVector1);
    //this->NodeInfo(10);
    energy1 = this->ElementTotalGetTotalEnergy();
    energy1 += this->ConstraintTotalGetTotalEnergy();
    //std::cout << "check stiffness:: energy1 "<<  energy1 << "\n";

    for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
        //for (int count=94; count<95; count++)
    {
        displacementsActiveDOFsCheck(count,0)+=interval;
        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        this->ElementTotalUpdateTmpStaticData();
        this->BuildGlobalGradientInternalPotentialVector(intForceVector2);
        //std::cout << "check stiffness:: intForceVector2"<< "\n";
        //intForceVector2.Trans().Info(10,6);
        //this->ConstraintInfo(10);
        energy2 = this->ElementTotalGetTotalEnergy();
        energy2 += this->ConstraintTotalGetTotalEnergy();
        stiffnessMatrixCSRVector2_CDF.SetColumn(count,(intForceVector2-intForceVector1)*(1./interval));
        intForceVectorCDF(count,0) = (energy2-energy1)/interval;
        displacementsActiveDOFsCheck(count,0)-=interval;
    }
    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    this->ElementTotalUpdateTmpStaticData();

    //if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max()>1e-1)
    if (true)
    {
        if (stiffnessMatrixCSRVector2Full.GetNumRows()<100)
        {
            mLogger << "globalStiffnessMatrix algo" << "\n";
            mLogger.Out(stiffnessMatrixCSRVector2Full,10,3,false);
            mLogger << "\n" << "globalStiffnessMatrix cdf" << "\n";
            mLogger.Out(stiffnessMatrixCSRVector2_CDF,10,3,false);
            mLogger<< "\n" << "error" << "\n";
            mLogger.Out((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full),10,3,false);
        }
        else
        {
            //extract the first 6x6 block
            NuTo::FullMatrix<double> blockAlgo(stiffnessMatrixCSRVector2Full.GetBlock(0,0,20,20));
            NuTo::FullMatrix<double> blockCDF(stiffnessMatrixCSRVector2_CDF.GetBlock(0,0,20,20));
            NuTo::FullMatrix<double> blockDelta(blockAlgo-blockCDF);

            mLogger << "block algo " << "\n";
            mLogger.Out(blockAlgo,10,3,false);

            mLogger << "block cdf " << "\n";
            mLogger.Out(blockCDF,10,3,false);

            mLogger << "block delta " << "\n";
            mLogger.Out(blockDelta,10,3,false);
        }

        double maxError;
        int row,col;
        maxError = (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max(row,col);
        mLogger << "maximum error stiffness is " << maxError << " at (" << row << "," << col << ") with abs value in correct matrix " << stiffnessMatrixCSRVector2Full(row,col) << "\n";

        if (stiffnessMatrixCSRVector2Full.GetNumRows()<100)
        {
            mLogger<< "\n" << "intForceVector algo" << "\n";
            mLogger.Out(intForceVector1.Trans(),10,3,false);
            mLogger<< "\n" << "intForceVector cdf" << "\n";
            mLogger.Out(intForceVectorCDF.Trans(),10,3,false);
            mLogger << "\n" << "error" << "\n";
            mLogger.Out((intForceVector1-intForceVectorCDF).Abs().Trans(),10,3);
        }
        else
        {
            //extract the first 6x6 block
            NuTo::FullMatrix<double> vecAlgo(intForceVector1.GetBlock(0,0,20,1));
            NuTo::FullMatrix<double> vecCDF(intForceVectorCDF.GetBlock(0,0,20,1));
            NuTo::FullMatrix<double> vecDelta(vecAlgo-vecCDF);

            mLogger << "res algo " << "\n";
            mLogger.Out(vecAlgo.Trans(),10,3,false);

            mLogger << "res cdf " << "\n";
            mLogger.Out(vecCDF.Trans(),10,3,false);

            mLogger << "res delta " << "\n";
            mLogger.Out(vecDelta.Trans(),10,3,false);
        }

        maxError = (intForceVector1-intForceVectorCDF).Abs().Max(row,col);
        mLogger << "maximum error resforce is " << maxError << " at (" << row << "," << col << ") with abs value in correct matrix " << intForceVector1(row,col) << "\n";

        //throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max()>1e-1)
        {
            //NodeInfo(10);
            mLogger << "stiffness is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<< "\n";
        }
        return false;
    }
    else
    {
        mLogger << "stiffness is OK "<< "\n";
        return true;
    }
}

//! @brief calculates the length of a crack in the fine scale
//! be careful, that the length, square/round, the angle and the offset of the finescale have to be set before
void NuTo::StructureMultiscale::CalculateAndSetCrackLengthFineScale()
{
    if (mSquareFineScaleModel)
    {
        //calculate direction vector
        double p1[2],p2[2];
        double u[2];
        double v[2];
        double w[2];
        v[0] = cos(mCrackAngle);
        v[1] = sin(mCrackAngle);
        double intersectionPoints[2][2];
        int numIntersectionPoints(0);
        for (int count=0; count<4; count++)
        {
            switch(count)
            {
            case 0:
                p1[0] = -mlFineScale*0.5;
                p1[1] = -mlFineScale*0.5;
                p2[0] =  mlFineScale*0.5;
                p2[1] = -mlFineScale*0.5;
                break;
            case 1:
                p1[0] =  mlFineScale*0.5;
                p1[1] = -mlFineScale*0.5;
                p2[0] =  mlFineScale*0.5;
                p2[1] =  mlFineScale*0.5;
                break;
            case 2:
                p1[0] =  mlFineScale*0.5;
                p1[1] =  mlFineScale*0.5;
                p2[0] = -mlFineScale*0.5;
                p2[1] =  mlFineScale*0.5;
                break;
            case 3:
                p1[0] = -mlFineScale*0.5;
                p1[1] =  mlFineScale*0.5;
                p2[0] = -mlFineScale*0.5;
                p2[1] = -mlFineScale*0.5;
                break;
            default:
                throw MechanicsException("[NuTo::StructureMultiscale::CalculateAndSetCrackLengthFineScale] this should never happen.");
            }
            w[0] = p1[0]-mShiftCenterDamage[0];
            w[1] = p1[1]-mShiftCenterDamage[1];
            u[0] = p2[0]-p1[0];
            u[1] = p2[1]-p1[1];

            double det = v[0]*u[1]-v[1]*u[0];
            if (fabs(det)>1e-10)
            {
                double s = (v[1]*w[0]-v[0]*w[1])/det;
                if (s>=0 && s<1)
                {
                    intersectionPoints[numIntersectionPoints][0] = p1[0] +s*u[0];
                    intersectionPoints[numIntersectionPoints][1] = p1[1] +s*u[1];
                    numIntersectionPoints++;
                }
            }
        }
        if (numIntersectionPoints!=2)
        {
            throw MechanicsException("[NuTo::StructureMultiscale::CalculateAndSetCrackLengthFineScale] Did not found two intersection points, there is something wrong.");
        }
        std::cout << "crack angle " << mCrackAngle*180./M_PI << ", mlFineScale " << mlFineScale << ", offset " << mShiftCenterDamage[0] << " " << mShiftCenterDamage[1] << "\n";
        std::cout << "intersection point 1 "<< intersectionPoints[0][0] << " " << intersectionPoints[0][1] << "\n";
        std::cout << "intersection point 2 "<< intersectionPoints[1][0] << " " << intersectionPoints[1][1] << "\n";

        //calculate length of crack and return
        mlFineScaleCrack = sqrt((intersectionPoints[0][0]-intersectionPoints[1][0])*(intersectionPoints[0][0]-intersectionPoints[1][0])+
                                (intersectionPoints[0][1]-intersectionPoints[1][1])*(intersectionPoints[0][1]-intersectionPoints[1][1]));
        std::cout << "crack length " << mlFineScaleCrack << "\n";
    }
    else
        throw MechanicsException("[NuTo::StructureMultiscale::CalculateAndSetCrackLengthFineScale] not implemented for round fine scale model.");
}


void NuTo::StructureMultiscale::CreateDamageDomainFromHomogeneousDomain()
{

    std::map<NodeBase*, NodeBase* > old2NewNodePtr;
    std::map<ElementBase*, ElementBase* > old2NewElementPtr;
    NuTo::FullMatrix<double> offset(2,1);
    offset(0,0) = mlFineScale*1.2;
    mCenterDamage[0]+=mCenterHomogeneous[0]+offset(0,0);
    mCenterDamage[1]=mCenterHomogeneous[1];
    CopyAndTranslate(offset,old2NewNodePtr,old2NewElementPtr);

    //create the groups with boundary nodes
    {
        mGroupBoundaryNodesDamage = GroupCreate("Nodes");
        boost::ptr_map<int,GroupBase>::iterator itGroupNew = mGroupMap.find(mGroupBoundaryNodesDamage);

        boost::ptr_map<int,GroupBase>::const_iterator itGroupOld = mGroupMap.find(mGroupBoundaryNodesHomogeneous);
        if (itGroupOld==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateDamageDomainFromHomogeneousDomain] Group boundary nodes homogeneous with the given identifier does not exist.");
        if (itGroupOld->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain] Group boundary nodes homogeneous is not a node group.");
        const Group<NodeBase> *nodeGroupOld = itGroupOld->second->AsGroupNode();
        Group<NodeBase> *nodeGroupNew = itGroupNew->second->AsGroupNode();
        for (Group<NodeBase>::const_iterator itNodeOld=nodeGroupOld->begin(); itNodeOld!=nodeGroupOld->end(); itNodeOld++)
        {
            std::map<NodeBase*, NodeBase* >::iterator itNodeMap = old2NewNodePtr.find(itNodeOld->second);
            if (itNodeMap==old2NewNodePtr.end())
                throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain]  old boundary node is not in the map, this is an implementation error.");
            NodeBase* newNodePtr = (*itNodeMap).second;
            int newNodeId = NodeGetId(newNodePtr);
            nodeGroupNew->AddMember(newNodeId,newNodePtr);
            //set the node to be in the cracked domain instead of the homogeneous domain
            switch(newNodePtr->GetNodeType())
            {
            case Node::NodeDisplacementsMultiscale2D:
            case Node::NodeCoordinatesDisplacementsMultiscale2D:
                newNodePtr->AsNodeDisplacementsMultiscale2D()->SetCrackedDomain(true);
                break;
            default:
                throw MechanicsException("[NuTo::StructureMultiscale::CreateDamageDomainFromHomogeneousDomain] the boundary nodes should be multiscale nodes.");
            }
        }
    }

    //create the groups with all nodes
    {
        mGroupNodesDamage = GroupCreate("Nodes");
        boost::ptr_map<int,GroupBase>::iterator itGroupNew = mGroupMap.find(mGroupNodesDamage);

        boost::ptr_map<int,GroupBase>::const_iterator itGroupOld = mGroupMap.find(mGroupNodesHomogeneous);
        if (itGroupOld==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateDamageDomainFromHomogeneousDomain] Group nodes homogeneous with the given identifier does not exist.");
        if (itGroupOld->second->GetType()!=NuTo::Groups::Nodes)
            throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain] Group nodes homogeneous is not a node group.");
        const Group<NodeBase> *nodeGroupOld = itGroupOld->second->AsGroupNode();
        Group<NodeBase> *nodeGroupNew = itGroupNew->second->AsGroupNode();
        for (Group<NodeBase>::const_iterator itNodeOld=nodeGroupOld->begin(); itNodeOld!=nodeGroupOld->end(); itNodeOld++)
        {
            std::map<NodeBase*, NodeBase* >::iterator itNodeMap = old2NewNodePtr.find(itNodeOld->second);
            if (itNodeMap==old2NewNodePtr.end())
                throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain] old node is not in the map, this is an implementation error.");
            NodeBase* newNodePtr = (*itNodeMap).second;
            int newNodeId = NodeGetId(newNodePtr);
            nodeGroupNew->AddMember(newNodeId,newNodePtr);
        }
    }

    //create the groups with all elements
    {
        mGroupElementsDamage = GroupCreate("Elements");
        boost::ptr_map<int,GroupBase>::iterator itGroupNew = mGroupMap.find(mGroupElementsDamage);

        boost::ptr_map<int,GroupBase>::const_iterator itGroupOld = mGroupMap.find(mGroupElementsHomogeneous);
        if (itGroupOld==mGroupMap.end())
            throw MechanicsException("[NuTo::StructureMultiscale::CreateDamageDomainFromHomogeneousDomain] Group elements homogeneous with the given identifier does not exist.");
        if (itGroupOld->second->GetType()!=NuTo::Groups::Elements)
            throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain] Group elements homogeneous is not an element group.");
        const Group<ElementBase> *elementGroupOld = itGroupOld->second->AsGroupElement();
        Group<ElementBase> *elementGroupNew = itGroupNew->second->AsGroupElement();
        for (Group<ElementBase>::const_iterator itElementOld=elementGroupOld->begin(); itElementOld!=elementGroupOld->end(); itElementOld++)
        {
            std::map<ElementBase*, ElementBase* >::iterator itElementMap = old2NewElementPtr.find(itElementOld->second);
            if (itElementMap==old2NewElementPtr.end())
                throw MechanicsException("[NuTo::StructureBase::CreateDamageDomainFromHomogeneousDomain] old element is not in the map, this is an implementation error.");
            ElementBase* newElementPtr = (*itElementMap).second;
            int newElementId = ElementGetId(newElementPtr);
            elementGroupNew->AddMember(newElementId,newElementPtr);
        }
    }
}





#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::StructureMultiscale)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
