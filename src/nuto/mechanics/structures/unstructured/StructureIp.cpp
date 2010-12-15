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
#define MAXNUMNEWTONITERATIONS 20
#define MIN_DELTA_STRAIN_FACTOR 1e-7
#define MYDEBUG false
#define tolCrackOpeningConstraint 1e-5
//! @brief constructor
//! @param mDimension number of nodes
NuTo::StructureIp::StructureIp ( int rDimension )  : Structure ( rDimension )
{
    if (rDimension!=2)
        throw MechanicsException("[NuTo::StructureIp::StructureIp] The concurrent multiscale model is only implemented for 2D.");
    mCrackAngle = M_PI*(-0.2);
    mDOFCrackAngle = -1;
    mCrackOpening[0] = 0.0;
    mCrackOpening[1] = 0.0;
    mDOFCrackOpening[0] = -1;
    mDOFCrackOpening[1] = -1;
    mEpsilonHom.mEngineeringStrain[0] = 0.;
    mEpsilonHom.mEngineeringStrain[1] = 0.;
    mEpsilonHom.mEngineeringStrain[2] = 0.;
    mConstraintFineScaleX = -1;
    mConstraintFineScaleY = -1;
    mConstraintAlpha = -1;
    mlX = 0.;
    mlY = 0.;
    mCrackTransitionZone = 0.;
    mBoundaryNodesTransformed = false;
    mIPName = std::string("fineScaleIp");
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
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Structure)
       & BOOST_SERIALIZATION_NVP(mCrackAngle)
       & BOOST_SERIALIZATION_NVP(mDOFCrackAngle)
       & BOOST_SERIALIZATION_NVP(mCrackOpening)
       & BOOST_SERIALIZATION_NVP(mDOFCrackOpening)
       & BOOST_SERIALIZATION_NVP(mEpsilonHom)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleX)
       & BOOST_SERIALIZATION_NVP(mConstraintFineScaleY)
       & BOOST_SERIALIZATION_NVP(mlX)
       & BOOST_SERIALIZATION_NVP(mlY)
       & BOOST_SERIALIZATION_NVP(mCrackTransitionZone)
       & BOOST_SERIALIZATION_NVP(mBoundaryNodesTransformed)
       & BOOST_SERIALIZATION_NVP(mIPName);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structureIp" << std::endl;
#endif
    std::cout << "save/restore mCrackAngle " << mCrackAngle << std::endl;
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


void NuTo::StructureIp::Save(std::stringstream& previousState)const
{
#ifdef ENABLE_SERIALIZATION
    boost::archive::binary_oarchive oba(previousState);
    oba << (*this);
#else
    throw MechanicsException("[NuTo::StructureIp::Save] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

void NuTo::StructureIp::Restore(std::stringstream& previousState)
{
#ifdef ENABLE_SERIALIZATION
    boost::archive::binary_iarchive iba(previousState);
    iba >> (*this);
#else
    throw MechanicsException("[NuTo::StructureIp::Restore] Serialization is required - switch it on in the CMakeFile.txt");
#endif //ENABLE_SERIALIZATION
}

//! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
//! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
void NuTo::StructureIp::TransformBoundaryNodes(int rBoundaryNodesId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rBoundaryNodesId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()==Groups::Nodes)
        TransformBoundaryNodes(dynamic_cast<Group<NodeBase>*>(itGroup->second));
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureIp::TransformBoundaryNodes] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

void NuTo::StructureIp::TransformBoundaryNodes(Group<NodeBase>* rBoundaryNodes)
{
    //copy the boundary nodes first, since the iterators are invalidated in NodeExchangePtr
    std::vector<NodeBase*> nodeVec(rBoundaryNodes->size());
    unsigned int countNode=0;
    for (Group<NodeBase>::iterator itNode=rBoundaryNodes->begin(); itNode!=rBoundaryNodes->end(); itNode++, countNode++)
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
    mConstraintFineScaleX = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(rBoundaryNodes, direction, 0);

    direction(0,0) = 0.;
    direction(1,0) = 1.;
    mConstraintFineScaleY = this->ConstraintLinearSetFineScaleDisplacementNodeGroup(rBoundaryNodes, direction, 0);
    mBoundaryNodesTransformed = true;

    if (minX!=0 || minY!=0)
    {
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] left lower corner has to be at the origin.");
    }
    if (minX==maxX || minY==maxY)
        throw MechanicsException("[NuTo::StructureIp::TransformBoundaryNodes] structure has zero width or height, either check your boundary group or the full structure.");
    mlX = maxX;
    mlY = maxY;
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
    // with respect to alpha, ux, uy
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,3> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);
    double bAlphaRow, bURow[2];

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
                globalRowDof -= this->mNumActiveDofs;
                assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                rDependentDofGradientVector(globalRowDof,0) += elementVector(rowCount,0);
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
            }
        }

        elementIter++;
    }
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
    // with respect to ehomxx, ehomyy, gammahomxy, ux, uy, alpha
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,3> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);
    double bAlphaRow, bURow[2], bAlphaCol, bUCol[2];

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

        // write element contribution to global matrix
        for (unsigned int rowCount = 0; rowCount < elementMatrixGlobalDofsRow.size(); rowCount++)
        {
            int globalRowDof = elementMatrixGlobalDofsRow[rowCount];
            assert(elementVectorGlobalDofs[rowCount]==globalRowDof);
            if (globalRowDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
                {
                    if (fabs(elementMatrix(rowCount, colCount))>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
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

                //add the influence of the second order derivatives
                if (mDOFCrackAngle < this->mNumActiveDofs)
                {
                    rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    else
                    {
                        rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        //rMatrixKJ.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                    else
                    {
                        rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        //rMatrixKJ.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                }
                else
                {
                    //rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        //rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    else
                    {
                        //rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        //rMatrixKK.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        //rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                    else
                    {
                        //rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        //rMatrixKK.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
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

                    if (globalRowDof < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackOpening[0], elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackOpening[0] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[0]);
                        }

                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackOpening[1], elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackOpening[1] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                    }
                    if (mappingDofMultiscaleNode[globalRowDof]!=-1)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackAngle,      mDOFCrackAngle,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }

                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[0], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1]-this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                            }
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[1], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                    }
                }
                if (mappingDofMultiscaleNode[globalRowDof]!=-1)
                {
                    if (globalColumnDof < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackAngle, globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackOpening[0], globalColumnDof, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackOpening[1], globalColumnDof, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                    }
                    else
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackAngle, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackOpening[0], globalColumnDof - this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackOpening[1], globalColumnDof - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
        elementIter++;
    }
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
    // with respect to ehomxx, ehomyy, gammahomxy, ux, uy, alpha
    std::vector<int> mappingDofMultiscaleNode;
    std::vector<std::array<double,3> > dDOF;
    std::vector<std::array<double,3> > dDOF2;
    CalculatedDispdGlobalDofs(mappingDofMultiscaleNode,dDOF,dDOF2);
    double bAlphaRow, bURow[2], bAlphaCol, bUCol[2];

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
                        rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));

                    }
                    else
                    {
                        rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                    }
                }
            }
            else
            {
                for (unsigned int colCount = 0; colCount < elementMatrixGlobalDofsColumn.size(); colCount++)
                {
                    int globalColumnDof = elementMatrixGlobalDofsColumn[colCount];
                    if (globalColumnDof < this->mNumActiveDofs)
                    {
                         rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                    }
                    else
                    {
                         rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
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

                //add the influence of the second order derivatives
                if (mDOFCrackAngle < this->mNumActiveDofs)
                {
                    rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    else
                    {
                        rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixKJ.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                    else
                    {
                        rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixKJ.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                }
                else
                {
                    rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][0]);
                    if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                    {
                        rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0], elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    else
                    {
                        rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                        rMatrixKK.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][1]);
                    }
                    if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                    {
                        rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1], elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                    }
                    else
                    {
                        rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
                        rMatrixKK.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementVector(rowCount, 0) * dDOF2[theDofMapRow][2]);
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

                    if (globalRowDof < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackOpening[0], elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackOpening[0] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[0]);
                        }

                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, mDOFCrackOpening[1], elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, mDOFCrackOpening[1] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                    }
                    else
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }
                        else
                        {
                            rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bAlphaCol);
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[0], elementMatrix(rowCount, colCount) * bUCol[0]);
                        }
                        else
                        {
                            rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[0]);
                        }

                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[1], elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                    }
                    if (mappingDofMultiscaleNode[globalRowDof]!=-1)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackAngle,      mDOFCrackAngle,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[0] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }

                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackAngle, mDOFCrackOpening[1] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }
                        else
                        {
                            rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs,      mDOFCrackAngle - this->mNumActiveDofs,      bAlphaRow * elementMatrix(rowCount, colCount) * bAlphaCol);
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }

                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1], bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }

                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[0], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[0], mDOFCrackOpening[1]-this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }
                        else
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackAngle -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            rMatrixKK.AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[0] -this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[1], bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackOpening[0] -this->mNumActiveDofs, mDOFCrackOpening[1]-this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                            }
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                            }
                            if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixJK.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[0] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            rMatrixJJ.AddEntry(mDOFCrackOpening[1], mDOFCrackOpening[1], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }
                        else
                        {
                            if (mDOFCrackAngle < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackAngle - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bAlphaCol);

                            }
                            if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                            {
                                rMatrixKJ.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackOpening[0], bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            else
                            {
                                rMatrixKK.AddEntry(mDOFCrackOpening[1]- this->mNumActiveDofs, mDOFCrackOpening[0] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[0]);
                            }
                            rMatrixKK.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, mDOFCrackOpening[1] - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount) * bUCol[1]);
                        }

                    }
                }
                if (mappingDofMultiscaleNode[globalRowDof]!=-1)
                {
                    if (globalColumnDof < this->mNumActiveDofs)
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackAngle,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKJ.AddEntry(mDOFCrackAngle - this->mNumActiveDofs,      globalColumnDof, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackOpening[0], globalColumnDof, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKJ.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, globalColumnDof, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(mDOFCrackOpening[1], globalColumnDof, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKJ.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, globalColumnDof, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                    }
                    else
                    {
                        if (mDOFCrackAngle < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackAngle, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKK.AddEntry(mDOFCrackAngle - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bAlphaRow * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[0] < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackOpening[0], globalColumnDof - this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKK.AddEntry(mDOFCrackOpening[0] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bURow[0]  * elementMatrix(rowCount, colCount));
                        }
                        if (mDOFCrackOpening[1] < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(mDOFCrackOpening[1], globalColumnDof - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKK.AddEntry(mDOFCrackOpening[1] - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, bURow[1]  * elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
        elementIter++;
    }
}

//! @brief calculate the derivative of the displacements at the nodes with respect to crack opening and crack orientation
//! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
//! @param rDOF return value, for each dof, the corresponding derivatives (alpha, ux, uy, )
//! @param r2DOF return value, for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy, )
void NuTo::StructureIp::CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<std::array<double,3> >& rDOF, std::vector<std::array<double,3> >& rDOF2)const
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
    if (MYDEBUG)
    {
        double interval(1e-5);
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
        std::cout << std::endl;
    }
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
            if (MYDEBUG)
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

            //derivative of displacement with respect to discontinuity (crack orientation)
            GetdDisplacementdCrackOrientation(coord,&(rDOF[countMultiscaleDofs][0]),&(rDOF[countMultiscaleDofs+1][0]));

            rDOF[countMultiscaleDofs][0]  +=  dDOFdEpsilonHomX[0]*bHomAlpha[0] +
                                              dDOFdEpsilonHomX[1]*bHomAlpha[1] +
                                              dDOFdEpsilonHomX[2]*bHomAlpha[2];
            rDOF[countMultiscaleDofs+1][0] += dDOFdEpsilonHomY[0]*bHomAlpha[0] +
                                              dDOFdEpsilonHomY[1]*bHomAlpha[1] +
                                              dDOFdEpsilonHomY[2]*bHomAlpha[2];
            if (MYDEBUG)
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


            //derivative of displacement with respect to discontinuity (crack opening)
            GetdDisplacementdCrackOpening(coord,&(rDOF[countMultiscaleDofs][1]),&(rDOF[countMultiscaleDofs+1][1]));
            rDOF[countMultiscaleDofs][1]   += dDOFdEpsilonHomX[0]*bHomU[0] +
                                              dDOFdEpsilonHomX[1]*bHomU[1] +
                                              dDOFdEpsilonHomX[2]*bHomU[2];
            rDOF[countMultiscaleDofs+1][1] += dDOFdEpsilonHomY[0]*bHomU[0] +
                                              dDOFdEpsilonHomY[1]*bHomU[1] +
                                              dDOFdEpsilonHomY[2]*bHomU[2];
            if (MYDEBUG)
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
            rDOF[countMultiscaleDofs][2]   += dDOFdEpsilonHomX[0]*bHomU[3] +
                                              dDOFdEpsilonHomX[1]*bHomU[4] +
                                              dDOFdEpsilonHomX[2]*bHomU[5];
            rDOF[countMultiscaleDofs+1][2] += dDOFdEpsilonHomY[0]*bHomU[3] +
                                              dDOFdEpsilonHomY[1]*bHomU[4] +
                                              dDOFdEpsilonHomY[2]*bHomU[5];
            if (MYDEBUG)
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
    if (MYDEBUG)
        std::cout << "End of routine" << std::endl << std::endl;

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
    rDisplacements[0] = mEpsilonHom.mEngineeringStrain[0] * rCoordinates[0] + 0.5 * mEpsilonHom.mEngineeringStrain[2] * rCoordinates[1];
    rDisplacements[1] = mEpsilonHom.mEngineeringStrain[1] * rCoordinates[1] + 0.5 * mEpsilonHom.mEngineeringStrain[2] * rCoordinates[0];
}

//! @brief ... calculate the displacement based on the crack opening
//! @param rCoordinates ... coordinates of the point
//! @param rCoarseDisplacements ... return array of displacements
void NuTo::StructureIp::GetDisplacementsCrack2D(double rCoordinates[2], double rDisplacements[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);

    if (d<-mCrackTransitionZone)
    {
        rDisplacements[0] = 0.;
        rDisplacements[1] = 0.;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rDisplacements[0] = mCrackOpening[0];
            rDisplacements[1] = mCrackOpening[1];
        }
        else
        {
            //smooth transition from cracking to none cracking
            double factor(0.5+0.5*sin(0.5*M_PI*d));
            rDisplacements[0] = factor*mCrackOpening[0];
            rDisplacements[1] = factor*mCrackOpening[1];

        }
    }
}

//! @brief derivative of displacement with respect to homogeneous strain
//! @param rdX_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
//! @param rdY_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
void NuTo::StructureIp::GetdDisplacementdEpsilonHom(double rCoordinates[2], double rdX_dEpsilonHom[3], double rdY_dEpsilonHom[3])const
{
    rdX_dEpsilonHom[0] = rCoordinates[0];
    rdX_dEpsilonHom[1] = 0.;
    rdX_dEpsilonHom[2] = 0.5 * rCoordinates[1];

    rdY_dEpsilonHom[0] = 0.;
    rdY_dEpsilonHom[1] = rCoordinates[1];
    rdY_dEpsilonHom[2] = 0.5 * rCoordinates[0];
}

//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
//! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to crack opening (ux, uy)
void NuTo::StructureIp::GetdDisplacementdCrackOpening(double rCoordinates[2], double rdX_dCrackOpening[2], double rdY_dCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);

    if (d<-mCrackTransitionZone)
    {
        rdX_dCrackOpening[0] = 0.;
        rdX_dCrackOpening[1] = 0.;
        rdY_dCrackOpening[0] = 0.;
        rdY_dCrackOpening[1] = 0.;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rdX_dCrackOpening[0] = 1.;
            rdX_dCrackOpening[1] = 0.;
            rdY_dCrackOpening[0] = 0.;
            rdY_dCrackOpening[1] = 1.;
        }
        else
        {
            //smooth transition from cracking to none cracking
            double factor(0.5+0.5*sin(0.5*M_PI*d));

            rdX_dCrackOpening[0] = factor;
            rdX_dCrackOpening[1] = 0.;
            rdY_dCrackOpening[0] = 0.;
            rdY_dCrackOpening[1] = factor;

        }
    }
}

//! @brief second derivative of displacement with respect to alpha and discontinuity (crack opening)
//! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to alpha and crack opening (ux, uy)
//! @param rdY_dCrackOpening[2] return value, derivative of y-displacement with respect to alpha and crack opening (ux, uy)
void NuTo::StructureIp::Getd2Displacementd2CrackOpening(double rCoordinates[2], double rdX_dAlphaCrackOpening[2], double rdY_dAlphaCrackOpening[2])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);

    if (d<-mCrackTransitionZone)
    {
        rdX_dAlphaCrackOpening[0] = 0.;
        rdX_dAlphaCrackOpening[1] = 0.;
        rdY_dAlphaCrackOpening[0] = 0.;
        rdY_dAlphaCrackOpening[1] = 0.;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rdX_dAlphaCrackOpening[0] = 0.;
            rdX_dAlphaCrackOpening[1] = 0.;
            rdY_dAlphaCrackOpening[0] = 0.;
            rdY_dAlphaCrackOpening[1] = 0.;
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double dFactor(0.5*cos(0.5*M_PI*d)*dDdAlpha);

            rdX_dAlphaCrackOpening[0] = dFactor;
            rdX_dAlphaCrackOpening[1] = 0.;
            rdY_dAlphaCrackOpening[0] = 0.;
            rdY_dAlphaCrackOpening[1] = dFactor;
        }
    }
}


//! @brief derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, derivative of x-displacement with respect to crack orientation (alpha)
void NuTo::StructureIp::GetdDisplacementdCrackOrientation(double rCoordinates[2], double rdX_dAlpha[1], double rdY_dAlpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);

    if (d<-mCrackTransitionZone)
    {
        rdX_dAlpha[0] = 0.;
        rdY_dAlpha[0] = 0.;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rdX_dAlpha[0] = 0.;
            rdY_dAlpha[0] = 0.;
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double dFactor(0.25*M_PI*cos(0.5*M_PI*d)*dDdAlpha);

            rdX_dAlpha[0] = dFactor*mCrackOpening[0];
            rdY_dAlpha[0] = dFactor*mCrackOpening[1];
        }
    }
}
//! @brief second derivative of displacement with respect to discontinuity (crack opening)
//! @param rdX_dAlpha[2] return value, second derivative of x-displacement with respect to crack orientation (alpha)
//! @param rdy_dAlpha[2] return value, second derivative of y-displacement with respect to crack orientation (alpha)
void NuTo::StructureIp::Getd2Displacementd2CrackOrientation(double rCoordinates[2], double rd2X_d2Alpha[1], double rd2Y_d2Alpha[1])const
{
    //calculate distance to crack
    double d = CalculateDistanceToCrack2D(rCoordinates);

    if (d<-mCrackTransitionZone)
    {
        rd2X_d2Alpha[0] = 0.;
        rd2Y_d2Alpha[0] = 0.;
    }
    else
    {
        if (d>mCrackTransitionZone)
        {
            rd2X_d2Alpha[0] = 0.;
            rd2Y_d2Alpha[0] = 0.;
        }
        else
        {
            //smooth transition from cracking to none cracking
            double dDdAlpha = CalculatedDistanceToCrack2DdAlpha(rCoordinates);
            double d2DdAlpha = Calculated2DistanceToCrack2Dd2Alpha(rCoordinates);
            double dFactor(0.25*M_PI*(0.5*M_PI*sin(0.5*M_PI*d)*dDdAlpha+cos(0.5*M_PI*d)*d2DdAlpha));

            rd2X_d2Alpha[0] = dFactor*mCrackOpening[0];
            rd2Y_d2Alpha[0] = dFactor*mCrackOpening[1];
        }
    }
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
    double factor=gamma/(2.*mlX);

    mEpsilonHom.mEngineeringStrain[0] = mEpsilonTot.mEngineeringStrain[0] + factor * sinAlpha * mCrackOpening[0];
    mEpsilonHom.mEngineeringStrain[1] = mEpsilonTot.mEngineeringStrain[1] - factor * cosAlpha * mCrackOpening[1];
    mEpsilonHom.mEngineeringStrain[2] = mEpsilonTot.mEngineeringStrain[2] + factor * ( sinAlpha * mCrackOpening[1] - cosAlpha * mCrackOpening[0]);
}

//calculate from the existing crack opening and orientation the cracking strain
void NuTo::StructureIp::SetTotalEngineeringStrain(EngineeringStrain2D& rTotalEngineeringStrain)
{
    mEpsilonTot = rTotalEngineeringStrain;
    CalculateHomogeneousEngineeringStrain();
}

//! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack orientation and crack opening
//! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
//! @parameter rbHomAlpha dHom wrt alpha
//! @paramter rbHomU[0-2] for wrt ux [3-5] for wrt uy
//! @parameter bHessian depsilondalpha2[0-2], depsilondalphadux[3-5], depsilondalphadux[6-8]
void NuTo::StructureIp::GetdEpsilonHomdCrack(double rbHomAlpha[3], double rbHomU[6], double rbHessian[9])const
{
    double sinAlpha(sin(mCrackAngle));
    double cosAlpha(cos(mCrackAngle));
    double gamma(1./std::max(fabs(sinAlpha),fabs(cosAlpha)));
    assert(mlX==mlY);
    double factor=gamma/(2.*mlX);

    rbHomU[0] = factor*sinAlpha;
    rbHomU[1] = 0.;
    rbHomU[2] = factor*(-cosAlpha);

    rbHomU[3] = 0.;
    rbHomU[4] = factor*(-cosAlpha);
    rbHomU[5] = factor*sinAlpha;

    if (fabs(sinAlpha)>fabs(cosAlpha))
    {
        rbHomAlpha[0] = 0.;
        rbHomAlpha[1] = mCrackOpening[1]/(2.*mlX*sinAlpha*fabs(sinAlpha));
        rbHomAlpha[2] = mCrackOpening[0]/(2.*mlX*sinAlpha*fabs(sinAlpha));

        if (rbHessian!=0)
        {
            rbHessian[0] = 0.;
            rbHessian[1] = -mCrackOpening[1]*cosAlpha/(mlX*sinAlpha*sinAlpha*fabs(sinAlpha));
            rbHessian[2] = -mCrackOpening[0]*cosAlpha/(mlX*sinAlpha*sinAlpha*fabs(sinAlpha));

            rbHessian[3] = 0.;
            rbHessian[4] = 0.;
            rbHessian[5] = 1/(2.*mlX*sinAlpha*fabs(sinAlpha));

            rbHessian[6] = 0.;
            rbHessian[7] = 1./(2.*mlX*sinAlpha*fabs(sinAlpha));
            rbHessian[8] = 0.;

        }
    }
    else
    {
        rbHomAlpha[0] = mCrackOpening[0]/(2.*mlX*cosAlpha*fabs(cosAlpha));
        rbHomAlpha[1] = 0.;
        rbHomAlpha[2] = mCrackOpening[1]/(2.*mlX*cosAlpha*fabs(cosAlpha));
        if (rbHessian!=0)
        {
            rbHessian[0] = mCrackOpening[0]*sinAlpha/(mlX*cosAlpha*cosAlpha*fabs(cosAlpha));
            rbHessian[1] = 0.;
            rbHessian[2] = mCrackOpening[1]*sinAlpha/(mlX*cosAlpha*cosAlpha*fabs(cosAlpha));

            rbHessian[3] = 1./(2.*mlX*cosAlpha*fabs(cosAlpha));
            rbHessian[4] = 0.;
            rbHessian[5] = 0.;

            rbHessian[6] = 0.;
            rbHessian[7] = 0.;
            rbHessian[8] = 1./(2.*mlX*cosAlpha*fabs(cosAlpha));
        }
    }
}

//! @brief add constraint equation for alpha in case of norm of crackopening less than a prescribed value
//! @return constraint equation
void NuTo::StructureIp::SetConstraintAlpha(EngineeringStrain2D& rTotalEngineeringStrain)
{
    if (mConstraintAlpha!=-1)
        this->mConstraintMap.erase(mConstraintAlpha);
    if (mCrackOpening[0]*mCrackOpening[0] + mCrackOpening[1]*mCrackOpening[1] > tolCrackOpeningConstraint)
        return;

    // find new constraint equation number
    mConstraintAlpha = 0;
    while (this->mConstraintMap.find(mConstraintAlpha) != this->mConstraintMap.end())
    {
        mConstraintAlpha++;
    }

    // create new constraint equation term
    ConstraintBase* constraintPtr = new ConstraintLinearGlobalCrackAngle(this);

    // insert constraint equation into map
    this->mConstraintMap.insert(mConstraintAlpha, constraintPtr);
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
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::StructureIp)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
