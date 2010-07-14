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

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/elements/ElementBase.h"

#include <ANN/ANN.h>
#include <set>

//! @brief constructor
//! @param mDimension number of nodes
NuTo::Structure::Structure ( int rDimension )  : StructureBase ( rDimension )
{
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Structure::Info()const
{
    StructureBase::Info();
    NodeInfo(mVerboseLevel);
    ElementInfo(mVerboseLevel);
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Structure::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Structure::serialize(Archive & ar, const unsigned int version)
{
    std::cout << "start serialization of structure" << std::endl;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase)
       & boost::serialization::make_nvp ("elementMap", mElementMap)
       & boost::serialization::make_nvp ("nodeMap", mNodeMap);
    std::cout << "finish serialization of structure" << std::endl;
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::Structure::Save] Error opening file.");
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
            throw MechanicsException ( "[NuTo::Structure::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::Structure::Save]File save exception in boost - " ) + std::string ( e.what() ) );
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
}

//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Restore (const std::string &filename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::Structure::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Structure::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Structure::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Structure::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::Structure::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::Structure::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
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
}
#endif // ENABLE_SERIALIZATION


// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::Structure::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementMatrix;
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
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
        }
        elementIter++;
    }
}

// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::Structure::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
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

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
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
        }
        elementIter++;
    }
}

// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::Structure::BuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == true);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementMatrix;
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());
        if(symmetryFlag == false)
        {
            throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatrices0Symmetric] element matrix is not symmetric (general sparse matrix required).");
        }

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
                        // add upper triangle and diagonal
                        if(globalColumnDof >= globalRowDof)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                    }
                    else
                    {
                        rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                    }
                }
            }
        }
        elementIter++;
    }
}

// based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::Structure::BuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKK) const
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == true);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);
    assert(rMatrixKK.IsSymmetric() == true);
    assert(rMatrixKK.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumEntries() == 0);

    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementMatrix;
    std::vector<int> elementMatrixGlobalDofsRow;
    std::vector<int> elementMatrixGlobalDofsColumn;

    // loop over all elements
    boost::ptr_map<int,ElementBase>::const_iterator elementIter = this->mElementMap.begin();
    while (elementIter != this->mElementMap.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->second->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());
        if(symmetryFlag == false)
        {
            throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatrices0Symmetric] element matrix is not symmetric (general sparse matrix required).");
        }

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
                        // add upper triangle and diagonal
                        if(globalColumnDof >= globalRowDof)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
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
                    if (globalColumnDof >= this->mNumActiveDofs)
                    {
                        // add upper triangle and diagonal
                        if(globalColumnDof >= globalRowDof)
                        {
                            rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
        elementIter++;
    }
}

void NuTo::Structure::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
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
        }
        elementIter++;
    }
}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(int rConstitutiveId)
{
    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveId);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::Structure::BuildNonlocalData] Constitutive law with the given identifier does not exist.");

    try
    {
    	BuildNonlocalData(itConstitutive->second);
    }
    catch(NuTo::MechanicsException e)
    {
        e.AddMessage("[NuTo::Structure::BuildNonlocalData] Error calculating nonlocal data.");
        throw e;
    }
    catch(...)
    {
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error calculating nonlocal data.");
    }
}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(const ConstitutiveBase* rConstitutive)
{
	double R(rConstitutive->GetNonlocalRadius());
	double R2(R*R);
    std::vector<ElementBase*> indexElement;
    std::vector<int> indexIp;
    std::vector<double> indexIpVolume;


	// build up search tree with all integration points
    boost::ptr_map<int,ElementBase>::iterator elementIter;
    for (elementIter = this->mElementMap.begin(); elementIter!= this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        if (elementPtr==0)
        	continue;

//TODO check element type
//TODO check constitutive type
        std::vector<double> ipVolume;
        elementPtr->GetIntegrationPointVolume(ipVolume);

        //calculate element contribution and iterate over all integration points
		for (int theIp = 0; theIp<elementPtr->GetNumIntegrationPoints();theIp++)
		{
			//theWeight = elementIter->second->GetGlobalIntegrationPointWeight(theIp);
			if (elementPtr->GetConstitutiveLaw(theIp)==rConstitutive)
			{
				indexElement.push_back(elementPtr);
				indexIp.push_back(theIp);
				indexIpVolume.push_back(ipVolume[theIp]);
			}
    	}
    }

    // build kd_tree
    ANNpointArray dataPoints;
    ANNkd_tree*   kdTree;

    if(mDimension == 2)
    {
    	dataPoints = annAllocPts(indexIp.size(),2);
        for(unsigned int count = 0; count < indexIp.size(); count++)
        {
            ANNpoint thePoint = dataPoints[count];
            double coordinates[3];
            indexElement[count]->GetGlobalIntegrationPointCoordinates(indexIp[count],coordinates);
            //the third parameter is probably zero, but in order to avoid writing another routine ...
            thePoint[0] = coordinates[0];
            thePoint[1] = coordinates[1];
        }
        kdTree = new ANNkd_tree(dataPoints,indexIp.size(),2);
    }
    else
    {
    	dataPoints = annAllocPts(indexIp.size(),3);
        for(unsigned int count = 0; count < indexIp.size(); count++)
        {
            ANNpoint thePoint = dataPoints[count];
            double coordinates[3];
            indexElement[count]->GetGlobalIntegrationPointCoordinates(indexIp[count],coordinates);
            thePoint[0] = coordinates[0];
            thePoint[1] = coordinates[1];
            thePoint[2] = coordinates[2];
        }
        kdTree = new ANNkd_tree(dataPoints,indexIp.size(),2);
    }

    // find the neighbors in radius R
    unsigned int allocatedResultPoints = 100;
    ANNidxArray nnIdx = new ANNidx[allocatedResultPoints];
    ANNdistArray dists = new ANNdist[allocatedResultPoints];

    for(unsigned int theIp = 0; theIp < indexIp.size(); theIp++)
    {
    	ElementBase* elementPtr = indexElement[theIp];
        int localIpNumber = indexIp[theIp];
        unsigned int numNeighborPoints = 0;
        do
        {
        	numNeighborPoints = kdTree->annkFRSearch(dataPoints[theIp],
                                                    R2,
                                                    allocatedResultPoints,
                                                    nnIdx,
                                                    dists,
                                                    0
                                                   );
            if(numNeighborPoints > allocatedResultPoints)
            {
                allocatedResultPoints = numNeighborPoints;
                delete [] nnIdx;
                delete [] dists;
                nnIdx = new ANNidx[allocatedResultPoints];
                dists = new ANNdist[allocatedResultPoints];
                numNeighborPoints = 0;
            }
        }
        while(numNeighborPoints == 0);

        //calculate total weight to scale the weights
        double totalVolume(0);
        for (unsigned int theNeighborPoint=0; theNeighborPoint<numNeighborPoints; theNeighborPoint++)
        {
        	if (dists[theNeighborPoint]<R2)
        	{
        		dists[theNeighborPoint] = (1-dists[theNeighborPoint]/R2);
        		dists[theNeighborPoint]*= dists[theNeighborPoint];
        		totalVolume+=indexIpVolume[nnIdx[theNeighborPoint]]*dists[theNeighborPoint];
        	}
        	else
        		dists[theNeighborPoint] = 0.;

        }
        totalVolume=1./totalVolume;

        //add all the nonlocal integration points to the nonlocal data of theIp
        for (unsigned int theNeighborPoint=0; theNeighborPoint<numNeighborPoints; theNeighborPoint++)
        {
        	int theNeighborIndex(nnIdx[theNeighborPoint]);
        	elementPtr->SetNonlocalWeight(localIpNumber, indexElement[theNeighborIndex],
        			indexIp[theNeighborIndex], indexIpVolume[theNeighborIndex]*dists[theNeighborPoint]*totalVolume);
        }
    }

    //just check the sum of the weights to be one
    for (elementIter = this->mElementMap.begin(); elementIter!= this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        if (elementPtr==0)
        	continue;

        for (int theIp=0; theIp< elementPtr->GetNumIntegrationPoints(); theIp++)
        {
			const std::vector<const NuTo::ElementBase*>& nonlocalElements(elementPtr->GetNonlocalElements());

			double sumW(0.);
			for (int countNonlocalElement=0; countNonlocalElement<(int)nonlocalElements.size(); countNonlocalElement++)
			{
				const std::vector<double>& weights(elementPtr->GetNonlocalWeights(theIp,countNonlocalElement));

				assert((int)weights.size()==nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints());

				//Go through all the integration points
				for (int theIP=0; theIP<(int)weights.size(); theIP++)
				{
					 //and add up to nonlocal equivalent plastic strain
					sumW+=weights[theIP];
				}
			}
        }
    }

    delete kdTree;
    annDeallocPts(dataPoints);

    delete [] nnIdx;
    delete [] dists;

    annClose();
/*
    if(EXIT_SUCCESS != 0)
    {
        INTERPRET_INTERN error_mess("ELEMENT_BUILD_NL_ELEMENTS_ANN: Error using ANN library.");
        return -1;
    }
*/
}
//! @brief import from gmsh
//! @param rFileName .. file name
//! @param rDOFs .. degrees of freedom for the nodes
//! @param rElementData .. element data for the elements to be created
//! @param rIPData .. ip data for the integration points to be created
void NuTo::Structure::ImportFromGmsh (const std::string& rFileName,
		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData)
{
    try
    {
    	std::set<int> groupIds;
    	ImportFromGmshAux(rFileName, rDOFs, rElementData, rIPData, false, groupIds);
    }
    catch(NuTo::MechanicsException e)
    {
        e.AddMessage("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
        throw e;
    }
    catch(...)
    {
    	throw NuTo::MechanicsException
    	   ("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
    }
}

void NuTo::Structure::ImportFromGmsh (const std::string& rFileName,
		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData,
		NuTo::FullMatrix<int>& rElementGroupIds)
{
    try
    {
    	std::set<int> groupIds;
    	ImportFromGmshAux(rFileName, rDOFs, rElementData, rIPData, true, groupIds);

    	rElementGroupIds.Resize(groupIds.size(),1);
    	int count(0);
    	for (std::set<int>::iterator it = groupIds.begin(); it != groupIds.end(); it++, count++)
    	{
    		rElementGroupIds(count,0) = *it;
    	}
    }
    catch(NuTo::MechanicsException e)
    {
        e.AddMessage("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
        throw e;
    }
    catch(...)
    {
    	throw NuTo::MechanicsException
    	   ("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
    }
}

class gmsh_node
{
public:
    gmsh_node():id(0)
    {
        this->Coordinates[0] = 0.;
        this->Coordinates[1] = 0.;
        this->Coordinates[2] = 0.;
    };
    unsigned int id;
    double Coordinates[3];
};

class gmsh_element
{
public:
    gmsh_element():id(0),type(0)
    {}
    ;
    unsigned int id;
    unsigned int type;
    std::vector<unsigned int> tags;
    std::vector<unsigned int> nodes;
};

//! @brief import from gmsh
//! @param rFileName .. file name
//! @param vector with the created groupes
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <iostream>
#include <string>

void NuTo::Structure::ImportFromGmshAux (const std::string& rFileName,
		const std::string& rDOFs, const std::string& rElementData, const std::string& rIPData,
		bool rAddGroups, std::set<int>& rElementGroupIds)
{
    const unsigned int num_elm_nodes[20] =
        {
            0,
            2,
            3,
            4,
            4,
            8,
            6,
            5,
            3,
            6,
            9,
            10,
            27,
            18,
            14,
            1,
            8,
            20,
            15,
            13
        };

    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    using qi::double_;
    using qi::uint_;
    using qi::phrase_parse;
    using ascii::space;
    using boost::phoenix::ref;
    using boost::phoenix::push_back;


    unsigned int num_nodes = 0;
    unsigned int num_elements = 0;
    std::ifstream file (rFileName, std::ios::in );
    if (file.is_open() == false)
    {
        throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error opening input file for read access.");
    }

    std::vector<gmsh_node> nodes;
    std::vector<gmsh_element> elements;

    // read first line
    std::string line;
    getline (file, line);

    // get file format
    std::getline (file, line);
    unsigned int mayor_version, minor_version, binary, double_size;

    std::string::iterator itFirst(line.begin());
    bool match (phrase_parse(itFirst,line.end(),
    		 uint_[ref(mayor_version)=qi::_1] >>'.' >> uint_[ref(minor_version)=qi::_1] >>
    		 uint_[ref(binary)=qi::_1] >>
    		 uint_[ref(double_size)=qi::_1],
             ascii::space));

    if (!match || itFirst!=line.end())
    {
    	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error reading file information.");
    }

    //std::cout << "version " << mayor_version <<"," << minor_version  << " binary " << binary << " double_size " << double_size << std::endl;

    if(mayor_version != 2)
    {
    	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid file format.");
    }


    if(binary == 0) // read ASCII - file
    {
        getline (file, line);
        if(line != "$EndMeshFormat")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndMeshFormat not found.");
        }

        // begin node section
        getline (file, line);
        if(line != "$Nodes")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Nodes not found.");
        }

        // read number of nodes
        getline (file, line);
        itFirst=line.begin();
        match = phrase_parse(itFirst,line.end(), uint_[ref(num_nodes)=qi::_1], ascii::space);

        if(!match || itFirst!=line.end())
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of nodes.");
        }

        // read node data
        nodes.resize(num_nodes);
        for(unsigned int node_count = 0; node_count < num_nodes; node_count++)
        {
            getline (file, line);
            itFirst=line.begin();
            match = phrase_parse(itFirst,line.end(),
            		uint_[ref(nodes[node_count].id)=qi::_1] >>
            		double_[ref(nodes[node_count].Coordinates[0])=qi::_1] >>
            		double_[ref(nodes[node_count].Coordinates[1])=qi::_1] >>
            		double_[ref(nodes[node_count].Coordinates[2])=qi::_1]
            		, ascii::space);

            if(!match || itFirst!=line.end())
            {
            	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading node data.");
            }
        }

        // end node section
        getline (file, line);
        if(line != "$EndNodes")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndNodes not found.");
        }

        // begin element section
        getline (file, line);
        if(line != "$Elements")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Elements not found.");
        }

        // read number of elements
        getline (file, line);
        itFirst=line.begin();
        match = phrase_parse(itFirst,line.end(), uint_[ref(num_elements)=qi::_1], ascii::space);

        if(!match || itFirst!=line.end())
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of elements.");
        }

        // read element data
        elements.resize(num_elements);
        std::vector<unsigned int> tmp_elem_data;
        for(unsigned int element_count = 0; element_count < num_elements; element_count++)
        {
            // read data to vector
            getline (file, line);
            itFirst=line.begin();
            match = phrase_parse(itFirst,line.end(),
            		*(uint_[push_back(boost::phoenix::ref(tmp_elem_data), qi::_1)]), ascii::space);

            if(!match || itFirst!=line.end())
            {
            	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading element data.");
            }
            std::vector<unsigned int>::iterator iter = tmp_elem_data.begin();

            // get id
            elements[element_count].id = *iter;
            iter++;
            // get type
            elements[element_count].type = *iter;
            iter++;

            // get number of tags and tags
            unsigned int num_tags = *iter;
            iter++;
            elements[element_count].tags.resize(num_tags);
            for(unsigned int tag_count = 0; tag_count < num_tags; tag_count++)
            {
                elements[element_count].tags[tag_count] = *iter;
                iter++;
            }

            // get number of nodes and nodes
            unsigned int cur_num_elm_nodes = num_elm_nodes[elements[element_count].type];
            elements[element_count].nodes.resize(cur_num_elm_nodes);
            for(unsigned int elem_node_count = 0; elem_node_count < cur_num_elm_nodes; elem_node_count++)
            {
                elements[element_count].nodes[elem_node_count]= *iter;
                iter++;
            }

            // check size
            if(iter != tmp_elem_data.end())
            {
            	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] invalid number of element data.");
            }
            tmp_elem_data.clear();
        }

        // end element section
        getline (file, line);
        if(line != "$EndElements")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndElements not found.");
        }
    }
/*
    else // binary format
    {
        if(double_size != sizeof(double))
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid size of double.");
        }

        // close file and open as binary
        file.close();
        file.open(file_name , std::ios::in | std::ios::binary);
        if (file.is_open() == false)
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error opening input file for read access.");
        }

        // read the first two lines again
        getline (file, line);
        getline (file, line);

        // check size of integer
        int one;
        file.read((char *)&one,sizeof(int));
        if(one != 1)
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid binary format.");
        }
        file.seekg(1,std::ios::cur);

        getline (file, line);
        if(line != "$EndMeshFormat")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndMeshFormat not found.");
        }

        // begin node section
        getline (file, line);
        if(line != "$Nodes")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Nodes not found.");
        }

        // read number of nodes
        getline (file, line);
        if(parse(line.c_str(),(uint_p[assign_a(num_nodes)]),space_p).full == false)
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of nodes.");
        }

        // read node data
        node_data.resize(num_nodes);
        for(unsigned int node_count = 0; node_count < num_nodes; node_count++)
        {
            file.read((char *)&node_data[node_count].id,sizeof(int));
            file.read((char *)node_data[node_count].coordinates, 3 * sizeof(double));
        }
        file.seekg(1,std::ios::cur);

        // end node section
        getline (file, line);
        if(line != "$EndNodes")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndNodes not found.");
        }

        // begin element section
        getline (file, line);
        if(line != "$Elements")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Elements not found.");
        }

        // read number of elements
        getline (file, line);
        if(parse(line.c_str(),(uint_p[assign_a(num_elements)]),space_p).full == false)
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of elements.");
        }
        element_data.resize(num_elements);

        // read element data
        unsigned int element_count = 0;
        do
        {
            // read element header
            unsigned int element_type;
            file.read((char *)&element_type,sizeof(unsigned int));
            unsigned int num_elm_follow;
            file.read((char *)&num_elm_follow,sizeof(unsigned int));
            unsigned int num_tags;
            file.read((char *)&num_tags,sizeof(unsigned int));

            unsigned int cur_num_elm_nodes = num_elm_nodes[element_type];
            for(unsigned int elm_count = 0; elm_count < num_elm_follow; elm_count++)
            {
                file.read((char *)&element_data[element_count].id,sizeof(int));
                element_data[element_count].type = element_type;
                element_data[element_count].tags.resize(num_tags);
                for(unsigned int tag_count = 0; tag_count < num_tags; tag_count++)
                {
                    file.read((char *)&element_data[element_count].tags[tag_count],sizeof(int));
                }
                element_data[element_count].nodes.resize(cur_num_elm_nodes);
                for(unsigned int elem_node_count = 0; elem_node_count < cur_num_elm_nodes; elem_node_count++)
                {
                    file.read((char *)&element_data[element_count].nodes[elem_node_count],sizeof(int));
                }
                element_count++;
            }
            break;
        }
        while(element_count != num_elements);
        file.seekg(1,std::ios::cur);

        // end element section
        getline (file, line);
        if(line != "$EndElements")
        {
        	throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndElements not found.");
        }
    }
*/

    //create the nodes
	NuTo::FullMatrix<double> coordinates;
	switch (mDimension)
	{
	case 2:
		coordinates.Resize(2,1);
		break;
	case 3:
		coordinates.Resize(3,1);
		break;
	default:
		throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Only implemented for 2D and 3D.");
	}
	std::map<int,int> newNodeNumber;
    for (unsigned int nodeCount=0; nodeCount<nodes.size(); nodeCount++)
    {
    	coordinates(0,0) = nodes[nodeCount].Coordinates[0];
    	coordinates(1,0) = nodes[nodeCount].Coordinates[1];
    	if (mDimension==3)
        	coordinates(2,0) = nodes[nodeCount].Coordinates[2];
    	newNodeNumber[ nodes[nodeCount].id] = NodeCreate(rDOFs, coordinates);
    }

	NuTo::FullMatrix<int> nodeNumbers;
    for (unsigned int elementCount=0; elementCount<elements.size(); elementCount++)
    {
    	nodeNumbers.Resize(elements[elementCount].nodes.size(),1);
    	for (unsigned int countNode=0; countNode< elements[elementCount].nodes.size(); countNode++)
    		nodeNumbers(countNode,0) = newNodeNumber[elements[elementCount].nodes[countNode]];
    	//std::cout << "element " << elementCount << " with nodes " << nodeNumbers.Trans() << std::endl;
    	int theElementId(-1);
    	switch (elements[elementCount].type)
    	{
    	case 2:
    		theElementId = ElementCreate("PLANE2D3N",nodeNumbers,rElementData, rIPData);
    		break;
    	case 8:
    		theElementId = ElementCreate("PLANE2D6N",nodeNumbers,rElementData, rIPData);
    		break;
    	default:
    		throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Element type not implemented in the import routine.");
    	}
    	if (rAddGroups)
    	{
			//add groups
			boost::ptr_map<int,GroupBase>::iterator itGroupMap(mGroupMap.find(elements[elementCount].tags[0]));
			if (itGroupMap==mGroupMap.end())
			{
				//create the element group
				GroupCreate(elements[elementCount].tags[0],NuTo::Groups::Elements);
			}
			rElementGroupIds.insert(elements[elementCount].tags[0]);
			GroupAddElement(elements[elementCount].tags[0],theElementId);
    	}
    }
}

