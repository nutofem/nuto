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

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"

#include <ANN.h>

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

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<const ElementBase*>& rElements) const
{
    boost::ptr_map<int,ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<ElementBase*>& rElements)
{
    boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

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
    std::vector<ElementWithDataBase*> indexElement;
    std::vector<int> indexIp;
    std::vector<double> indexIpVolume;


	// build up search tree with all integration points
    boost::ptr_map<int,ElementBase>::iterator elementIter;
    for (elementIter = this->mElementMap.begin(); elementIter!= this->mElementMap.end(); elementIter++)
    {
        ElementWithDataBase* elementPtr = dynamic_cast<ElementWithDataBase*>(elementIter->second);
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
    	std::cout << "ip " << theIp << std::endl;
    	ElementWithDataBase* elementPtr = indexElement[theIp];
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
        	elementPtr->SetNonlocalWeight(localIpNumber, rConstitutive, indexElement[theNeighborIndex],
        			indexIp[theNeighborIndex], dists[theNeighborPoint]*totalVolume);
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
