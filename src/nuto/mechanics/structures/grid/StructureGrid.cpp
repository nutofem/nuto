// $Id$
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/elements/Voxel8N.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements3D.h"


NuTo::StructureGrid::StructureGrid(int rDimension) : StructureBase(rDimension)
{
	mVoxelLocation = 0;
	mDofIsNotConstraint = 0;
}

NuTo::StructureGrid::~StructureGrid()
{
	if(mVoxelLocation)
		delete mVoxelLocation;
	mVoxelLocation = 0;
//	if(mDofIsNotConstraint) //not needed
		delete [] mDofIsNotConstraint;
	mDofIsNotConstraint = 0;
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::StructureGrid::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::StructureGrid::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::StructureGrid::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::StructureGrid::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::StructureGrid::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::StructureGrid::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::StructureGrid::serialize(Archive & ar, const unsigned int version)
{
  //  & BOOST_SERIALIZATION_NVP(mDimension);
    std::cout << "start serialization of grid structure" << std::endl;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase)
       & BOOST_SERIALIZATION_NVP (mElementVec)
       & BOOST_SERIALIZATION_NVP (mNodeVec)
       & BOOST_SERIALIZATION_NVP (mNumVoxel)
       & BOOST_SERIALIZATION_NVP (mVoxelSpacing)
       & BOOST_SERIALIZATION_NVP (mGridDimension)
       & BOOST_SERIALIZATION_NVP (mNumMaterials)
       & BOOST_SERIALIZATION_NVP (mLocalCoefficientMatrix0)
       & BOOST_SERIALIZATION_NVP (mVoxelLocation);
     std::cout << "finish serialization of grid structure" << std::endl;
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureGrid::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::StructureGrid::Save] Error opening file.");
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
            throw MechanicsException ( "[NuTo::StructureGrid::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureGrid::Save]File save exception in boost - " ) + std::string ( e.what() ) );
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
void NuTo::StructureGrid::Restore (const std::string &filename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::StructureGrid::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureGrid::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureGrid::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
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
            throw MechanicsException ( "[NuTo::StructureGrid::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureGrid::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
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


//! @brief  store all elements of a structure in a vector
void NuTo::StructureGrid::GetElementsTotal(std::vector<ElementBase*>& rElements)
{
    boost::ptr_vector<ElementBase>::iterator ElementIter = this->mElementVec.begin();
    while (ElementIter != this->mElementVec.end())
    {
        rElements.push_back(&(*ElementIter));
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::StructureGrid::GetElementsTotal(std::vector<std::pair<int, ElementBase*> >& rElements)
{
	rElements.reserve(mElementVec.size());
    for (unsigned int id=0; id< this->mElementVec.size();id++)
    {
    	rElements.push_back(std::pair<int, ElementBase*>(id,&(mElementVec[id])));
    }
}

//! @brief  store all elements of a structure in a vector
void NuTo::StructureGrid::GetElementsTotal(std::vector<const ElementBase*>& rElements) const
{
    boost::ptr_vector<ElementBase>::const_iterator ElementIter = this->mElementVec.begin();
    while (ElementIter != this->mElementVec.end())
    {
        rElements.push_back(&(*ElementIter));
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::StructureGrid::GetElementsTotal(std::vector<std::pair<int, const ElementBase*> >& rElements) const
{
	rElements.reserve(mElementVec.size());
    for (unsigned int id=0; id< this->mElementVec.size();id++)
    {
    	rElements.push_back(std::pair<int, const ElementBase*>(id,&(mElementVec[id])));
    }
}

//! @brief  store all nodes of a structure in a vector
void NuTo::StructureGrid::GetNodesTotal(std::vector<NodeBase*>& rNodes)
{
	rNodes.reserve(mNodeVec.size());
    boost::ptr_vector<NodeBase>::iterator NodeIter = this->mNodeVec.begin();
    while (NodeIter != this->mNodeVec.end())
    {
    	rNodes.push_back(&(*NodeIter));
    	NodeIter++;
    }
}

// store all nodes of a structure in a vector
void NuTo::StructureGrid::GetNodesTotal(std::vector<std::pair<int, const NodeBase*> >& rNodes) const
{
	rNodes.reserve(mNodeVec.size());
    for (unsigned int id=0; id< this->mNodeVec.size();id++)
    {
    	rNodes.push_back(std::pair<int, const NodeBase*>(id,&(mNodeVec[id])));
    }
}


//! @brief  store all nodes of a structure in a vector
void NuTo::StructureGrid::GetNodesTotal(std::vector<const NodeBase*>& rNodes) const
{
	rNodes.reserve(mNodeVec.size());
    boost::ptr_vector<NodeBase>::const_iterator NodeIter = this->mNodeVec.begin();
    while (NodeIter != this->mNodeVec.end())
    {
    	rNodes.push_back(&(*NodeIter));
    	NodeIter++;
    }
}

// store all nodes of a structure in a vector
void NuTo::StructureGrid::GetNodesTotal(std::vector<std::pair<int,NodeBase*> >& rNodes)
{
	rNodes.reserve(mNodeVec.size());
    for (unsigned int id=0; id< this->mNodeVec.size();id++)
    {
    	rNodes.push_back(std::pair<int, NodeBase*>(id,&(mNodeVec[id])));
    }
}

//! @brief returns number of Voxels
//! @return number of Voxels
int NuTo::StructureGrid::GetNumVoxels() const
{
    return mNumVoxel;
}

//! @brief returns  VoxelSpacing
//! @return VoxelSpacing
const double* NuTo::StructureGrid::GetVoxelSpacing() const
{
    return mVoxelSpacing;
}

//! @brief returns GridOrigin
 //! @return GridOrigin
const double* NuTo::StructureGrid::GetGridOrigin() const
{
     return mGridOrigin;
}

 //! @brief returns GridDimension
 //! @return GridDimension
const int* NuTo::StructureGrid::GetGridDimension() const
{
     return mGridDimension;
}
//! @brief Get NumMaterials
//! @return NumMaterial
const int NuTo::StructureGrid::GetNumMaterials() const
{
    return mNumMaterials;
}

void NuTo::StructureGrid::ImportFromVtkASCIIFileHeader(const char* rFileName)
{
    try
    {
        using namespace boost::spirit::classic;
#define LocToGlob
        // open file
        std::ifstream file(rFileName, std::ios::in);
        if (file.is_open() == false)
        {
            throw MechanicsException("[StructureGrid::ImportFromVtkASCIIFile] error opening file.");
        }
        mImageDataFile=rFileName;

        // read header
        assert(file.is_open());

        // read first four lines
        std::string line;
        getline (file, line);
        getline (file, line);
        getline (file, line);
        getline (file, line);
        // read dimension
        getline (file, line);
        if (parse(line.c_str(),("DIMENSIONS " >> uint_p[assign_a(this->mGridDimension[0])] >> ' '
                                           >> uint_p[assign_a(this->mGridDimension[1])] >> ' '
                                           >> uint_p[assign_a(this->mGridDimension[2])] >> *space_p)).full == false)
       {
            throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading dimension.");
       }
        // read spacing
        getline (file, line);
        if (parse(line.c_str(),("SPACING ">> real_p[assign_a(this->mVoxelSpacing[0])] >> ' '
                                >> real_p[assign_a(this->mVoxelSpacing[1])] >> ' '
                                >> real_p[assign_a(this->mVoxelSpacing[2])] >> *space_p)).full == false)
        {
            throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading spacing.");
        }
        // read origin
        getline (file, line);
        if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(this->mGridOrigin[0])] >> ' '
                                 >> real_p[assign_a(this->mGridOrigin[1])] >> ' '
                                 >> real_p[assign_a(this->mGridOrigin[2])] >> *space_p)).full == false)
        {
             throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading origin.");
        }
        // read number of entries
        getline (file, line);
        if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(this->mNumVoxel)] >>  *space_p)).full == false)
               {
                   throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading number of entries.");
               }
        // read data type
        getline (file, line);
        // read empty line
        getline (file, line);

        // close file
       file.close();
       // test of saved data
       if (mNumVoxel<1)
       {
           throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader] error number of entries is negative or zero.");
       }
       if (mGridDimension[0]<1||mGridDimension[1]<1||mGridDimension[2]<1)
       {
           throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader] error dimension is negative or zero.");
       }
   }
    catch (MechanicsException &e)
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
}


//! @brief Get voxel number and location for all elements
//! @return FullMatrix columns elements, rows voxel number and number in x, y,z direction
NuTo::FullMatrix<int>* NuTo::StructureGrid::GetVoxelNumAndLocMatrix()
{
// EXIST matrix ?
// yes: return pointer
// no: allocate matrix, calculate matrix and return pointer
	if (mVoxelLocation)
		return mVoxelLocation;
	else
	{
		mVoxelLocation = new NuTo::FullMatrix<int>(GetNumElements(),4);
		CalculateVoxelNumAndLocMatrix(mVoxelLocation);
		return mVoxelLocation;
	}
}

//! @brief Calculate VoxeLNumAndLocMatrix
//! @return FullMatrix columns elements, rows voxel number and number in x, y,z direction
void NuTo::StructureGrid::CalculateVoxelNumAndLocMatrix(FullMatrix<int> *rVoxelLocation)
{
    std::cout<<__FILE__<<" "<<__LINE__<<"in CalculateVoxelNumAndLocMatrix()"<<std::endl;
	NuTo::Voxel8N* thisElement;
	int actElemNum =0;
	int loc[4];
	for(int dim2 =0; dim2<mGridDimension[2];++dim2)
	{
		loc[3]=dim2;

		for(int dim1 =0; dim1<mGridDimension[1];++dim1)
		{
			loc[2]=dim1;
			for(int dim0 =0; dim0<mGridDimension[0];++dim0)
			{
				//std::cout<<__FILE__<<" "<<__LINE__<<"in Schleife"<<std::endl;
				loc[1]=dim0;
				loc[0]=(dim0+1)+((dim1*mGridDimension[0])+1)+((dim2*mGridDimension[0]*mGridDimension[1])+1) -3;
				thisElement= static_cast<Voxel8N*>(ElementGetElementPtr(actElemNum));
				//std::cout<<__FILE__<<" "<<__LINE__<<" numElem - Voxel -  voxel "<< actElemNum <<"  "<<thisElement->GetVoxelID() <<" "<<loc[0]<<std::endl;
				if (thisElement->GetVoxelID() == loc[0])
				{
				    for (int count=0; count<4; ++count)
					{
						(*rVoxelLocation)(actElemNum,count) = loc[count];
					}
				    //std::cout<<__FILE__<<" "<<__LINE__<<" "<<rVoxelLocation(actElemNum,0) <<" "<<rVoxelLocation(actElemNum,1) <<" "<<rVoxelLocation(actElemNum,2) <<" "<<rVoxelLocation(actElemNum,3) <<std::endl;
					actElemNum++;
					if (actElemNum==GetNumElements())
					{
						dim0=mGridDimension[0];
						dim1=mGridDimension[1];
						dim2=mGridDimension[2];
					}
				}
			}
		}
	}

}

//! @brief Get voxels corner numbers from bottom to top counter-clockwise
//! @return array of number of corners with corner numbers
void NuTo::StructureGrid::GetCornersOfVoxel(int rElementNumber,int *rVoxLoc, int *corners)
{
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::GetCornerOfVoxel] error dimension must be 3.");
 //   int corners[8];
	if (mVerboseLevel>4)
		std::cout<<__FILE__<<" "<<__LINE__<<"  " <<" voxel location x: "<<  rVoxLoc[1]<<" y: "  <<rVoxLoc[2]<<" z: "<<rVoxLoc[3]<<" nbr: "<<rVoxLoc[0]<<std::endl;

	corners[0] = rVoxLoc[3]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[2]     * (mGridDimension[1]+1) + rVoxLoc[1];
	corners[1] = rVoxLoc[3]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[2]     * (mGridDimension[1]+1) + rVoxLoc[1]+1;
	corners[2] = rVoxLoc[3]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[2]+1) * (mGridDimension[1]+1) + rVoxLoc[1] +1;
	corners[3] = rVoxLoc[3]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[2]+1) * (mGridDimension[1]+1) + rVoxLoc[1];

	corners[4] = (rVoxLoc[3]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[2]     * (mGridDimension[1]+1) + rVoxLoc[1];
	corners[5] = (rVoxLoc[3]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[2]     * (mGridDimension[1]+1) + rVoxLoc[1]+1;
	corners[6] = (rVoxLoc[3]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[2]+1) * (mGridDimension[1]+1) + rVoxLoc[1]+1;
	corners[7] = (rVoxLoc[3]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[2]+1) * (mGridDimension[1]+1) + rVoxLoc[1];
	if (mVerboseLevel>4)
	{
	    std::cout<<__FILE__<<" "<<__LINE__<<" Ecken des "<<rElementNumber<< ". Elements: "<< corners[0]<<"  "<<corners[1]<<"  "<<corners[2]<<"  "<<corners[3]<<"  "<<std::endl;
	    std::cout<<__FILE__<<" "<<__LINE__<<"  "<< corners[4]<<"  "<<corners[5]<<"  "<<corners[6]<<"  "<<corners[7]<<"  "<<std::endl;
	}
}

//! @brief Get LocalCoefficientMatrix0
//! @param NumLocalCoefficientMatrix0 number of stiffness matrix
NuTo::FullMatrix<double>* NuTo::StructureGrid::GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0)
{
 	if (rNumLocalCoefficientMatrix0<0 || rNumLocalCoefficientMatrix0>=GetNumMaterials())
        throw MechanicsException("[NuTo::StructureGrid::GetLocalCoefficientMatrix0] No valid material number.");
// 	std::cout<<__FILE__<<" "<<__LINE__<<"  "<< &mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0] <<std::endl;
// 	std::cout<<__FILE__<<" "<<__LINE__<<"  "<< mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0] <<std::endl;
    return &mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0];
}

void NuTo::StructureGrid::BuildLocalCoefficientMatrix0() const
{
// calculate stiffness for given Modulus
    throw MechanicsException("[StructureGrid::BuildLocalCoefficientMatrix0] method is not yet implemented.");
}

//! @brief  based on the global dofs build submatrices of the global coefficent matrix0
void NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
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
    boost::ptr_vector<ElementBase>::const_iterator elementIter = this->mElementVec.begin();
    while (elementIter != this->mElementVec.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
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
void NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK) const
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
    boost::ptr_vector<ElementBase>::const_iterator elementIter = this->mElementVec.begin();
    while (elementIter != this->mElementVec.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
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
void NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK) const
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
    boost::ptr_vector<ElementBase>::const_iterator elementIter = this->mElementVec.begin();
    while (elementIter != this->mElementVec.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());
        if(symmetryFlag == false)
        {
            throw MechanicsException("[NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0Symmetric] element matrix is not symmetric (general sparse matrix required).");
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
void NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKK) const
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
    boost::ptr_vector<ElementBase>::const_iterator elementIter = this->mElementVec.begin();
    while (elementIter != this->mElementVec.end())
    {
        // calculate element contribution
        bool symmetryFlag = false;
        elementIter->CalculateCoefficientMatrix_0(elementMatrix, elementMatrixGlobalDofsRow, elementMatrixGlobalDofsColumn, symmetryFlag);
        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementMatrixGlobalDofsRow.size());
        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementMatrixGlobalDofsColumn.size());
        if(symmetryFlag == false)
        {
            throw MechanicsException("[NuTo::StructureGrid::BuildGlobalCoefficientSubMatrices0Symmetric] element matrix is not symmetric (general sparse matrix required).");
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

void NuTo::StructureGrid::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
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
    boost::ptr_vector<ElementBase>::const_iterator elementIter = this->mElementVec.begin();
    while (elementIter != this->mElementVec.end())
    {
        // calculate element contribution
        elementIter->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
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
