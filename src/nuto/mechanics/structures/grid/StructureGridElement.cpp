#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/Voxel8N.h"
#include <sstream>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumElements() const
{
    return mElementVec.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::ElementBase* NuTo::StructureGrid::ElementGetElementPtr(int rElementNumber)
{
    if (rElementNumber<0 || rElementNumber>=GetNumElements())
        throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] Conversion from string to int did not yield valid element number.");
    return &mElementVec[rElementNumber];
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
const NuTo::ElementBase* NuTo::StructureGrid::ElementGetElementPtr(int rElementNumber) const
{
    if (rElementNumber<0 || rElementNumber>=GetNumElements())
        throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] Conversion from string to int did not yield valid element number.");
    return &mElementVec[rElementNumber];
}

//! @brief gives the identifier of a element
//! @param reference to a element
//! @return identifier
int NuTo::StructureGrid::ElementGetId(const ElementBase* rElement)const
{
    int elementNumber(0);
    boost::ptr_vector<ElementBase>::const_iterator it;
    for (it = mElementVec.begin(); it!= mElementVec.end(); it++,elementNumber++)
    {
        if (&(*it)==rElement)
        {
            break;
        }
    }
    if (it==mElementVec.end())
        throw MechanicsException("[NuTo::StructureGrid::GetElementId] Element does not exist.");
    return elementNumber;
}

//! @brief info about the elements in the Structure
void NuTo::StructureGrid::ElementInfo(int mVerboseLevel)const
{
    std::cout<<"number of elements: " << mElementVec.size() <<std::endl;
}

//! @brief create element grid without data free elements
//! @param reference to a base coefficient matrix, to a ColorToMaterialMatrix and to an element type
void NuTo::StructureGrid::CreateElementGrid( NuTo::SparseMatrixCSRGeneral<double>& rBaseCoefficientMatrix0,
const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType)
{
    unsigned int numElements=0;       //counter for created elements
    NuTo::FullMatrix<int> imageValues (mNumVoxel,1);         //Color value for each voxel
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    typedef NuTo::SparseMatrixCSRGeneral<double> sparseMat ;
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixHelp ;
    std::vector<sparseMat> stiffnessMatrix(1);
    int numCoeffMat=0;   //material counter
    std::vector<double> youngsModulus(1);
    youngsModulus[0]=0;
    int matFlag=1;
    if (rColorToMaterialData.GetNumColumns()==1) //Only Young's modulus as changing parameter
    {
        for(unsigned int countVoxels =0; countVoxels<mNumVoxel;countVoxels++)//countVoxels correspond to VoxelID
        {
            if (rColorToMaterialData(imageValues(countVoxels,0),0)>0.1) //if Modul is> zero
            {
                if (youngsModulus[0]==0) //no coefficient matrix yet
                {
                    //set youngsModulus and add on material on counter
                    youngsModulus[numCoeffMat]=rColorToMaterialData(imageValues(countVoxels,0),0);
                    //stiffnessMatrixHelp = rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat];
                    stiffnessMatrix[numCoeffMat] = rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat];
                    NuTo::StructureGrid::ElementCreate(stiffnessMatrix[numCoeffMat],numElements,countVoxels,rElementType);//element number, element id, attr.
                    numElements++;
                    numCoeffMat++;
                    matFlag=0; //matrix already added
                }
                else
                {
                    for(int countMat=0;countMat<numCoeffMat;countMat++)
                    {
                        if (rColorToMaterialData(imageValues(countVoxels,0),0)==youngsModulus[countMat]) //same modulus already used
                        {
                            NuTo::StructureGrid::ElementCreate(stiffnessMatrix[countMat],numElements,countVoxels,rElementType);//element number, element id, attr.
                            numElements++;
                            countMat=numCoeffMat;
                            matFlag=0; //do not add a new coefficient matrix
                        }
                    }
                }
                if (matFlag==1)    //no equal coeff matrix found, create new
                {
                    youngsModulus.push_back(rColorToMaterialData(imageValues(countVoxels,0),0));
                    stiffnessMatrix.push_back(rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat]) ;
                    NuTo::StructureGrid::ElementCreate(stiffnessMatrix[numCoeffMat],numElements,countVoxels,rElementType);//element number, element id, attr.
                    numElements++;
                    numCoeffMat++;
                }
            }
            matFlag=1; //initialize matFlag for next step
      }
    }

}
//! @TODO ElementCreate without elementNumber

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::StructureGrid::ElementCreate ( NuTo::SparseMatrixCSRGeneral<double>& rCoefficientMatrix0,unsigned int rElementNumber, unsigned int rElementID, const std::string& rElementType)
{
    ElementCreate(rCoefficientMatrix0,rElementNumber,rElementID,rElementType,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA"));
}
//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::StructureGrid::ElementCreate(NuTo::SparseMatrixCSRGeneral<double>& rCoefficientMatrix0,unsigned int rElementNumber, unsigned int rElementID,  const std::string& rElementType,
        const std::string& rElementDataType, const std::string& rIpDataType)
{
    // get element type
    std::string upperCaseElementType;
    std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

    Element::eElementType elementType;
    if (upperCaseElementType=="VOXEL8N")
    {
        elementType = NuTo::Element::VOXEL8N;
    }
    else
    {
        throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Element type "+upperCaseElementType +" does not exist for structured grids.");
    }
    // check element data
     std::string upperCaseElementDataType;
     std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

     ElementData::eElementDataType elementDataType;
     if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
     {
     	elementDataType = NuTo::ElementData::CONSTITUTIVELAWIP;
     }
     else if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
    	{
 		elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL;
    	}
 	else
 	{
 		throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Element data type "+upperCaseElementDataType +" does not exist for structured grids.");
 	}

     // check ip data
     std::string upperCaseIpDataType;
     std::transform(rIpDataType.begin(), rIpDataType.end(), std::back_inserter(upperCaseIpDataType), (int(*)(int)) toupper);

     IpData::eIpDataType ipDataType;
     if (upperCaseIpDataType=="NOIPDATA")
     {
     	ipDataType = NuTo::IpData::NOIPDATA;
     }
     else if (upperCaseIpDataType=="STATICDATA")
    	{
     	ipDataType = NuTo::IpData::STATICDATA;
    	}
     else if (upperCaseIpDataType=="STATICDATANONLOCAL")
    	{
     	ipDataType = NuTo::IpData::STATICDATANONLOCAL;
    	}
 	else
 	{
 		throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Element data type "+upperCaseIpDataType +" does not exist for structured grids.");
 	}

     this->ElementCreate(rCoefficientMatrix0,rElementNumber, rElementID, elementType,elementDataType, ipDataType);
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes

void NuTo::StructureGrid::ElementCreate (NuTo::SparseMatrixCSRGeneral<double>& rCoefficientMatrix0,unsigned int rElementNumber,
		unsigned int rElementID, Element::eElementType rElementType, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
    // const IntegrationTypeBase *ptrIntegrationType;
    ElementBase* ptrElement;
    switch (rElementType)
    {
     case NuTo::Element::VOXEL8N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Voxel8N::GetStandardIntegrationType());
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Voxel8N is a 3D element.");
        }
        ptrElement = new NuTo::Voxel8N(this,rElementID,rCoefficientMatrix0,rElementDataType,rIpDataType);

        //ptrElement = new NuTo::Voxel8N(this,&rLocalCoefficientMatrix0, rElementDataType);
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::StructureGrid::ElementCreate] Invalid element type.");
    }
    this->mElementVec.push_back(ptrElement);
}

