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

//! @brief gives the identifier of a node
//! @param reference to a node
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
void NuTo::StructureGrid::CreateElementGrid(const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType)
{
    int numVoxel=mGridDimension[0]*mGridDimension[1]*mGridDimension[2];
    unsigned int numElements;
    NuTo::FullMatrix<int> imageValues (numVoxel,1);
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    if(rColorToMaterialData.GetNumRows()==1)
    {
        for(int countVoxels =0; countVoxels<numVoxel;countVoxels++)//countVoxels correspond to VoxelID
        {
             std::cout<<"in CreateElementGrid routine!"<<std::endl;
             if (rColorToMaterialData(imageValues(countVoxels,0),0)>0.1) //if Modul is zero
             {
                 NuTo::StructureGrid::ElementCreate(numElements,countVoxels,rElementType);//element number, element id, attr.
                 numElements++;
                 std::cout<<"Element erstellt Nr.:"<<numElements-1<<std::endl;
             }
         }
    }
}
//! @TODO ElementCreate without elementNumber

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::StructureGrid::ElementCreate (unsigned int rElementNumber, unsigned int rElementID, const std::string& rElementType)
{
    ElementCreate(rElementNumber,rElementID,rElementType,std::string("CONSTITUTIVELAWELEMENT_NOSTATICDATA"));
}
//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::StructureGrid::ElementCreate(unsigned int rElementNumber, unsigned int rElementID,  const std::string& rElementType,
        const std::string& rElementDataType)
{
    // get element type
    std::string upperCaseElementType;
    std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

    ElementBase::eElementType elementType;
    if (upperCaseElementType=="VOXEL8N")
    {
        elementType = NuTo::ElementBase::VOXEL8N;
    }
    else
    {
        throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Element type "+upperCaseElementType +" does not exist for structured grids.");
    }
    // check element data
     std::string upperCaseElementDataType;
     std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

     NuTo::ElementDataBase::eElementDataType elementDataType;
     if (upperCaseElementDataType=="CONSTITUTIVELAWELEMENT_NOSTATICDATA")
     {
         elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_NOSTATICDATA;
     }
     else if (upperCaseElementDataType=="CONSTITUTIVELAWELEMENT_STATICDATA")
     {
         elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_STATICDATA;
     }
     else if (upperCaseElementType=="CONSTITUTIVELAWIP_NOSTATICDATA")
     {
         elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWIP_NOSTATICDATA;
     }
     else if (upperCaseElementType=="CONSTITUTIVELAWIP_STATICDATA")
     {
         elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWIP_STATICDATA;
     }
     else
     {
         throw MechanicsException("[NuTo::Structure::ElementCreate] Element data type "+upperCaseElementDataType +" does not exist.");
     }

     this->ElementCreate(rElementNumber, rElementID, elementType,elementDataType);
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes

void NuTo::StructureGrid::ElementCreate (unsigned int rElementNumber, unsigned int rElementID, ElementBase::eElementType rElementType, ElementDataBase::eElementDataType rElementDataType)
{
    // const IntegrationTypeBase *ptrIntegrationType;
    ElementBase* ptrElement;
    switch (rElementType)
    {
     case NuTo::ElementBase::VOXEL8N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Voxel8N::GetStandardIntegrationType());
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Voxel8N is a 3D element.");
        }
        ptrElement = new NuTo::Voxel8N(this,rElementID, rElementDataType);
        //ptrElement = new NuTo::Voxel8N(this,rElementDataType);
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::StructureGrid::ElementCreate] Invalid element type.");
    }
    this->mElementVec.push_back(ptrElement);


}

