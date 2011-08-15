// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/elements/Voxel8N.h"
#include <sstream>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumElements() const
{
    return mElementVec.size();
}

//! @brief a reference to a element
//! @param identifier
//! @return reference to a element
NuTo::Voxel8N* NuTo::StructureGrid::ElementGetPtr(int rElementNumber)
{
    if (rElementNumber<0 || rElementNumber>=GetNumElements())
    {
    	std::cout<<" ElementNumber "<<rElementNumber<<std::endl;
    	throw MechanicsException("[NuTo::StructureGrid::ElementGetPtr] Conversion from string to int did not yield valid element number.");
    }
    return mElementVec[rElementNumber];
}

//! @brief a reference to a element
//! @param identifier
//! @return reference to a element
const NuTo::Voxel8N* NuTo::StructureGrid::ElementGetPtr(int rElementNumber) const
{
    if (rElementNumber<0 || rElementNumber>=GetNumElements())
        throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] Conversion from string to int did not yield valid element number.");
    return mElementVec[rElementNumber];
}

//! @brief gives the identifier of a element
//! @param reference to a element
//! @return identifier
int NuTo::StructureGrid::ElementGetId(Voxel8N* rElement)
{
    int elementNumber(0);
    std::vector<Voxel8N*>::iterator it;
    for (it = mElementVec.begin(); it!= mElementVec.end(); it++,elementNumber++)
    {
        if ((*it)==rElement)
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
/*
//! @brief create element grid without data free elements
//! @param reference to a base coefficient matrix, to a ColorToMaterialMatrix and to an element type
void NuTo::StructureGrid::CreateElementGrid( NuTo::FullMatrix<double>& rBaseCoefficientMatrix0,
const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType)
{
    //int numElements=0;       //counter for created elements
    NuTo::FullMatrix<int> imageValues (mNumVoxel,1);         //Color value for each voxel
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    NuTo::StructureGrid::FullMat stiffnessMatrixHelp ;
    int numCoeffMat=0;   //material counter
    std::vector<double> youngsModulus(1);
    youngsModulus[0]=-1;
    int matExistsAlready=1;
    bool elemExists = true; //element exists
    if (rColorToMaterialData.GetNumColumns()!=1) //Only Young's modulus as changing parameter
      	throw MechanicsException("[NuTo::StructureGrid::CreateElementGrid] Wrong number of changing material parameters.");

	for(int countVoxels =0; countVoxels<mNumVoxel;countVoxels++)//countVoxels correspond to VoxelID
	{
		if (rColorToMaterialData(imageValues(countVoxels,0),0)>0) //if Modul is> zero
			elemExists = false;

		if (youngsModulus[0]==-1) //no coefficient matrix yet
		{
			//set youngsModulus and add on material on counter
			youngsModulus[numCoeffMat]=rColorToMaterialData(imageValues(countVoxels,0),0);
			std::cout<<__FILE__<<" " <<__LINE__<<" "<<" Young's modulus "<<numCoeffMat<<" : "<<youngsModulus[numCoeffMat]<<std::endl;
			if (!elemExists) //Young's Modulus equal zero
			{
				stiffnessMatrixHelp=NULL;
			}
			else
				stiffnessMatrixHelp = rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat];
			//std::cout<<__FILE__<<" " <<__LINE__<<" "<< stiffnessMatrixHelp <<std::endl;
		   this->mLocalCoefficientMatrix0.push_back(stiffnessMatrixHelp);
			NuTo::StructureGrid::ElementCreate(elemExists,numCoeffMat,countVoxels,rElementType);//mat number, element id, attr.
			//numElements++;
			numCoeffMat++;
			this->mNumMaterials = numCoeffMat;
			matExistsAlready=0; //matrix already added
		}
		else
		{
			for(int countMat=0;countMat<numCoeffMat;countMat++)
			{
				if (rColorToMaterialData(imageValues(countVoxels,0),0)==youngsModulus[countMat]) //same modulus already used
				{
					NuTo::StructureGrid::ElementCreate(elemExists,countMat,countVoxels,rElementType);//mat number, element id, attr.
					//numElements++;
					countMat=numCoeffMat;
					matExistsAlready=0; //do not add a new coefficient matrix
				}
			}
		}
		if (matExistsAlready==1)    //no equal coeff matrix found, create new
		{
			youngsModulus.push_back(rColorToMaterialData(imageValues(countVoxels,0),0));
			stiffnessMatrixHelp = rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat];
			this->mLocalCoefficientMatrix0.push_back(stiffnessMatrixHelp);
			NuTo::StructureGrid::ElementCreate(elemExists,numCoeffMat,countVoxels,rElementType);//mat number, element id, attr.
			//numElements++;
			numCoeffMat++;
			this->mNumMaterials = numCoeffMat;
		}
		matExistsAlready=1; //initialize matExistsAlready for next step
	}
}
*/
//! @brief create element grid without data free elements
//! @param reference to a base coefficient matrix, to a ColorToMaterialMatrix and to an element type
void NuTo::StructureGrid::CreateElementGrid( NuTo::FullMatrix<double>& rBaseCoefficientMatrix0,
const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType)
{
    //int numElements=0;       //counter for created elements
    NuTo::FullMatrix<int> imageValues (mNumVoxel,1);         //Color value for each voxel
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    NuTo::StructureGrid::FullMat stiffnessMatrixHelp ;
    int numCoeffMat=0;   //material counter
    std::vector<double> youngsModulus(1);
    youngsModulus[0]=0;
    bool matExistsAlready= true; //material exists
    bool elemExists = true; //element exists
    if (rColorToMaterialData.GetNumColumns()!=1) //Only Young's modulus as changing parameter
      	throw MechanicsException("[NuTo::StructureGrid::CreateElementGrid] Wrong number of changing material parameters.");

	for(int countVoxels =0; countVoxels<mNumVoxel;countVoxels++)//countVoxels correspond to VoxelID
	{
		if (rColorToMaterialData(imageValues(countVoxels,0),0)>0) //if Modul is> zero
			elemExists = true;
		else
		{
			elemExists = false;
			matExistsAlready = true; //no new stiffness matrix for non existent element
		}
		if (elemExists)
		{
			if (youngsModulus[0]==0) //no coefficient matrix yet
				matExistsAlready=false;
			else
			{
				for(int countMat=0;countMat<numCoeffMat;countMat++)
				{
					if (rColorToMaterialData(imageValues(countVoxels,0),0)==youngsModulus[countMat]) //same modulus already used
					{
						matExistsAlready = true;
						numCoeffMat=countMat;
					}
					else
						matExistsAlready=false;
				}
			}
		}

		if (!matExistsAlready)
		{
			//set youngsModulus and add on material on counter
			youngsModulus[numCoeffMat]=rColorToMaterialData(imageValues(countVoxels,0),0);
			stiffnessMatrixHelp = rBaseCoefficientMatrix0 * youngsModulus[numCoeffMat];
			std::cout<<__FILE__<<" " <<__LINE__<<" "<<" Young's modulus "<<numCoeffMat<<" : "<<youngsModulus[numCoeffMat]<<std::endl;
			//std::cout<<__FILE__<<" " <<__LINE__<<" "<< stiffnessMatrixHelp <<std::endl;
			this->mLocalCoefficientMatrix0.push_back(stiffnessMatrixHelp);
		}
		//std::cout<<" element create nr." << countVoxels<< "matnr "<<numCoeffMat<<std::endl;
		ElementCreate(elemExists,numCoeffMat,countVoxels,rElementType);//mat number, element id, attr.
		if (!matExistsAlready)
		{
			numCoeffMat++;
			this->mNumMaterials = numCoeffMat;
			matExistsAlready=true; //matrix already added
		}
	}
}
//! @TODO ElementCreate without elementNumber

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
void NuTo::StructureGrid::ElementCreate(bool flag, int rNumCoefficientMatrix0,int rElementID,  const std::string& rElementType)
{
    // get element type
    std::string upperCaseElementType;
    std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

    if (upperCaseElementType!="VOXEL8N")
    	throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Element type "+upperCaseElementType +" does not exist for structured grids.");
    this->ElementCreate(flag, rNumCoefficientMatrix0, rElementID);
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes

void NuTo::StructureGrid::ElementCreate (bool flag, int rNumCoefficientMatrix0,int rElementID)
{
    Voxel8N* ptrElement=0;
    if (flag)
    	ptrElement = new NuTo::Voxel8N(this,rNumCoefficientMatrix0);
   this->mElementVec.push_back(ptrElement);
}

//! @brief Deletes an element
//! @param rElementIdent identifier for the element
void NuTo::StructureGrid::ElementDelete(int rElementNumber)
{

	// @TODO [NuTo::StructureGrid::ElementDelete] has to be implemented
    throw MechanicsException("[NuTo::StructureGrid::ElementDelete] Not implemented yet!!!");

}



