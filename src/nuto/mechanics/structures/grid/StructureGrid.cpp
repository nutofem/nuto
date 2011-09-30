// $Id$
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif

# ifdef _OPENMP
    #include <omp.h>
# endif


#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/elements/Voxel8N.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements3D.h"

NuTo::StructureGrid::StructureGrid(int rDimension)
{
   if (rDimension!=3)
	{
		throw MechanicsException("[StructureGrid::StructureGrid] The dimension of the grid structure is either so far 3.");
	}
    mDimension = rDimension;
 	mVoxelLocation = 0;
	mDofIsNotConstraint = 0;
	mCalcVoxelLocation = 0;
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

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureGrid::Info()const
{
   std::cout << "dimension : " << mDimension << "\n";

   std::cout  << "num dofs : " << (mGridDimension[0]+1)*(mGridDimension[1]+1)* (mGridDimension[2]+1)<< "\n";
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
    std::cout << "start serialization of grid structure" << std::endl;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP (mDimension)
       & BOOST_SERIALIZATION_NVP (mElementVec)
       & BOOST_SERIALIZATION_NVP (mNodeVec)
       & BOOST_SERIALIZATION_NVP (mDimension)
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

//! @brief import routine for basic grid data without StructureGrid data space
void NuTo::StructureGrid::ImportFromVtkASCIIFileHeader(const char* rFileName,int *rGridDimension,double *rVoxelSpacing,double *rGridOrigin, int &rNumVoxel)
{
    std::cout<<__FILE__<<" "<<__LINE__<<" in ImportFromVtkASCIIFileHeader without StructureGrid data \n ";
    try
    {
        using namespace boost::spirit::classic;
        // open file
        std::ifstream file(rFileName, std::ios::in);
        if (file.is_open() == false)
        {
            throw MechanicsException("[StructureGrid::ImportFromVtkASCIIFile] error opening file.");
        }
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
       if (parse(line.c_str(),("DIMENSIONS " >> int_p[assign_a(rGridDimension[0])] >> ' '
                                           >> int_p[assign_a(rGridDimension[1])] >> ' '
                                           >> int_p[assign_a(rGridDimension[2])] >> *space_p)).full == false)
       {
            throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading dimension.");
       }
       // read spacing
        getline (file, line);
        if (parse(line.c_str(),("SPACING ">> real_p[assign_a(rVoxelSpacing[0])] >> ' '
                                >> real_p[assign_a(rVoxelSpacing[1])] >> ' '
                                >> real_p[assign_a(rVoxelSpacing[2])] >> *space_p)).full == false)
        {
            throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading spacing.");
        }
        // read origin
        getline (file, line);
        if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(rGridOrigin[0])] >> ' '
                                 >> real_p[assign_a(rGridOrigin[1])] >> ' '
                                 >> real_p[assign_a(rGridOrigin[2])] >> *space_p)).full == false)
        {
             throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading origin.");
        }
        // read number of entries
        getline (file, line);
        if (parse(line.c_str(),("POINT_DATA ">> int_p[assign_a(rNumVoxel)] >>  *space_p)).full == false)
	   {
		   throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader]error reading number of entries.");
	   }
        // close file
       file.close();
       // test of saved data
      if (rNumVoxel<1)
       {
           throw MechanicsException("[StructureGrid::importFromVtkASCIIFileReadHeader] error number of entries is negative or zero.");
       }
       if (rGridDimension[0]<1||rGridDimension[1]<1||rGridDimension[2]<1)
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


//! @brief import routine for basic grid data with StructureGrid data space
void NuTo::StructureGrid::ImportFromVtkASCIIFileHeader(const char* rFileName)
{
    try
    {
        using namespace boost::spirit::classic;
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

void NuTo::StructureGrid::ImportFromVtkASCIIFile(const char* rFileName,std::vector<int> &rData)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MechanicsException("[StructureGrid::ImportFromVtkASCIIFile] error opening file.");
    }
     // read first four lines
    unsigned int numEntries(0);

    std::string line;
    for (int count=0;count<7;count++)
    {
        getline (file, line);
    }

    // read number of entries
    getline (file, line);
    if (parse(line.c_str(),("POINT_DATA ">> int_p[assign_a(numEntries)] >>  *space_p)).full == false)
           {
               throw MechanicsException("[StructureGrid::importFromVtkASCIIFile]error reading number of entries.");
           }
    // read data type
    getline (file, line);
    // read empty line
    getline (file, line);

  // read entries
    if(numEntries!=rData.size())
        throw MechanicsException("[StructureGrid::importFromVtkASCIIFile] number of entries is not equal to vector size.");

    std::vector<int> imageValues(0);
    int value;

    while(getline(file,line))
    {
    	std::istringstream iss(line);
    	while(iss >> value)
    	{
			imageValues.push_back(value);
    	}
    }
    rData=imageValues;
     // close file
   file.close();
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



//! @brief create node data without StructureGrid
//! @brief set bool for node
void NuTo::StructureGrid::CreateGrid(int rThresholdMaterialValue, std::vector<int>& imageValues ,const std::vector<double>& rColorToMaterialData,int* rGridDimension,boost::dynamic_bitset<> &rNodeExist,boost::dynamic_bitset<> &rElemExist,std::vector<double>& youngsModulus,std::vector<int>& materialOfElem)
{
    int numGridNodes=(rGridDimension[0]+1)*(rGridDimension[1]+1)*(rGridDimension[2]+1);//all nodes of the grid
    //int numMatNodes=0;//all existing nodes (with material)
    int * coincidentVoxels=new int[8];
    bool flag=false;
    for (int countNodes =0; countNodes<numGridNodes;++countNodes)//countNodes correspond to nodeID
    {
         //get coincident voxels for each node, check if one voxel has material, then create node
         coincidentVoxels=GetCoincidenceVoxelIDs(countNodes);
         flag=false;
         for (int count =0; count<8; count++)
         {
        	 // voxel exist (for boundary nodes)
             if (coincidentVoxels[count]>-1)
             {
                 // voxel has material
            	 // color value 0 is material, 255 is air
            	 // material value smaller than thresholdvalue
				 std::cout<<"  "<<countNodes<<" voxel "<<coincidentVoxels[count]<<" image "<< imageValues[coincidentVoxels[count]];
                 if(imageValues[coincidentVoxels[count]]<rThresholdMaterialValue)
                 {
                     flag=true; //node exists
                     count=8;
                    rNodeExist.set(countNodes,true);
					 //numMatNodes++;
                 }
             }
         }
    }
	std::cout<<__FILE__<<" nodeExist " <<rNodeExist<<"\n";

    int numCoeffMat =0;
    bool matExistsAlready= false; //material exists
    //std::bitset<numVoxel> rElementExist; //0 = false, all 0 here
    int numVoxel=rGridDimension[0]*rGridDimension[1]*rGridDimension[2];//all nodes of the grid
  	for(int countVoxels =0; countVoxels<numVoxel;++countVoxels)//countVoxels correspond to VoxelID
 	{
 		if (rColorToMaterialData[imageValues[countVoxels]]>0) //if Modul is> zero
 		{
 			rElemExist.flip(countVoxels); //set to 1 = true, element exists

 			for(int countMat=0;countMat<numCoeffMat;countMat++)
			{
				if (rColorToMaterialData[imageValues[countVoxels]]==youngsModulus.at(countMat)) //same modulus already used
				{
					materialOfElem[countVoxels]=countMat;
					countMat=numCoeffMat+1;
					matExistsAlready=true;
				}
				else
					matExistsAlready=false;
			}
 			if (!matExistsAlready)
 	 		{
 	 			//set youngsModulus and add on material on counter
 	 			youngsModulus.push_back(rColorToMaterialData[imageValues[countVoxels]]);
 	  			numCoeffMat++;
 	 			matExistsAlready=true; //matrix already added
 	 		}
  		}
 		else
  			materialOfElem[countVoxels]=-1;
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
		CalculateVoxelLocations();
		return mVoxelLocation;
	}
}

//! @brief Calculate VoxelLocations for all elements
//! @param int[3] number in x, y,z direction
void NuTo::StructureGrid::CalculateVoxelLocations()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
	if(mVerboseLevel>3)
   std::cout<<__FILE__<<" "<<__LINE__<<"in  StructureGrid::CalculateVoxelLocations()"<< std::endl;
    int numElems= GetNumElements();
	mVoxelLocation = new NuTo::FullMatrix<int>(GetNumElements(),3);
	Voxel8N* thisElement=0;
    for (int element=0;element<numElems;++element)
    {
        thisElement= ElementGetPtr(element);
        if (thisElement)
        {
			int numDimxy=element/((mGridDimension[0])*(mGridDimension[1]));
			int numDimx=0;
			int residual1=element%((mGridDimension[0])*(mGridDimension[1]));
			int residual2=0;
			numDimx=residual1/(mGridDimension[0]);
			residual2=residual1%(mGridDimension[0]);
			int rVoxelLocation[3]={0};
			rVoxelLocation[0]=residual2;
			rVoxelLocation[1]=numDimx;
			rVoxelLocation[2]=numDimxy;
			thisElement->SetVoxelLocation(rVoxelLocation);
			//also for global vector
			(*mVoxelLocation)(element,0)=rVoxelLocation[0];
			(*mVoxelLocation)(element,1)=rVoxelLocation[1];
			(*mVoxelLocation)(element,2)=rVoxelLocation[2];
			if(mVerboseLevel>4)
				std::cout<<__FILE__<<" "<<__LINE__<<"loc "<<rVoxelLocation[0]<<", "<<rVoxelLocation[1]<<", "<<rVoxelLocation[2]<<std::endl;

        }
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
    	std::cout<<"[NuTo::StructureGrid::CalculateVoxelLocations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
    	std::cout<<"[NuTo::StructureGrid::CalculateVoxelLocations] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
}

//! @brief Get voxels corner numbers from bottom to top counter-clockwise
//! @return array of number of corners with corner numbers
void NuTo::StructureGrid::GetCornersOfVoxel(int rElementNumber,int *rVoxLoc, int *corners)
{
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::GetCornerOfVoxel] error dimension must be 3.");
 //   int corners[8];
	if (mVerboseLevel>4)
		std::cout<<__FILE__<<" "<<__LINE__<<"  " <<" voxel location x: "<<  rVoxLoc[0]<<" y: "  <<rVoxLoc[1]<<" z: "<<rVoxLoc[2]<<std::endl;

	corners[0] = rVoxLoc[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[1]     * (mGridDimension[1]+1) + rVoxLoc[0];
	corners[1] = rVoxLoc[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[1]     * (mGridDimension[1]+1) + rVoxLoc[0]+1;
	corners[2] = rVoxLoc[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[1]+1) * (mGridDimension[1]+1) + rVoxLoc[0] +1;
	corners[3] = rVoxLoc[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[1]+1) * (mGridDimension[1]+1) + rVoxLoc[0];

	corners[4] = (rVoxLoc[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[1]     * (mGridDimension[1]+1) + rVoxLoc[0];
	corners[5] = (rVoxLoc[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxLoc[1]     * (mGridDimension[1]+1) + rVoxLoc[0]+1;
	corners[6] = (rVoxLoc[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[1]+1) * (mGridDimension[1]+1) + rVoxLoc[0]+1;
	corners[7] = (rVoxLoc[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxLoc[1]+1) * (mGridDimension[1]+1) + rVoxLoc[0];
	if (mVerboseLevel>4)
	{
	    std::cout<<__FILE__<<" "<<__LINE__<<" Ecken des "<<rElementNumber<< ". Elements: "<< corners[0]<<"  "<<corners[1]<<"  "<<corners[2]<<"  "<<corners[3]<<"  "<<std::endl;
	    std::cout<<__FILE__<<" "<<__LINE__<<"  "<< corners[4]<<"  "<<corners[5]<<"  "<<corners[6]<<"  "<<corners[7]<<"  "<<std::endl;
	}
}
//! @brief NodeGetConstraintSwitch
//! @param rGlobalDof
//! @return switch for constraint
bool* NuTo::StructureGrid::GetConstraintSwitch()
{
	return mDofIsNotConstraint;
}


//! @brief Set NodeIds for all nodes at all elements
void NuTo::StructureGrid::SetAllNodeIds()
{
#ifdef SHOW_TIME
    //std::clock_t start,end;
    //start=clock();
	timespec startn,endn,diffn,diffsum;
	diffsum.tv_sec=0.;
	diffsum.tv_nsec=0.;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::GetCornerOfVoxel] error dimension must be 3.");
	Voxel8N* thisElement;
	int nodeIds[8];
	for (int element=0;element<GetNumElements();++element)
	{
		thisElement= ElementGetPtr(element);
		if (thisElement)
		{
			int * locVoxLoc=thisElement->GetVoxelLocation();
			//get grid corners of the voxel
			GetCornersOfVoxel(element, locVoxLoc, nodeIds);
			//std::cout<<"SetNodeIds : element "<<element<<"nodes "<<		nodeIds[0] <<" ,"<<	nodeIds[1]<<" ,"<<	nodeIds[2]<<" ,"<<	nodeIds[3] <<" ,"<<			nodeIds[4] <<" ,"<<		nodeIds[5]<<" ,"<<			nodeIds[6] <<" ,"<<		nodeIds[7]<<" ,"<<	std::endl;
			thisElement->SetNodeIds(nodeIds);
		}
	}
#ifdef SHOW_TIME
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
	diffn=diff(startn,endn);
	diffsum.tv_sec+=diffn.tv_sec;
	diffsum.tv_nsec+=diffn.tv_nsec;
	//end=clock();
	if (mShowTime)
        //std::cout<<"[NuTo::StructureGrid::SetAllNodeIds] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
        std::cout<<"[NuTo::StructureGrid::SetAllNodeIds] " << diffsum.tv_sec <<" sec: "<<diffsum.tv_nsec/1000000.<<" msec" << std::endl;
#endif
}

//! @brief Set NodeIds for all nodes at all nodes
void NuTo::StructureGrid::SetAllNodeIdsAtNode()
{
#ifdef SHOW_TIME
    //std::clock_t start,end;
    //start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::SetAllNodeIdsAtNode] error dimension must be 3.");
	NodeGrid3D* thisNode;
	Voxel8N* thisElement;
	int* localNodeIds=0;
	//local node ids sorted after element and node of element
	int orderNode[8][8]={
		{0,1,4,3,9,10,13,12},
		{1,2,5,4,10,11,14,13},
		{4,5,8,7,13,14,17,16},
		{3,4,7,6,12,13,16,15},
		{9,10,13,12,18,19,22,21},
		{10,11,14,13,19,20,23,22},
		{13,14,17,16,22,23,26,25},
		{12,13,16,15,21,22,25,24},
	};
	//field of 27 nodes which contain element and node of element
	std::vector<FullMatrix<double>*> partMatrix;
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in SetAllNodeIdsAtNode() "<<std::endl;

	for (int nodeNumber=0;nodeNumber<GetNumNodes();++nodeNumber)
	{
		thisNode=NodeGridGetNodePtr(nodeNumber);
		if (thisNode)
		{
			int nodeIds[27];
			for (int i=0;i<27;++i)
				nodeIds[i]=-1;
			int* elementIds=thisNode->GetElementIds();
			for (int element=0;element<8;++element)
			{
				if(elementIds[element]>-1)
				{
					thisElement= ElementGetPtr(elementIds[element]);
					if(thisElement)
					{
						localNodeIds=thisElement->GetNodeIds();
						for(int node=0;node<8;++node)
						{
							nodeIds[orderNode[element][node]]=localNodeIds[node];// overhead,set id 4 times, but slower with if clause
						}
					}
				}
			}
			thisNode->SetNodeIds(nodeIds);
		}
	}
#ifdef SHOW_TIME
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
    diffn=diff(startn,endn);
    //end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureGrid::SetAllNodeIdsAtNodes] " << diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec" <<std::endl;
#endif
}

//! @brief Set partCoefficientmatrix for all nodes
void NuTo::StructureGrid::SetAllPartCoefficientMatrix0()
{
#ifdef SHOW_TIME
    //std::clock_t start,end;
    //start=clock();
	timespec startn,endn,diffn;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::SetAllPartCoefficientMatrix0] error dimension must be 3.");
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in SetAllPartCoefficientMatrix0 "<<std::endl;
	// first field of 8 local nodes of one element
	// second field of 8 local nodes of one element
	// contains NBN ordering local node id
	int orderNode[8][8]={
		{13,14,17,16,22,23,26,25},
		{12,13,16,15,21,22,25,24},
		{9,10,13,12,18,19,22,21},
		{10,11,14,13,19,20,23,22},
		{4,5,8,7,13,14,17,16},
		{3,4,7,6,12,13,16,15},
		{0,1,4,3,9,10,13,12},
		{1,2,5,4,10,11,14,13},
	};

	for (int element=0;element<GetNumElements();++element)
	{
		Voxel8N* thisElement;
		thisElement= ElementGetPtr(element);
		if (thisElement)
		{
			int numLocalStiffness =thisElement->GetNumLocalStiffnessMatrix();
			FullMatrix<double>* matrix= GetLocalCoefficientMatrix0(numLocalStiffness);

			int* nodeIds=thisElement->GetNodeIds();
			NodeGrid3D* thisNodeI;
			NodeGrid3D* thisNodeJ;
			int k=0;
			for (int i=0;i<8;++i)
			{
				thisNodeI=NodeGridGetNodePtr(nodeIds[i]);
				if (!thisNodeI)
					std::cout<<" StrGrid  test  nod ptr print if 0"<<std::endl;
				for(int j=k;j<8;++j)
				{
					thisNodeJ=NodeGridGetNodePtr(nodeIds[j]);
					if (!thisNodeJ)
						std::cout<<" StrGrid  test  nod ptr print if 0"<<std::endl;
					//How does this work without copying a 3x3 matrix?
					//Adress of a part of the matrix not possible
					FullMatrix<double> partMatrix(3,3);
					partMatrix=(matrix->GetBlock(i*3,j*3,3,3));
					//thisNodeI->SetPartCoefficientMatrix0(orderNode[i][j],partMatrix);
					thisNodeI->SetPartCoefficient0(orderNode[i][j],partMatrix);
					if (i!=j)
					{
						partMatrix=(matrix->GetBlock(j*3,i*3,3,3));
						//thisNodeJ->SetPartCoefficientMatrix0(orderNode[j][i],partMatrix);
						thisNodeJ->SetPartCoefficient0(orderNode[j][i],partMatrix);
					}
				}
				++k;
			}
		}
	}

#ifdef SHOW_TIME
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
	diffn=diff(startn,endn);
	//end=clock();
	if (mShowTime)
	    std::cout<<"[NuTo::StructureGrid::SetAllPartCoefficientMatrix0] " << diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec" << std::endl;
        //std::cout<<"[NuTo::StructureGrid::SetAllPartCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

}

//! @brief Set ElementIds for the elements at each nodes
void NuTo::StructureGrid::SetAllElementIds()
{
#ifdef SHOW_TIME
    //std::clock_t start,end;
	//start=clock();
	timespec startn,endn,diffn,diffsum;
    diffsum.tv_sec=0.;
    diffsum.tv_nsec=0.;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
	if (mDimension != 3)
        throw MechanicsException("[StructureGrid::SetAllElementIds] error dimension must be 3.");
	NodeGrid3D* thisNode;
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<"in SetAllElementIds() "<<std::endl;
	int *elementIds=0;
	for (int node=0;node<GetNumNodes();++node)
	{
       thisNode=NodeGridGetNodePtr(node);
       if (thisNode)
       {
    	   elementIds=GetCoincidenceVoxelIDs(node);
    	   thisNode->SetElementIds(elementIds);
       }
	}
/*
#ifdef SHOW_TIME
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
    diffn=diff(startn,endn);
    diffsum.tv_sec+=diffn.tv_sec;
    diffsum.tv_nsec+=diffn.tv_nsec;
     if (mShowTime)
        std::cout<<"[NuTo::StructureGrid::SetAllElementIds: before element get id] " << diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec" << std::endl;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&startn);
#endif
*/
#ifdef SHOW_TIME
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&endn);
	diffn=diff(startn,endn);
    diffsum.tv_sec=diffn.tv_sec;
    diffsum.tv_nsec=diffn.tv_nsec;
	//end=clock();
	if (mShowTime)
		//std::cout<<"[NuTo::StructureGrid::SetAllElementIds] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//		std::cout<<"[NuTo::StructureGrid::SetAllElementIds] " << diffn.tv_sec <<" sec: "<<diffn.tv_nsec/1000000.<<" msec" << std::endl;
		std::cout<<"[NuTo::StructureGrid::SetAllElementIds] " << diffsum.tv_sec <<" sec: "<<diffsum.tv_nsec/1000000.<<" msec" << std::endl;
#endif
}

//! @brief Get LocalCoefficientMatrix0
//! @param NumLocalCoefficientMatrix0 number of stiffness matrix
NuTo::FullMatrix<double>* NuTo::StructureGrid::GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0)
{
 	if (rNumLocalCoefficientMatrix0<0 || rNumLocalCoefficientMatrix0>=GetNumMaterials())
        throw MechanicsException("[NuTo::StructureGrid::GetLocalCoefficientMatrix0] No valid material number.");
    return &mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0];
}
