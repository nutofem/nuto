// $Id $
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif
# ifdef _OPENMP
    #include <omp.h>
# endif


#include "nuto/mechanics/structures/grid/OctreeGrid.h"
#include "nuto/mechanics/structures/grid/MultiGridStructure.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/MortonOrder.h"
#include "nuto/mechanics/constitutive/mechanics/LinearElasticEngineeringStress.h"
#ifdef ENABLE_OPTIMIZE
#include "nuto/optimize/MisesWielandt.h"
#endif //ENABLE_OPTIMIZE
#include <algorithm>

//#define PLANESTRESS
#ifdef ENABLE_OPTIMIZE
NuTo::OctreeGrid::OctreeGrid(int rDimension):  CallbackHandlerGrid ()
#else //ENABLE_OPTIMIZE
NuTo::OctreeGrid::OctreeGrid(int rDimension)
#endif //ENABLE_OPTIMIZE
{
   if (rDimension!=3)
	{
		throw MechanicsException("[OctreeGrid::OctreeGrid] The dimension of the grid structure is so far 3.");
	}
    mDimension = rDimension;
    mNumVoxel=0;  //number of voxels
    mGridDimension.resize(rDimension);
    mVoxelSpacing.resize(rDimension);	//spacing between center of neighbor voxels / dimension of each voxel
    mGridOrigin.resize(rDimension);		// origin of the model , in the center of the first voxel
    mOctreeDimension.resize(rDimension);
    mMatrixFreeMethod=0;
    mUseDiagHessian =true;
   	mUseMisesWielandt =false;
    mImageDataFile="InputFile";
    mMortonOrderIsTrue=true;
    mNumMaterials=0;
    mNumBasisMaterials=1;
    mNumConstraintDofs=0;
    mC11=0,mC12=0,mC44=0;
    mCurrentGridNumber=0;
    mWeightingFactor=0.0;
    mpFineGrid=0;
    mpCoarseGrid=0;
    mNumLevels=UINT32_MAX;
}

NuTo::OctreeGrid::~OctreeGrid()
{	// coarse and fine grid pointer have not to be deleted, will be deleted in MultiGrid
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::OctreeGrid::Info()const
{
	std::cout<<"Octree grid Info\n "
			"-----------------------------------------------------------------------------------\n";
	std::cout<<"Dimension ...................................... " << mDimension << "\n";
	std::cout<<"GridNumber ..................................... " <<mCurrentGridNumber <<"\n";
	std::cout<<"GridDimension .................................. " <<mGridDimension[0]<<" "<<mGridDimension[1]<<" "<<mGridDimension[2]<<"\n";
	std::cout<<"Number of voxels ............................... "<<GetNumVoxels()<<"\n";
	std::cout<<"Number of elements ............................. "<<GetNumElements()<<"\n";
	std::cout<<"Number of dofs (free/constraint) ............... "<<GetNumNodes()*3<<" ("<<GetNumNodes()*3-GetNumConstraints()<<"/"<<GetNumConstraints()<<")\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
	std::cout<<"\n Solver configuration\n"
			"-----------------------------------------------------------------------------------\n";
	std::cout<<"MatrixFreeMethod is ............................ ";
	if (mMatrixFreeMethod) 								std::cout<<" NBN\n";
	else 												std::cout<<" EBE\n";
	std::cout<<"Preconditioner is .............................. ";
	if (mUseDiagHessian) 								std::cout<<" DIAGONAL D_ii\n";
	else 												std::cout<<" MULTIGRID\n";
	std::cout<<"Scaling of Preconditioner is ................... ";
	if (mUseMisesWielandt) 								std::cout<<" EIGENVALUE - 1/maxLambda\n";
	else 												std::cout<<" 1\n";
	if(GetWeightingFactor()!=0)
	std::cout<<"Value of scaling factor is ..................... "<<mWeightingFactor<<"\n";

//	else 												std::cout<<" 1/NumberParameters\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
	std::cout<<"Number of different materials .................. "<<GetNumMaterials()<<"\n";
//	for(int i=0;i<GetNumMaterials();++i)
//	{
		std::cout<<"Material "<<0<<"\n";
		std::cout<<"-   Youngs modulus .............................. "<<mData.find(0)->second.weight<<"\n";
//	}
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::OctreeGrid::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::OctreeGrid::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::OctreeGrid::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::OctreeGrid::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::OctreeGrid::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::OctreeGrid::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::OctreeGrid::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of grid structure" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP (mDimension)
       & BOOST_SERIALIZATION_NVP (mNumVoxel)
       & BOOST_SERIALIZATION_NVP (mGridDimension)
       & BOOST_SERIALIZATION_NVP (mVoxelSpacing)
       & BOOST_SERIALIZATION_NVP (mGridOrigin)
       & BOOST_SERIALIZATION_NVP (mOctreeDimension)
       & BOOST_SERIALIZATION_NVP (mImageDataFile)
       & BOOST_SERIALIZATION_NVP (mNumMaterials)
       & BOOST_SERIALIZATION_NVP (mNumBasisMaterials)
       & BOOST_SERIALIZATION_NVP (mNumConstraintDofs)
       & BOOST_SERIALIZATION_NVP (mLocalCoefficientMatrix0)
       & BOOST_SERIALIZATION_NVP (mBasisEdgeCoefficientMatrix0)
       & BOOST_SERIALIZATION_NVP (mLocalDerivativeShapeFunctions)
//	serialization does not work
//       & BOOST_SERIALIZATION_NVP (mDofIsConstraint)
       & BOOST_SERIALIZATION_NVP (mDisplacements)
       & BOOST_SERIALIZATION_NVP (mLinearElasticEngineeringStrains)
       & BOOST_SERIALIZATION_NVP (mLinearElasticEngineeringStresses)
       & BOOST_SERIALIZATION_NVP (mMatrixFreeMethod)
       & BOOST_SERIALIZATION_NVP(mUseDiagHessian)
       & BOOST_SERIALIZATION_NVP(mUseMisesWielandt)
       & BOOST_SERIALIZATION_NVP (mFineEdgeId)
       & BOOST_SERIALIZATION_NVP (mWeightingFactor)
       & BOOST_SERIALIZATION_NVP (mpFineGrid)
       & BOOST_SERIALIZATION_NVP (mpCoarseGrid)
       & BOOST_SERIALIZATION_NVP (mNumLevels);
#ifdef DEBUG_SERIALIZATION
     std::cout << "finish serialization of grid structure" << std::endl;
#endif
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::OctreeGrid::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::OctreeGrid::Save] Error opening file.");
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
            throw MechanicsException ( "[NuTo::OctreeGrid::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception &e )
    {
        std::string s ( std::string ( "[NuTo::OctreeGrid::Save]File save exception in boost - " ) + std::string ( e.what() ) );
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
void NuTo::OctreeGrid::Restore (const std::string &filename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::OctreeGrid::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::OctreeGrid::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::OctreeGrid::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
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
            throw MechanicsException ( "[NuTo::OctreeGrid::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception &e )
    {
        std::string s ( std::string ( "[NuTo::OctreeGrid::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
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
#endif

//! @brief import routine for basic grid data without OctreeGrid data space
void NuTo::OctreeGrid::ImportFromVtkASCIIFileHeader(std::string rFileName,size_t *rGridDimension,double *rVoxelSpacing,double *rGridOrigin, size_t rNumVoxel)
{
    std::cout<<__FILE__<<" "<<__LINE__<<" in ImportFromVtkASCIIFileHeader without OctreeGrid data \n ";
    try
    {
        using namespace boost::spirit::classic;
        // open file
        std::ifstream file(rFileName, std::ios::in);
        if (file.is_open() == false)
        {
            throw MechanicsException("[OctreeGrid::ImportFromVtkASCIIFile] error opening file.");
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
       if (parse(line.c_str(),("DIMENSIONS " >> uint_p[assign_a(rGridDimension[0])] >> ' '
                                           >> uint_p[assign_a(rGridDimension[1])] >> ' '
                                           >> uint_p[assign_a(rGridDimension[2])] >> *space_p)).full == false)
       {
            throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading dimension.");
       }
       // read spacing
        getline (file, line);
        if (parse(line.c_str(),("SPACING ">> real_p[assign_a(rVoxelSpacing[0])] >> ' '
                                >> real_p[assign_a(rVoxelSpacing[1])] >> ' '
                                >> real_p[assign_a(rVoxelSpacing[2])] >> *space_p)).full == false)
        {
            throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading spacing.");
        }
        // read origin
        getline (file, line);
        if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(rGridOrigin[0])] >> ' '
                                 >> real_p[assign_a(rGridOrigin[1])] >> ' '
                                 >> real_p[assign_a(rGridOrigin[2])] >> *space_p)).full == false)
        {
             throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading origin.");
        }
        // read number of entries
        getline (file, line);
        if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(rNumVoxel)] >>  *space_p)).full == false)
	   {
		   throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading number of entries.");
	   }
        // close file
       file.close();
       // test of saved data
      if (rNumVoxel<1)
       {
           throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader] error number of entries is negative or zero.");
       }
       if (rGridDimension[0]<1||rGridDimension[1]<1||rGridDimension[2]<1)
       {
           throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader] error dimension is negative or zero.");
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


//! @brief import routine for basic grid data with OctreeGrid data space
void NuTo::OctreeGrid::ImportFromVtkASCIIFileHeader(std::string rFileName)
{
    try
    {
        using namespace boost::spirit::classic;
        // open file
        std::ifstream file(rFileName, std::ios::in);
        if (file.is_open() == false)
        {
            throw MechanicsException("[OctreeGrid::ImportFromVtkASCIIFile] error opening file.");
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
            throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading dimension.");
       }
        // read spacing
        getline (file, line);
        if (parse(line.c_str(),("SPACING ">> real_p[assign_a(this->mVoxelSpacing[0])] >> ' '
                                >> real_p[assign_a(this->mVoxelSpacing[1])] >> ' '
                                >> real_p[assign_a(this->mVoxelSpacing[2])] >> *space_p)).full == false)
        {
            throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading spacing.");
        }
        // read origin
        getline (file, line);
        if (parse(line.c_str(),("ORIGIN ">> real_p[assign_a(this->mGridOrigin[0])] >> ' '
                                 >> real_p[assign_a(this->mGridOrigin[1])] >> ' '
                                 >> real_p[assign_a(this->mGridOrigin[2])] >> *space_p)).full == false)
        {
             throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading origin.");
        }
        // read number of entries
        getline (file, line);
        if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(this->mNumVoxel)] >>  *space_p)).full == false)
               {
                   throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader]error reading number of entries.");
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
           throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader] error number of entries is negative or zero.");
       }
       if (mGridDimension[0]<1||mGridDimension[1]<1||mGridDimension[2]<1)
       {
           throw MechanicsException("[OctreeGrid::importFromVtkASCIIFileReadHeader] error dimension is negative or zero.");
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

void NuTo::OctreeGrid::ImportFromVtkASCIIFile(const std::string rFileName,std::vector<int> &rData)
{
    using namespace boost::spirit::classic;

    // open file
    std::ifstream file(rFileName, std::ios::in);
    if (file.is_open() == false)
    {
        throw MechanicsException("[OctreeGrid::ImportFromVtkASCIIFile] error opening file.");
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
               throw MechanicsException("[OctreeGrid::importFromVtkASCIIFile]error reading number of entries.");
           }
    // read data type
    getline (file, line);
    // read empty line
    getline (file, line);

  // read entries
    if(numEntries!=rData.size())
        throw MechanicsException("[OctreeGrid::importFromVtkASCIIFile] number of entries is not equal to vector size.");

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
//! @brief create grid data in octree
//! @param rThresholdMaterialValue ... threshold between material one and two
//! @param imageValues ... vector of image data
//! @param rColorToMaterialData ... vector of material data (Young's Modulus) mapped to color points

void NuTo::OctreeGrid::CreateOctree(int rThresholdMaterialValue, std::string fileName,
		 std::vector<double>& rColorToMaterialData)
{
	bool fineResolutionBoundary=true;

	mImageDataFile=fileName;
	std::vector<int> imageValues (mNumVoxel);
	ImportFromVtkASCIIFile( fileName,imageValues);

	size_t numElements=0;

	uint32_t key=0;
	std::map<uint32_t,data> rData;

	for (size_t element=0;element<mNumVoxel;++element)
	{
		if(rColorToMaterialData[imageValues[element]]>0) // mData only for existing elements
		{
			uint32_t z=(element/(mGridDimension[0]*mGridDimension[1]))-1; // subtract frame element
			uint32_t residual1=element%((mGridDimension[0])*(mGridDimension[1]));
			uint32_t y=residual1/(mGridDimension[0])-1;
			uint32_t x=residual1%(mGridDimension[0])-1;
			key =MortonOrder::EncodeMorton3D(x,y,z);
//			std::cout<<" x "<<x<<" y "<<y<<" z "<<z<<" key "<<key<<"\n";
			data myData;
			myData.id=0;
			myData.level=0;
			myData.weight=rColorToMaterialData[imageValues[element]];
			rData[key]=myData;
			++numElements;

		}
	}
	// octree
	std::map<double,size_t> youngsModulus;
	size_t numMats=0;
	youngsModulus[rData.begin()->second.weight]=++numMats;
	// count number of materials and set in mNumMaterials
	for(std::map<uint32_t,data>::iterator it=rData.begin();it!=rData.end();++it)
	{
    	if(youngsModulus.find(it->second.weight)==youngsModulus.end()) // E does not yet exist
    		youngsModulus[it->second.weight]=++numMats;
	}
    mNumMaterials=numMats;

//	std::cout<<"Pow bitwise 8^2 "<<( 8<<3 )<<" 8^3 "<<(8<<6)<<" 8^4 "<< (8<<9)<<"\n";

	uint32_t neighbor=0,
			node=0;
	double weight=0.;
	uint32_t x=0,y=0,z=0,
			 level=0,
			 fac=0;
	uint32_t rNumLevels=1; //one level is minimum
	size_t numElementsOld=0;

	while(numElements!=numElementsOld && rNumLevels<mNumLevels)
    {
//    	std::cout<<"Coarsening of level "<<i<<"\n";
    	numElementsOld=numElements;
//    	std::cout<<"Start one level of coarsening\n ";
    	for(std::map<uint32_t,data>::iterator it=rData.begin();it!=rData.end();++it)
		{
    		// check material in neighbor nodes if key%8==0
			weight=it->second.weight;
			key=it->first;
			level=it->second.level,
			//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
			fac=0;
			//for testing only 3 levels
			if(it->second.level>2)
				break;
			if(key%(8<<(3*level))==0 && weight>0)
			{
//				std::cout<<"Element with key "<<key<<" is checked for coarsing.  \n";
				int count=1;
				bool equal=true;
				if(level==0)
					fac=1;
				else
					fac=(2<<(level-1));
				do
				{
					x=MortonOrder::DecodeMorton3X(count)*fac;
					y=MortonOrder::DecodeMorton3Y(count)*fac;
					z=MortonOrder::DecodeMorton3Z(count)*fac;
//					std::cout<<"Key "<<key<<" neighbor x "<<x<<" y "<<y<<" z "<<z<<"\n";
					neighbor=MortonOrder::Neighbor3D(key,x,y,z);

					if(rData.find(neighbor)!=rData.end())
					{
//						std::cout<<"Base E "<<weight<<" ? "<<rData.find(neighbor)->second.weight<<" \n";
						if((rData.find(neighbor)->second).level!=level ||(rData.find(neighbor)->second).weight!=weight)
						{
							equal=false;
							break;
						}
					}
					else //neighbor does not exist
					{
						equal=false;
						break;
					}

					// second check: coarse element neighbors have to be at most one level difference
					// reuse variables
					x=MortonOrder::DecodeMorton3X(count)*2*fac;
					y=MortonOrder::DecodeMorton3Y(count)*2*fac;
					z=MortonOrder::DecodeMorton3Z(count)*2*fac;
					neighbor=MortonOrder::Neighbor3D(key,x,y,z);
//					std::cout<<"Key "<<key<<" neighbor x "<<x<<" y "<<y<<" z "<<z<<"\n";
					uint32_t shift=count^7;
					if(rData.find(neighbor)!=rData.end())
					{
						//check also if neighbors of the coarse have the same base level 2:1 grid
						if((rData.find(neighbor)->second).level!=level)
						{
							equal=false;
//							std::cout<<"second check level of neighbor element returns false.\n";
							break;
						}

						if((shift&3)==3||(shift&5)==5||(shift&6)==6)
						{
							//plane
							x=MortonOrder::DecodeMorton3X(shift)*fac;
							y=MortonOrder::DecodeMorton3Y(shift)*fac;
							z=MortonOrder::DecodeMorton3Z(shift)*fac;
							// key of hanging node
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							if(rData.find(node)==rData.end())
							{
								equal=false;
//								std::cout<<"second check level of neighbor-neighbor element returns false.\n";
								break;
							}
						}
//						 shift in one direction with one
						// for all nodes on edges
						if((shift&1)==1)//shift&001==1 -> do shift x
						{
							x=MortonOrder::DecodeMorton3X(shift)*fac;
							y=0;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							if(rData.find(node)==rData.end())
							{
								equal=false;
//								std::cout<<"second check level of neighbor-neighbor element returns false.\n";
								break;
							}
						}
						if((shift&2)==2)
						{
							x=0;
							y=MortonOrder::DecodeMorton3Y(shift)*fac;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							if(rData.find(node)==rData.end())
							{
								equal=false;
//								std::cout<<"second check level of neighbor-neighbor element returns false.\n";
								break;
							}
						}
						if((shift&4)==4)//shift&100!=0 -> do shift z
						{
							x=0;
							y=0;
							z=MortonOrder::DecodeMorton3Z(shift)*fac;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							if(rData.find(node)==rData.end())
							{
								equal=false;
//								std::cout<<"second check level of neighbor-neighbor element returns false.\n";
								break;
							}

						}

						x=MortonOrder::DecodeMorton3X(count)*2*fac;
						y=MortonOrder::DecodeMorton3Y(count)*2*fac;
						z=MortonOrder::DecodeMorton3Z(count)*2*fac;
						neighbor=MortonOrder::Neighbor3D(key,x,y,z);
	//					std::cout<<"Key "<<key<<" neighbor x "<<x<<" y "<<y<<" z "<<z<<"\n";
						if(rData.find(neighbor)!=rData.end())
						{
							//check also if neighbors of the coarse have the same base level 2:1 grid
							if((rData.find(neighbor)->second).level!=level)
							{
								equal=false;
//								std::cout<<"second check level of neighbor element returns false.\n";
								break;
							}
						}
					}
					////////////////////////////////////////////////
					// added for fine resolution boundaries
					////////////////////////////////////////////////
					else
					{
						if(fineResolutionBoundary)
						{
							equal=false;
	//						std::cout<<"second check level of neighbor element returns false (no neighbor = fines scale).\n";
							break;
						}
					}

					// third check: neg. neighbors have to be at most one level difference if they exist
					x=MortonOrder::DecodeMorton3X(count)*2*fac;
					y=MortonOrder::DecodeMorton3Y(count)*2*fac;
					z=MortonOrder::DecodeMorton3Z(count)*2*fac;
//					std::cout<<"Key "<<key<<" neighbor x "<<x<<" y "<<y<<" z "<<z<<"\n";
					neighbor=MortonOrder::NegNeighbor3D(key,-1* (int32_t) x,-1* (int32_t) y,-1* (int32_t) z);
//					std::cout<<" neighbor "<<neighbor<<"\n";
					if(rData.find(neighbor)!=rData.end())
					{
//						std::cout<<"NegNeighbor "<<neighbor<<" of "<<key<<" exist. Levels "<<(rData.find(neighbor)->second).level<<" "<<level<<"\n";
						//check also if neighbors of the coarse have the same base level 2:1 grid
						if((rData.find(neighbor)->second).level<level)
						{
							equal=false;
//							std::cout<<"Third check level of neighbor element returns false.\n";
							break;
						}
					}
					else
					{
						// element does not exist, but maybe smaller
//						std::cout<<" big neg. neighbor does not exist.\n";
						x=MortonOrder::DecodeMorton3X(count)*fac;
						y=MortonOrder::DecodeMorton3Y(count)*fac;
						z=MortonOrder::DecodeMorton3Z(count)*fac;
//						std::cout<<"Key "<<key<<" neighbor x "<<x<<" y "<<y<<" z "<<z<<"\n";
						neighbor=MortonOrder::NegNeighbor3D(key,-1* (int32_t) x,-1* (int32_t) y,-1* (int32_t) z);
//						std::cout<<" neighbor "<<neighbor<<"\n";
						if(rData.find(neighbor)!=rData.end()) // find neg. neighbor of smaller level in same direction
						{
//							std::cout<<"New neighbor "<<neighbor<<" of "<<key<<" exist. Levels "<<(rData.find(neighbor)->second).level<<" "<<level<<"\n";
							//check also if neighbors of the coarse have the same base level 2:1 grid
							if((rData.find(neighbor)->second).level<level)
							{
								equal=false;
//									std::cout<<"Third check level of smaller negative neighbor element returns false.\n";
								break;
							}
						}
						////////////////////////////////////////////////
						// added for fine resolution boundaries
						////////////////////////////////////////////////
						else
						{
							if(fineResolutionBoundary)
							{
								equal=false;
//								std::cout<<"second check level of neighbor element returns false (no neighbor = fines scale).\n";
								break;
							}
						}
					}
					++count;
				}
				while (equal && count<7);

				if(equal)
				{
//					std::cout<<" Create coarse element.\n";
					it->second.level++;
					// weight of  coarse element
//					it->second.weight*=2;
					it->second.weight*=0.5;
					// delete fine elements 1-7 - set all zero
					for(int i=1;i<8;++i)
					{
						x=MortonOrder::DecodeMorton3X(i)*fac;
						y=MortonOrder::DecodeMorton3Y(i)*fac;
						z=MortonOrder::DecodeMorton3Z(i)*fac;
						neighbor=MortonOrder::Neighbor3D(key,x,y,z);
//						std::cout<<"Erase neighbor with key "<<neighbor<<"\n";
						rData.erase(neighbor);
					}

					numElements-=7;
//					for(int i=0;i<3;++i)
//						mOctreeDimension[i]--;
				}
			}
		}
		++rNumLevels;
	}

    mNumElements=numElements;
//    std::cout<<"Number of elements "<<mNumElements<<"\n";

    // create nodes
//	std::cout<<"Pow bitwise 2^2 "<<( 2<<1 )<<" 2^3 "<<(2<<2)<<" 2^4 "<< (2<<3)<<"\n";
	for(std::map<uint32_t,data>::iterator it=rData.begin();it!=rData.end();++it)
	{
		uint32_t key=it->first;
		uint32_t x=0,y=0,z=0,
				 level=it->second.level,
				 fac=0,
				 neighbor;
		//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
		if(level==0)
			fac=1;
		else
			fac=(2<<(level-1));
		if(it->second.weight>0)
		{
			for(uint32_t i=1;i<8;++i)
			{

				// compute the x,y,z shift for each node (local)
				x=MortonOrder::DecodeMorton3X(i)*(fac);
				y=MortonOrder::DecodeMorton3Y(i)*(fac);
				z=MortonOrder::DecodeMorton3Z(i)*(fac);
//				std::cout<<"Key "<<key<<" level "<<level<<"Shift -  x: "<<x<<" y: "<<y<<" z: "<<z<<" \n";
				neighbor=MortonOrder::Neighbor3D(key,x,y,z);

				if(rData.find(neighbor)==rData.end()) //node is not there
				{
					data myData;
	//					myData.key=neighbor;
					myData.id=0;
					myData.level=level;
					myData.weight=-1;
					rData[neighbor]=myData;
					// for  test purpose
					x=MortonOrder::DecodeMorton3X(key);
					y=MortonOrder::DecodeMorton3Y(key);
					z=MortonOrder::DecodeMorton3Z(key);
//					std::cout<<"Key "<<key<<" level "<<level<<" x: "<<x<<" y: "<<y<<" z: "<<z<<" \n";
					x=MortonOrder::DecodeMorton3X(neighbor);
					y=MortonOrder::DecodeMorton3Y(neighbor);
					z=MortonOrder::DecodeMorton3Z(neighbor);
//					std::cout<<"Add node "<<neighbor<<" x: "<<x<<" y: "<<y<<" z: "<<z<<" \n";
				}
			}
		}
	}

	// set ids for all data
	uint32_t countId=0;
	for(std::map<uint32_t,data>::iterator it=rData.begin();it!=rData.end();++it)
	{
		it->second.id=countId++;
	}
//	std::cout<<" nbr of ids "<<countId<<" numData "<<rData.size()<<"\n";
	assert(countId==rData.size());

	mData=rData;


	mNumBasisMaterials=1;

	mGridDimension[0]-=2; //subtract frame elements not needed for morton order
	mGridDimension[1]-=2; //subtract frame elements not needed for morton order
	mGridDimension[2]-=2; //subtract frame elements not needed for morton order

	// set hanging nodes with direction of constraint nodes
	HangingNodesSearch();

	mDisplacements.resize(3*(mData.size()+1),0.0);// initialized with zero
	mResidual.resize(3*(mData.size()+1),0.0);// initialized with zero
	mExtForces.resize(3*(mData.size()+1),0.0);// initialized with zero

}

//! @brief search for hanging nodes and saves key of nodes to which it is constraints
void NuTo::OctreeGrid::HangingNodesSearch()
{
	for(std::map<uint32_t,data>::iterator it=mData.begin();it!=mData.end();++it)
	{
		if(it->second.weight>0) // is an element
		{
			uint32_t neighbor=0,
					 node=0,
					 node2=0;
			uint32_t key=it->first;
			uint32_t x=0,y=0,z=0,
					 level=it->second.level,
					 fac=0;
			int32_t dx=0,dy=0,dz=0;
			if(level==0)
				fac=1;
			else
				fac=(2<<(level-1));
			uint32_t count=1;

//			std::cout<<"element "<<element++<<" key "<<key<<"\n";
			bool equal=true;
			do
			{
				// go through all nodes of element
				x=MortonOrder::DecodeMorton3X(count)*fac;
				y=MortonOrder::DecodeMorton3Y(count)*fac;
				z=MortonOrder::DecodeMorton3Z(count)*fac;
				neighbor=MortonOrder::Neighbor3D(key,x,y,z);
				std::map<uint32_t,data>::iterator it_neighbor=mData.find(neighbor);
//				std::cout<<" node nbr "<<count<<" key "<<it_neighbor->first<<"\n";

				if(it_neighbor->second.weight>0)//check if neighbor element exist
				{
					uint32_t fac_neighbor=0;
					if(it_neighbor->second.level==0)
						fac_neighbor=1;
					else
						fac_neighbor=(2<<(it_neighbor->second.level-1));

//					std::cout<<"node is element\n";
					//count^111 -> only one(s) in direction of shift
					uint32_t shift=count^7;
//					std::cout<<"shift "<<shift<<"\n";
					if(it_neighbor->second.level < level)//neighbor element level is smaller
					{
//						std::cout<<"neighbor element level is smaller (factor="<<fac_neighbor<<")\n";
						// set max three hanging nodes
						// shift two directions, 7 has no neighbor nodes to the it element
						//shift&011!=0 -> do shift xy
						if((shift&3)==3||(shift&5)==5||(shift&6)==6)
						{

							//three hanging nodes - here the one in a plane
							x=MortonOrder::DecodeMorton3X(shift)*fac_neighbor;
							y=MortonOrder::DecodeMorton3Y(shift)*fac_neighbor;
							z=MortonOrder::DecodeMorton3Z(shift)*fac_neighbor;
							// key of hanging node
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=shift; // takes ones of both in result
//							std::cout<<"hanging node key "<<node<<" const "<<shift <<"\n";
//							std::cout<<" x "<<MortonOrder::DecodeMorton3X(key)<<" "<<mGridDimension[0]<<" fac "<<fac<<"\n";
							if(MortonOrder::DecodeMorton3X(key)+fac==mGridDimension[0])
							{
								// double fac
								x=MortonOrder::DecodeMorton3X(shift)*fac;
								y=MortonOrder::DecodeMorton3Y(shift)*fac_neighbor;
								z=MortonOrder::DecodeMorton3Z(shift)*fac_neighbor;
								node=MortonOrder::Neighbor3D(neighbor,x,y,z);
//								std::cout<<"maxNode "<<node<<" count "<<count<<" shift "<<shift<<" Shift^1 "<<(shift^1)<<"\n";
								mData.find(node)->second.constraint|=(shift^1); // takes ones of both in result
//								std::cout<<"hanging node key "<<node<<" const "<<(shift^1)<<"\n";

							}
//							std::cout<<" y "<<MortonOrder::DecodeMorton3Y(key)<<" "<<mGridDimension[1]<<"\n";
							if(MortonOrder::DecodeMorton3Y(key)+fac==mGridDimension[1])
							{
								// double fac
								x=MortonOrder::DecodeMorton3X(shift)*fac_neighbor;
								y=MortonOrder::DecodeMorton3Y(shift)*fac;
								z=MortonOrder::DecodeMorton3Z(shift)*fac_neighbor;
								node=MortonOrder::Neighbor3D(neighbor,x,y,z);
//								std::cout<<"maxNode "<<node<<" count "<<count<<" shift "<<shift<<" Shift^2 "<<(shift^2)<<"\n";
								mData.find(node)->second.constraint|=(shift^2); // takes ones of both in result
//								std::cout<<"hanging node key "<<node<<" const "<<(shift^2)<<"\n";

							}
//							std::cout<<" z "<<MortonOrder::DecodeMorton3Z(key)<<" "<<mGridDimension[2]<<"\n";
							if(MortonOrder::DecodeMorton3Z(key)+fac==mGridDimension[2])
							{
								// double fac
								z=MortonOrder::DecodeMorton3Z(shift)*fac;
								x=MortonOrder::DecodeMorton3X(shift)*fac_neighbor;
								y=MortonOrder::DecodeMorton3Y(shift)*fac_neighbor;
								node=MortonOrder::Neighbor3D(neighbor,x,y,z);
//								std::cout<<"maxNode "<<node<<" count "<<count<<" shift "<<shift<<" Shift^4 "<<(shift^4)<<"\n";
								mData.find(node)->second.constraint|=(shift^4); // takes ones of both in result
//								std::cout<<"hanging node key "<<node<<" const "<<(shift^4)<<"\n";
							}


						}

						// shift in one direction with one
						// for all nodes on edges
						if((shift&1)!=0)//shift&001!=0 -> do shift x
						{
							x=MortonOrder::DecodeMorton3X(shift)*fac_neighbor;
							y=0;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=1;// save 001 - x-dir, takes ones of both in result
//							std::cout<<"hanging node key "<<node<<" const 1\n";
						}
						if((shift&2)!=0)//shift&010!=0 -> do shift y
						{
							x=0;
							y=MortonOrder::DecodeMorton3Y(shift)*fac_neighbor;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=2; // save 010 - y-dir
//							std::cout<<"hanging node key "<<node<<" const 2\n";
						}
						if((shift&4)!=0)//shift&100!=0 -> do shift z
						{
							x=0;
							y=0;
							z=MortonOrder::DecodeMorton3Z(shift)*fac_neighbor;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=4; // save 100 - z-dir
//							std::cout<<"hanging node key "<<node<<" const 4\n";
						}



					}
					//neighbor level is bigger
					//hanging nodes of it element, not of neighbor here
					else if (it_neighbor->second.level > level)
					{
//						std::cout<<"neighbor level is bigger \n";
						//fac stays the same
						//count^111 -> only one in direction of shift
						if((shift&3)==3||(shift&5)==5||(shift&6)==6)
						{
							//three hanging nodes - here the one in a plane
							x=MortonOrder::DecodeMorton3X(shift)*fac;
							y=MortonOrder::DecodeMorton3Y(shift)*fac;
							z=MortonOrder::DecodeMorton3Z(shift)*fac;
							// key of hanging node
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=shift; // takes ones of both in result
//							std::cout<<"hanging node key "<<node<<" const "<<shift <<"\n";

//							std::cout<<" x "<<MortonOrder::DecodeMorton3X(node)<<" "<<mGridDimension[0]<<" fac "<<fac<<"\n";
							if(MortonOrder::DecodeMorton3X(node)+fac==mGridDimension[0])
							{
								node2=MortonOrder::Neighbor3D(node,fac,0,0);
//								std::cout<<"shift "<<shift<<"\n";
								mData.find(node2)->second.constraint|=(shift^1);
//								std::cout<<"hanging node key "<<node2<<" const "<<(shift^1)<<"\n";
							}
//							std::cout<<" y "<<MortonOrder::DecodeMorton3Y(node)<<" "<<mGridDimension[1]<<"\n";
							if(MortonOrder::DecodeMorton3Y(node)+fac==mGridDimension[1])
							{
								node2=MortonOrder::Neighbor3D(node,0,fac,0);
//								std::cout<<"shift "<<shift<<"\n";
								mData.find(node2)->second.constraint|=(shift^2);
//								std::cout<<"hanging node key "<<node2<<" const "<<(shift^2)<<"\n";
							}
//							std::cout<<" z "<<MortonOrder::DecodeMorton3Z(node)<<" "<<mGridDimension[2]<<"\n";
							if(MortonOrder::DecodeMorton3Z(node)+fac==mGridDimension[2])
							{
								node2=MortonOrder::Neighbor3D(node,0,0,fac);
//								std::cout<<"shift "<<shift<<"\n";
								mData.find(node2)->second.constraint|=(shift^4);
//								std::cout<<"hanging node key "<<node2<<" const "<<(shift^4)<<"\n";
							}
						}

						// shift in one direction with one
						// for all nodes on edges
						if((shift&1)!=0)//shift&001!=0 -> do shift x
						{
							x=MortonOrder::DecodeMorton3X(shift)*fac;
							y=0;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=1;// save 001
//							std::cout<<"hanging node key "<<node<<" const 1\n";
						}
						if((shift&2)!=0)//shift&010!=0 -> do shift y
						{
							x=0;
							y=MortonOrder::DecodeMorton3Y(shift)*fac;
							z=0;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=2; // save 010
//							std::cout<<"hanging node key "<<node<<" const 2\n";

						}
						if((shift&4)!=0)//shift&100!=0 -> do shift z
						{
							x=0;
							y=0;
							z=MortonOrder::DecodeMorton3Z(shift)*fac;
							node=MortonOrder::Neighbor3D(neighbor,x,y,z);
							mData.find(node)->second.constraint|=4; // save 100
//							std::cout<<"hanging node key "<<node<<" const 4\n";
						}
					}
				}
				else //node is not element
				{
//					std::cout<<"node is not element\n";
//					std::cout<<"neighbor "<<neighbor<<"\n";
////					std::cout<<" x "<<MortonOrder::DecodeMorton3X(neighbor)<<" "<<mGridDimension[0]<<"\n";
//					std::cout<<" y "<<MortonOrder::DecodeMorton3Y(neighbor)<<" "<<mGridDimension[1]<<"\n";
//					std::cout<<" z "<<MortonOrder::DecodeMorton3Z(neighbor)<<" "<<mGridDimension[2]<<"\n";
					if(MortonOrder::DecodeMorton3X(neighbor)%2!=0)// x is uneven
					{
						dx=-fac;
					}
					else
					{
						dx=0;
					}
					if(MortonOrder::DecodeMorton3Y(neighbor)%2!=0)//y is uneven
					{
						dy=-fac;
					}
					else
					{
						dy=0;
					}
					if(MortonOrder::DecodeMorton3Z(neighbor)%2!=0)//z is uneven
						dz=-fac;
					else
						dz=0;

					node=MortonOrder::NegNeighbor3D(neighbor,dx,dy,dz);
					// constraint direction
					node2=MortonOrder::EncodeMorton3D(abs(dx),abs(dy),abs(dz));
					std::map<uint32_t,data>::iterator it_node=mData.find(node);
					if(it_node!=mData.end())
					{
//						std::cout<<" node exist\n";
						//node is element and level is bigger
						if(it_node->second.weight>0)
						{
							if(it_node->second.level>it_neighbor->second.level)
							{
								it_neighbor->second.constraint|=node2;
//								std::cout<<"hanging node key "<<neighbor<<" const "<<node2<<"\n";

							}
						}
						else //node exist, but is not element
						{
//							std::cout<<" node exist, but is not element\n";
							bool max_true=false;
							if(MortonOrder::DecodeMorton3X(neighbor)==mGridDimension[0])
							{
								dx-=2*fac;
								max_true=true;
							}
							if(MortonOrder::DecodeMorton3Y(neighbor)==mGridDimension[1])
							{
								dy-=2*fac;
								max_true=true;
							}
							if(MortonOrder::DecodeMorton3Z(neighbor)==mGridDimension[2])
							{
								dz-=2*fac;
								max_true=true;
							}
							if(max_true)
							{
								node=MortonOrder::NegNeighbor3D(neighbor,dx,dy,dz);
//								std::cout<<"help node "<<node<<" neighbor "<<neighbor;
//								std::cout<<" x "<<dx<<" y "<<dy<<" z "<<dz<<"\n";
								std::map<uint32_t,data>::iterator it_node=mData.find(node);
								if(it_node!=mData.end())
								{
									//node is element and level is bigger
									if(it_node->second.weight>0 && it_node->second.level>it_neighbor->second.level)
									{
										it_neighbor->second.constraint|=node2;
//										std::cout<<"hanging node key "<<neighbor<<" const "<<node2<<"\n";
									}
								}
							}
						}
					}
				}
			++count;
			}
			while (equal && count<7);
		}
	}
//	std::cout<<"const ";
//	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
//		std::cout<<(int) it->second.constraint<<" ";
//	std::cout<<"\n";
}

//! @brief correct solution for hanging nodes
//! @param displacement solution
 void NuTo::OctreeGrid::HangingNodesCorrection(std::vector<double>& u)const
{
	assert(u.size()-3==mData.size()*3||u.size()==mData.size()*3);
	// loop over "nodes" for check of constraint
	// do it on displacements
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{

		//check for hanging node and calculate by median of constraint nodes
		// for all six cases of constraint nodes
		// for all nodes not free
		uint32_t constraint=it->second.constraint;
		uint32_t id=it->second.id;
		uint32_t key=it->first;
		std::map<uint32_t,data>::const_iterator it_node;
		uint32_t neighbor,
				 level=it->second.level;
		uint32_t fac=0;
		//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
		if(level==0)
			fac=1;
		else
			fac=(2<<(level-1));

		int32_t xc=0,yc=0,zc=0;
		int32_t signA =-1;
		// for all three cases of constraint nodes on planes extra
		if((constraint&3)==3 || (constraint&5)==5 ||	(constraint&6)==6)
		{
//			std::cout<<"\nplane - key "<<key<<"(constraint= "<<constraint<<") : ";
			int32_t signB =-1;
			int32_t null=0;
			std::vector<int32_t *> signPlane(3);


			if((constraint&1)==1  ) // x
			{
				signPlane[0]=&signA;
				xc=fac;
				if ((constraint&2)==2 ) // y
				{
					signPlane[1]=&signB;
					signPlane[2]=&null;
					yc=fac;zc=0;
				}
				else //z
				{
					signPlane[1]=&null;
					signPlane[2]=&signB;
					yc=0;zc=fac;

				}
			}
			else
			{
				signPlane[0]=&null;
				signPlane[1]=&signA;
				signPlane[2]=&signB;
				xc=0;yc=fac;zc=fac;
			}
			if((constraint&8)!= 8) //   x-dir. is not constraint
				u[3*id  ]=0.;
			if((constraint&16)!=16) //  y-dir. is not constraint
				u[3*id+1]=0.;
			if((constraint&32)!=32) //  z-dir. is not constraint
				u[3*id+2]=0.;

			for(int helpA=0;helpA<2;++helpA) // two times loop
			{
				for(int helpB=0;helpB<2;++helpB) // two times loop
				{
					neighbor=MortonOrder::NegNeighbor3D(key,(*signPlane[0])*xc,*signPlane[1]*yc,*signPlane[2]*zc);
//					std::cout<<neighbor<<" ";
					it_node=mData.find(neighbor);
					if(it_node!=mData.end())
					{
						// overwrite local residual
						if((constraint&8)!= 8) //   x-dir. is not constraint
							u[3*id  ]+=u[3*(it_node->second.id)]/4.;
						if((constraint&16)!=16) //  y-dir. is not constuaint
							u[3*id+1]+=u[3*(it_node->second.id)+1]/4.;
						if((constraint&32)!=32) //  z-dir. is not constraint
							u[3*id+2]+=u[3*(it_node->second.id)+2]/4.;
					}
					//change sign
					signB*=-1;
				}
				//change sign
				signA*=-1;
			}
		}
		// for all three cases of constraint nodes on edges
		else if((constraint&1)==1 || (constraint&2)==2 || (constraint&4)==4) //  110 | 010 =  110 != 010 -> result only equal value when constraint equal value
		{
//			std::cout<<"\nedge - key "<<key<<"(constraint= "<<constraint<<") : ";
			if((constraint&1)==1)
			{
				xc=fac;yc=0;zc=0;
			}
			else if((constraint&2)==2)
			{
				xc=0;yc=fac;zc=0;
			}
			else if((constraint&4)==4)
			{
				xc=0;yc=0;zc=fac;
			}
			if((constraint&8)!= 8) //   x-dir. is not constraint
				u[3*id  ]=0.;
			if((constraint&16)!=16) //  y-dir. is not constraint
				u[3*id+1]=0.;
			if((constraint&32)!=32) //  z-dir. is not constraint
				u[3*id+2]=0.;


			for(int help=0;help<2;++help) // two times loop
			{
				// only one value differs from zero
				neighbor=MortonOrder::NegNeighbor3D(key,signA*xc,signA*yc,signA*zc);
//				std::cout<<neighbor<<" ";

				it_node=mData.find(neighbor);
				if(it_node!=mData.end())
				{
					if((constraint&8)!= 8) //   x-dir. is not constraint
						u[3*id  ]+=u[3*(it_node->second.id)  ]/2.;
					if((constraint&16)!=16) //  y-dir. is not constraint
						u[3*id+1]+=u[3*(it_node->second.id)+1]/2.;
					if((constraint&32)!=32) //  z-dir. is not constraint
						u[3*id+2]+=u[3*(it_node->second.id)+2]/2.;
				}
				//change sign
				signA*=-1;
			}

		}
	}

}
//! @brief returns number of Voxels
//! @return number of Voxels
size_t NuTo::OctreeGrid::GetNumVoxels() const
{
    	return mNumVoxel;
}

//! @brief returns  VoxelSpacing
//! @return VoxelSpacing
const std::vector<double> NuTo::OctreeGrid::GetVoxelSpacing() const
{
    return mVoxelSpacing;
}

//! @brief returns GridOrigin
 //! @return GridOrigin
const std::vector<double> NuTo::OctreeGrid::GetGridOrigin() const
{
     return mGridOrigin;
}

 //! @brief returns GridDimension
 //! @return GridDimension
const std::vector<size_t> NuTo::OctreeGrid::GetGridDimension() const
{
     return mGridDimension;
}
//! @brief returns GridDimension
 //! @return GridDimension
void  NuTo::OctreeGrid::SetGridDimension(std::vector<size_t> &rGridDimension)
{
	mGridDimension=rGridDimension;
}
//! @brief returns the number of nodes
//! @return number of nodes
size_t NuTo::OctreeGrid::GetNumNodes() const
{
	return mData.size();
}
//! @brief returns the number of elements
//! @return number of elements
int NuTo::OctreeGrid::GetNumElements() const
{
	return mNumElements;
}

//! @brief Get NumMaterials
//! @return NumMaterial
const int NuTo::OctreeGrid::GetNumMaterials() const
{
    return mNumMaterials;
}
//! @brief Get NumBasisMaterials
//! @return NumBasisMaterial
const int NuTo::OctreeGrid::GetNumBasisMaterials() const
{
    return mNumBasisMaterials;
}

//! @brief Set NumBasisMaterials
void NuTo::OctreeGrid::SetNumBasisMaterials(int rNumBasisMaterials)
{
	mNumBasisMaterials=rNumBasisMaterials;
}
//! @brief Get CurrentGridNumber
//! @return rCurrentGridNumber
const int NuTo::OctreeGrid::GetCurrentGridNumber() const
{
	return mCurrentGridNumber;
}

//! @brief Set CurrentGridNumber
void NuTo::OctreeGrid::SetCurrentGridNumber(int rCurrentGridNumber)
{
	mCurrentGridNumber=rCurrentGridNumber;
}

//! @brief Get number of Constraints
 //! @return NumConstraints
const size_t NuTo::OctreeGrid::GetNumConstraints() const
{
	 return mNumConstraintDofs;
}

void NuTo::OctreeGrid::SetMatrixFreeMethod(bool rMatrixFreeMethod)
{
	mMatrixFreeMethod=rMatrixFreeMethod;
	if(mVerboseLevel>5)
	{
		std::cout<<"MatrixFreeMethod is ";
		if (mMatrixFreeMethod)
			std::cout<<" NBN \n";
		else
			std::cout<<" EBE \n";
	}
}

bool NuTo::OctreeGrid::GetMatrixFreeMethod()
{
	return mMatrixFreeMethod;
}

//! @brief set basis element stiffness
//! @param rVoxelSpacing ... element length,
//! @param rPoissonsRatio
//! @return stiffnesMatrix ... basis element stiffness matrix
void NuTo::OctreeGrid::SetBasisElementStiffnessMatrix(double rPoissonsRatio,int rBasisMaterialNum)
{
	try
	{
		if (mDimension!=3)
		   throw MechanicsException("[OctreeGrid::SetBasisElementStiffnessMatrix] Only 3D is implemented.");
		if(rBasisMaterialNum>=mNumBasisMaterials)
			   throw MechanicsException("[OctreeGrid::SetBasisElementStiffnessMatrix] Basis Material number too large.");

		if(rBasisMaterialNum<(int) mLocalCoefficientMatrix0.size())
		   std::cout<<"[OctreeGrid::SetBasisElementStiffnessMatrix] This basis element stiffness matrix for material "<<rBasisMaterialNum<<" will be replaces.\n";

		mPoissonsRatio=rPoissonsRatio;
		NuTo::Structure myHelpStruc(3);

		myHelpStruc.SetVerboseLevel(0);
#ifdef SHOW_TIME
		myHelpStruc.SetShowTime(false);
#endif
		// create material law
		int myMat=myHelpStruc.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
		myHelpStruc.ConstitutiveLawSetPoissonsRatio(myMat, rPoissonsRatio);
		myHelpStruc.ConstitutiveLawSetYoungsModulus(myMat, 1.0);

		// create nodes
		NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(3);
		NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(8);
		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(0)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(1)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(3)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(2)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(4)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(5)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(7)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(6)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

		// first element create

	   // elementIncidence.Info();
		// elementIncidence has no influence on stiffness matrix
		// first element create
	    int myInterpolationType = myHelpStruc.InterpolationTypeCreate("Brick3D");
	    myHelpStruc.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
	    myHelpStruc.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

	    // elementIncidence.Info();
		int myHelpElement=myHelpStruc.ElementCreate(myInterpolationType, elementIncidence);
		//myHelpStruc.ElementTotalConvertToInterpolationType(1.e-6, 3);
		myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);
		int mySection1 = myHelpStruc.SectionCreate("VOLUME");
		myHelpStruc.ElementSetSection(myHelpElement,mySection1);

		// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
		NuTo::FullVector<int,Eigen::Dynamic> rows;
		NuTo::FullVector<int,Eigen::Dynamic> coluums;
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrix;
		myHelpStruc.ElementStiffness(0,stiffnessMatrix,rows,coluums );
//		if(mVerboseLevel>3)
//		myHelpStruc.BuildGlobalCoefficientMatrix0();
		// change node 2 <-> 3 and 6 <-> 7

		//  6 - 7
		//  |   |
		//  4 - 5
		//
		//  2 - 3
		//  |   |
		//  0 - 1

		//switch columns of node 2 and 3
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> helpMat=stiffnessMatrix.GetColumn(3*2);
		stiffnessMatrix.SetColumn(3*2,stiffnessMatrix.GetColumn(3*3));
		stiffnessMatrix.SetColumn(3*3,helpMat);
		helpMat=stiffnessMatrix.GetColumn(3*2+1);
		stiffnessMatrix.SetColumn(3*2+1,stiffnessMatrix.GetColumn(3*3+1));
		stiffnessMatrix.SetColumn(3*3+1,helpMat);
		helpMat=stiffnessMatrix.GetColumn(3*2+2);
		stiffnessMatrix.SetColumn(3*2+2,stiffnessMatrix.GetColumn(3*3+2));
		stiffnessMatrix.SetColumn(3*3+2,helpMat);
		//switch row of node 2 and 3
		helpMat=stiffnessMatrix.GetRow(3*2);
		stiffnessMatrix.SetRow(3*2,stiffnessMatrix.GetRow(3*3));
		stiffnessMatrix.SetRow(3*3,helpMat);
		helpMat=stiffnessMatrix.GetRow(3*2+1);
		stiffnessMatrix.SetRow(3*2+1,stiffnessMatrix.GetRow(3*3+1));
		stiffnessMatrix.SetRow(3*3+1,helpMat);
		helpMat=stiffnessMatrix.GetRow(3*2+2);
		stiffnessMatrix.SetRow(3*2+2,stiffnessMatrix.GetRow(3*3+2));
		stiffnessMatrix.SetRow(3*3+2,helpMat);

		//switch columns of node 6 and 7
		helpMat=stiffnessMatrix.GetColumn(3*6);
		stiffnessMatrix.SetColumn(3*6,stiffnessMatrix.GetColumn(3*7));
		stiffnessMatrix.SetColumn(3*7,helpMat);
		helpMat=stiffnessMatrix.GetColumn(3*6+1);
		stiffnessMatrix.SetColumn(3*6+1,stiffnessMatrix.GetColumn(3*7+1));
		stiffnessMatrix.SetColumn(3*7+1,helpMat);
		helpMat=stiffnessMatrix.GetColumn(3*6+2);
		stiffnessMatrix.SetColumn(3*6+2,stiffnessMatrix.GetColumn(3*7+2));
		stiffnessMatrix.SetColumn(3*7+2,helpMat);
		//switch row of node 6 and 7
		helpMat=stiffnessMatrix.GetRow(3*6);
		stiffnessMatrix.SetRow(3*6,stiffnessMatrix.GetRow(3*7));
		stiffnessMatrix.SetRow(3*7,helpMat);
		helpMat=stiffnessMatrix.GetRow(3*6+1);
		stiffnessMatrix.SetRow(3*6+1,stiffnessMatrix.GetRow(3*7+1));
		stiffnessMatrix.SetRow(3*7+1,helpMat);
		helpMat=stiffnessMatrix.GetRow(3*6+2);
		stiffnessMatrix.SetRow(3*6+2,stiffnessMatrix.GetRow(3*7+2));
		stiffnessMatrix.SetRow(3*7+2,helpMat);


		stiffnessMatrix.WriteToFile("stiffness","   ");
		if((int) mLocalCoefficientMatrix0.size()==rBasisMaterialNum)
		{
			std::vector<double> stiffness(24*24);
			for (int i=0;i<24;++i)
			{
				for (int j=0;j<24;++j)
				{
					stiffness[(24*i)+j]=stiffnessMatrix(i,j); //row based saved
				}
			}

			mLocalCoefficientMatrix0.push_back(stiffness);

			//calculate shape derivative functions
			const double rLocalCoordinates[3]={0,0,0};
			mLocalDerivativeShapeFunctions.resize(24);

			const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = myHelpStruc.ElementGetElementPtr(myHelpElement)->GetInterpolationType()->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(0);
			//std::cout << "derivativeShapeFunctionsGeometryNatural " << derivativeShapeFunctionsGeometryNatural.rows() << " " << derivativeShapeFunctionsGeometryNatural.cols() << std::endl;
			mLocalDerivativeShapeFunctions.resize(derivativeShapeFunctionsGeometryNatural.rows()*derivativeShapeFunctionsGeometryNatural.cols());
			for (int node=0; node<derivativeShapeFunctionsGeometryNatural.rows();node++)
			{
				for (int dim=0; dim<derivativeShapeFunctionsGeometryNatural.cols();dim++)
				{
					mLocalDerivativeShapeFunctions[node*derivativeShapeFunctionsGeometryNatural.cols()+dim] = derivativeShapeFunctionsGeometryNatural(node,dim);
				}
			}

			//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" mLocalDerivativeShapeFunctions "<<mLocalDerivativeShapeFunctions[0]<<" "<<mLocalDerivativeShapeFunctions[1]<<" "<<mLocalDerivativeShapeFunctions[2]<<"\n";
//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" mLocalDerivativeShapeFunctions "<<mLocalDerivativeShapeFunctions[21]<<" "<<mLocalDerivativeShapeFunctions[22]<<" "<<mLocalDerivativeShapeFunctions[23]<<"\n";
            LinearElasticEngineeringStress myMaterial;
			myMaterial.SetPoissonsRatio(rPoissonsRatio);
            myMaterial.SetYoungsModulus(1.0);
#ifdef PLANESTRESS
            myMaterial.CalculateCoefficients2DPlainStress(mC11,mC12,mC44);
			//static_cast<NuTo::LinearElasticEngineeringStress*> (myElementPointer->GetConstitutiveLaw(0))->CalculateCoefficients2DPlainStress(mC11,mC12,mC44);
#else
            myMaterial.CalculateCoefficients3D(mC11,mC12,mC44);
			//static_cast<NuTo::LinearElasticEngineeringStress*> (myElementPointer->GetConstitutiveLaw(0))->CalculateCoefficients3D(mC11,mC12,mC44);
#endif //PLANESTRESS
//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" stiffnessTensorCoefficients "<<mC11<<" "<<mC12<<" "<<" "<<mC44<<"\n";
		}
		else
			throw MechanicsException("[OctreeGrid::SetBasisElementStiffnessMatrix] Basis material number is wrong.");
		// for strain and stress computations
		// get shape function derivatives for element center
	}
	catch( ... )
	{
		throw MechanicsException("[OctreeGrid::SetBasisElementStiffnessMatrix] Error.");
	}


}

//! @brief create grid data
//! @param rThresholdMaterialValue ... threshold between material one and two
//! @param imageValues ... vector of image data
//! @param rColorToMaterialData ... vector of material data (Young's Modulus) mapped to color points

//! @brief set displacement boundary conditions
//! @brief all nodes at one plane perpendicular to direction are constrained
//! @param rDirection ... 0 = x,1 = y,2 = z
//! @param rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
//! @param rValue ... value of boundary displacement
//! @param rDisplVector ... output initial displacement vector
void NuTo::OctreeGrid::SetDisplacementConstraints(size_t rDirection,size_t *rGridLocation,double rValue,std::vector<double> &rDisplVector)
{
	//with Morton order
	//free node constraint=0
	//displacement constraint node constraint=8
	//constraint=1-7 hanging node with direction of constaint nodes

	uint32_t key=0;
//	if (mVerboseLevel>1)
//	{
//		std::cout<<"[OctreeGrid::SetDisplacementConstraints] constraint regions\n x: "<<rGridLocation[0]<<" - "<<rGridLocation[1];
//		std::cout<<" y: "<<rGridLocation[2]<<" - "<<rGridLocation[3];
//		std::cout<<" z: "<<rGridLocation[4]<<" - "<<rGridLocation[5]<<"\n";
//	}
	for(size_t dim2=rGridLocation[4];dim2<=rGridLocation[5];++dim2)
	{
		for(size_t dim1=rGridLocation[2];dim1<=rGridLocation[3];++dim1)
		{
			for(size_t dim0=rGridLocation[0];dim0<=rGridLocation[1];++dim0)
			{
				key= MortonOrder::EncodeMorton3D(dim0,dim1,dim2);
				std::map<uint32_t,data>::iterator it=mData.find(key);
				if (it!=mData.end()) // key exist
				{
					it->second.constraint|=(uint32_t)(8<<rDirection); //x - 8; y- 16, z - 32
//					std::cout<<" key "<<it->first<<" constraint "<<it->second.constraint<<"\n";
					mDisplacements[(it->second.id)*3+rDirection]=rValue;
					++mNumConstraintDofs;

				}
			}
		}
	}
	rDisplVector=mDisplacements;
}

//! @brief Approximate condition number through max and min diagonal singular values
//! @return approximate condition number
double NuTo::OctreeGrid::ApproximateSystemConditionNumber()
{
	size_t numNodes=mData.size();
	#pragma data present(mData,mLocalCoefficientMatrix0)
	#pragma acc parallel loop
	std::vector<double> rProduct(3*numNodes);
	std::vector<size_t> nodeNumbers(8);

	#pragma acc parallel loop
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		double elemStiff=it->second.weight;
		if(elemStiff>0)
		{
			elemStiff*=mLocalCoefficientMatrix0[0][0];
			// get nodes

			uint32_t key=it->first;
			uint32_t x=0,y=0,z=0,
					 level=it->second.level,
					 neighbor,
					 fac=0;
			//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
			if(level==0)
				fac=1;
			else
				fac=(2<<(level-1));

			rProduct[3*(it->second.id)]+=elemStiff;
			rProduct[3*(it->second.id)+1]+=elemStiff;
			rProduct[3*(it->second.id)+2]+=elemStiff;
			for(uint32_t node=1;node<8;++node)
			{
				x=MortonOrder::DecodeMorton3X(node)*fac;
				y=MortonOrder::DecodeMorton3Y(node)*fac;
				z=MortonOrder::DecodeMorton3Z(node)*fac;
				neighbor=MortonOrder::Neighbor3D(key,x,y,z);
				std::map<uint32_t,data>::const_iterator it_neighbor=mData.find(neighbor);
				rProduct[3*(it_neighbor->second.id)]+=elemStiff;
				rProduct[3*(it_neighbor->second.id)+1]+=elemStiff;
				rProduct[3*(it_neighbor->second.id)+2]+=elemStiff;
			}
		}
	}
	double maxValue=0.;
	double minValue=UINT32_MAX;

	size_t numDOFs=mData.size()*3;
	for (size_t i=0;i<numDOFs;++i)
	{
		if(maxValue<(rProduct[i]*rProduct[i]))
			maxValue=rProduct[i]*rProduct[i];
		if(minValue>(rProduct[i]*rProduct[i]))
			minValue=(rProduct[i]*rProduct[i]);
	}
	maxValue=sqrt(maxValue);
	minValue=sqrt(minValue);
	std::cout<<"[OctreeGrid] Max. scalar diagonal value of matrix  : "<<maxValue<<"\n";
	std::cout<<"[OctreeGrid] Min. scalar diagonal value of matrix  : "<<minValue<<"\n";
	std::cout<<"[OctreeGrid] Approximate condition number of matrix: "<<maxValue/minValue<<"\n";
	return (maxValue/minValue);
}

//! @brief Get LocalCoefficientMatrix0
//! @param NumLocalCoefficientMatrix0 number of stiffness matrix
const std::vector<double>* NuTo::OctreeGrid::GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0) const
{
 	if (rNumLocalCoefficientMatrix0<0 || rNumLocalCoefficientMatrix0>=GetNumBasisMaterials())
        throw MechanicsException("[NuTo::OctreeGrid::GetLocalCoefficientMatrix0] No valid material number.");
    return &mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0];
}
//! @brief calculate gradient
void NuTo::OctreeGrid::Gradient (std::vector<double>& rValues,std::vector<double>& rGradient)
{
		CalculateMatrixVectorProductEBE(rValues,rGradient);
}
//! @brief ... calculate matrix-vector-product in element-by-element way
//! @brief ... with local vectors
//! @param ... u - parameters input, r - gradient output
void NuTo::OctreeGrid::CalculateMatrixVectorProductEBE(std::vector<double>& u,std::vector<double>& r)const
{
	size_t numParas=r.size();
	int dofsElem=24;

	// make sure residual vector is reset with zero //
#pragma acc parallel loop present(r[0:numParas])
	for(size_t i=0;i<numParas;++i)
		r[i]=0;
	//hanging node correction
	HangingNodesCorrection(u);

	//----------------------------------------------//
	// global external force vector (active dofs)
    //	std::vector<double>  force(mNumParameters,1);

	//move local data in loop for acc ?
	#pragma acc firstprivate(residual[0,dofsElem],displacement[0,dofsElem])
	std::vector<double> residual (dofsElem);
	std::vector<double> displacement(dofsElem);
	uint32_t numNodes=8;
	std::vector<size_t> constraint(numNodes);
	std::vector<uint32_t> id(numNodes);
	//loop over all elements
	#pragma acc parallel loop present(u[0:numParas],r[0:numParas],mData[0:numElems],mLocalCoefficientMatrix0[0:576])
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		uint32_t key=it->first;
		uint32_t level=it->second.level;
		double youngsModulus=it->second.weight;

		if(youngsModulus>0) //element exist
		{
			//get nodes
			uint32_t x=0,y=0,z=0,
					nodeId=it->second.id,
					node,
					fac=0;
			//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
			if(level==0)
				fac=1;
			else
				fac=(2<<(level-1));
			//first node: id equal element
			id[0]=nodeId;
//			std::cout<<"element "<<node<<" weight "<<youngsModulus<<" nodes "<<node<<"\n";
//			key_nodes[0]=key;
			displacement[0]=u[3*nodeId];
			displacement[1]=u[3*nodeId+1];
			displacement[2]=u[3*nodeId+2];
			constraint[0]=it->second.constraint;

			std::map<uint32_t,data>::const_iterator it_node;

			#pragma acc parallel loop
//			std::cout<<"id : "<<id[0]<<" key "<<key<<" key nodes: \n";
			for(uint32_t count=1;count<numNodes;++count)
			{
				x=MortonOrder::DecodeMorton3X(count)*fac;
				y=MortonOrder::DecodeMorton3Y(count)*fac;
				z=MortonOrder::DecodeMorton3Z(count)*fac;
				node=MortonOrder::Neighbor3D(key,x,y,z);
				it_node=mData.find(node);
				id[count]=it_node->second.id;
//				std::cout<<it_node->first<<" ";
				constraint[count]=it_node->second.constraint;
//				key_nodes[count]=it_node->first;
				displacement[3*count+0]=u[3*id[count]+0];
				displacement[3*count+1]=u[3*id[count]+1];
				displacement[3*count+2]=u[3*id[count]+2];
			}

			#pragma acc parallel loop //needed or for each thread initialized already with zero
			for(int i=0;i<dofsElem;++i)
			{
				residual[i]=0;
				#pragma acc parallel loop
				for(int j=0;j<dofsElem;++j)
					residual[i]+=youngsModulus*mLocalCoefficientMatrix0[0][dofsElem*i+j]*displacement[j];

			}
			// if there is no displacement constraint write residual to r
			//the value saved in constraint is the direction in which the constraint nodes lay
			for(uint32_t count=0;count<numNodes;++count)
			{
				// not an hanging node in any direction
				if((constraint[count]&1)!= 1 && (constraint[count]&2)!= 2 && (constraint[count]&4)!= 4)
				{
					if((constraint[count]&8)!= 8) //   x-dir. is not constraint
						r[3*id[count]  ]+=residual[3*count  ];
					if((constraint[count]&16)!=16) //  y-dir. is not constraint
						r[3*id[count]+1]+=residual[3*count+1];
					if((constraint[count]&32)!=32) //  z-dir. is not constraint
						r[3*id[count]+2]+=residual[3*count+2];
				}
			}
		}
	}
}

void NuTo::OctreeGrid::BuildGlobalCoefficientMatrix(std::vector<double>& rKglob,std::vector<double>& rVector)const
{
	int dofsElem=24;

	std::vector<double> residual (dofsElem);
	uint32_t numNodes=8;
	std::vector<size_t> constraint(numNodes);
	std::vector<uint32_t> id(numNodes);
	int dofs=3*mData.size();//Nbr of dofs=3*nodes
	int freeDofs=0;
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		if((it->second.constraint&1)!= 1 && (it->second.constraint&2)!= 2 && (it->second.constraint&4)!= 4)
		{
			if((it->second.constraint&8)!= 8) //   x-dir. is not constraint
				++freeDofs;
			if((it->second.constraint&16)!=16) //  y-dir. is not constraint
				++freeDofs;
			if((it->second.constraint&32)!=32) //  z-dir. is not constraint
				++freeDofs;
		}
	}
	rKglob.resize(dofs);
	//loop over all elements
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		uint32_t key=it->first;
		uint32_t level=it->second.level;
		double youngsModulus=it->second.weight;

		if(youngsModulus>0) //element exist
		{
			//get nodes
			uint32_t x=0,y=0,z=0,
					nodeId=it->second.id,
					node,
					fac=0;
			//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
			if(level==0)
				fac=1;
			else
				fac=(2<<(level-1));
			//first node: id equal element
			id[0]=nodeId;
//			std::cout<<"element "<<node<<" weight "<<youngsModulus<<" nodes "<<node<<"\n";
//			key_nodes[0]=key;
			constraint[0]=it->second.constraint;

			std::map<uint32_t,data>::const_iterator it_node;

			#pragma acc parallel loop
//			std::cout<<"id : "<<id[0]<<" key "<<key<<" key nodes: \n";
			for(uint32_t count=1;count<numNodes;++count)
			{
				x=MortonOrder::DecodeMorton3X(count)*fac;
				y=MortonOrder::DecodeMorton3Y(count)*fac;
				z=MortonOrder::DecodeMorton3Z(count)*fac;
				node=MortonOrder::Neighbor3D(key,x,y,z);
				it_node=mData.find(node);
				id[count]=it_node->second.id;
//				std::cout<<it_node->first<<" ";
				constraint[count]=it_node->second.constraint;
//				key_nodes[count]=it_node->first;
			}

			#pragma acc parallel loop //needed or for each thread initialized already with zero
			for(uint32_t i=0;i<numNodes;++i)
			{
				// not an hanging node in any direction
//				if((constraint[i]&1)!= 1 && (constraint[i]&2)!= 2 && (constraint[i]&4)!= 4)
//				{
					#pragma acc parallel loop
					for(uint32_t j=0;j<numNodes;++j)
					{
//						if((constraint[i]&8)!= 8) //   x-dir. is not constraint
//						{
							for(int k=0;k<3;++k)
								rKglob[dofs*3*id[i]+3*id[j]+k]+=youngsModulus*mLocalCoefficientMatrix0[0][dofsElem*3*i+3*j+k];
//						}
//						if((constraint[i]&16)!=16) //  y-dir. is not constraint
//						{
							for(int k=0;k<3;++k)
								rKglob[dofs*3*id[i]+1+3*id[j]+k]+=youngsModulus*mLocalCoefficientMatrix0[0][dofsElem*3*i+1+3*j+k];
//						}
//						if((constraint[i]&32)!=32) //  z-dir. is not constraint
//						{
							for(int k=0;k<3;++k)
							{
								for(int k=0;k<3;++k)
									rKglob[dofs*3*id[i]+3+3*id[j]+k]+=youngsModulus*mLocalCoefficientMatrix0[0][dofsElem*3*i+2+3*j+k];
							}
//						}
					}
//				}
			}
		}
	}
//	// new loop to add constraint stiffness parts and delete constraint dofs
//	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
//	{
//		//check for hanging node and store stiffness part at free dofs
//		// for all six cases of constraint nodes
//		// for all nodes not free
//		uint32_t constraint0 = it->second.constraint;
//		uint32_t id=it->second.id;
//		std::map<uint32_t,data>::const_iterator it_node;
//		uint32_t neighbor;
//		uint32_t fac=0;
//		//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
//		if(level==0)
//			fac=1;
//		else
//			fac=(2<<(level-1));
//
//		int32_t xc=0,yc=0,zc=0;
//		int32_t signA =-1;
//		int32_t signB =-1;
//		int32_t null=0;
//		std::vector<int32_t *> signPlane(3);
//
//		// for all three cases of constraint nodes on planes extra
//		if((constraint0&3)==3 || (constraint0&5)==5 ||	(constraint0&6)==6)
//		{
////			std::cout<<"\nplane - key "<<key<<"(constraint0= "<<constraint0<<") : ";
//
//
//			if((constraint0&1)==1  ) // x
//			{
//				signPlane[0]=&signA;
//				xc=fac;
//				if ((constraint0&2)==2 ) // y
//				{
//					signPlane[1]=&signB;
//					signPlane[2]=&null;
//					yc=fac;zc=0;
//				}
//				else //z
//				{
//					signPlane[1]=&null;
//					signPlane[2]=&signB;
//					yc=0;zc=fac;
//
//				}
//			}
//			else
//			{
//				signPlane[0]=&null;
//				signPlane[1]=&signA;
//				signPlane[2]=&signB;
//				xc=0;yc=fac;zc=fac;
//			}
//			//hanging node can not have boundary conditions
//
//			for(int helpA=0;helpA<2;++helpA) // two times loop
//			{
//				for(int helpB=0;helpB<2;++helpB) // two times loop
//				{
//					neighbor=MortonOrder::NegNeighbor3D(key,(*signPlane[0])*xc,*signPlane[1]*yc,*signPlane[2]*zc);
////					std::cout<<neighbor<<" ";
//					it_node=mData.find(neighbor);
//					if(it_node!=mData.end())
//					{
//							rKglob[3*(it_node->second.id)*dofs+ 3*(it_node->second.id)]
//							       +=rKglob[dofsElem*3*(it_node->key%((it_node->second.level+1)*8))+2+3*j+k]
//							                                   [3*(it_node->second.id)]/4.;
//					}
//					//change sign
//					signB*=-1;
//				}
//				//change sign
//				signA*=-1;
//			}
//		}
//		// for all three cases of constraint0 nodes on edges
//		else if((constraint0&1)==1 || (constraint0&2)==2 || (constraint0&4)==4) //  110 | 010 =  110 != 010 -> result only equal value when constraint0 equal value
//		{
////			std::cout<<"\nedge - key "<<key<<"(constraint0= "<<constraint0<<") : ";
//			if((constraint0&1)==1)
//			{
//				xc=fac;yc=0;zc=0;
//			}
//			else if((constraint0&2)==2)
//			{
//				xc=0;yc=fac;zc=0;
//			}
//			else if((constraint0&4)==4)
//			{
//				xc=0;yc=0;zc=fac;
//			}
//			if((constraint0&8)!= 8) //   x-dir. is not constraint0
//				u[3*id  ]=0.;
//			if((constraint0&16)!=16) //  y-dir. is not constraint0
//				u[3*id+1]=0.;
//			if((constraint0&32)!=32) //  z-dir. is not constraint0
//				u[3*id+2]=0.;
//
//
//			for(int help=0;help<2;++help) // two times loop
//			{
//				// only one value differs from zero
//				neighbor=MortonOrder::NegNeighbor3D(key,signA*xc,signA*yc,signA*zc);
////				std::cout<<neighbor<<" ";
//
//				it_node=mData.find(neighbor);
//				if(it_node!=mData.end())
//				{
//					if((constraint0&8)!= 8) //   x-dir. is not constraint
//						u[3*id  ]+=u[3*(it_node->second.id)  ]/2.;
//					if((constraint0&16)!=16) //  y-dir. is not constraint
//						u[3*id+1]+=u[3*(it_node->second.id)+1]/2.;
//					if((constraint0&32)!=32) //  z-dir. is not constraint
//						u[3*id+2]+=u[3*(it_node->second.id)+2]/2.;
//				}
//				//change sign
//				signA*=-1;
//			}
//
//		}
//	}

	// test
	std::vector<int> myvector;

	// set some values (from 1 to 10)
	for (int i=1; i<=10; i++) myvector.push_back(i);

	// erase the 6th element
	myvector.erase (myvector.begin()+5);

	// erase the first 3 elements:
	myvector.erase (myvector.begin(),myvector.begin()+3);

	std::cout << "myvector contains:";
	for (unsigned i=0; i<myvector.size(); ++i)
	std::cout << ' ' << myvector[i];
	std::cout << '\n';


	  // new loop to set stiffness parts of constraint dofs =0
	// revers iterator
	for(std::map<uint32_t,data>::const_reverse_iterator it=mData.rbegin();it!=mData.rend();++it)
	{
		uint32_t constraint0 = it->second.constraint;
		uint32_t id=it->second.id;
		std::cout<<"node "<<id<<"\n";
		if((constraint0&32)== 32) //   z-dir. is constraint
		{
			rKglob.erase(rKglob.begin()+(3*id+2)*dofs,rKglob.begin()+(3*id+2)*dofs+dofs);
			for(int i=dofs-1;i>=0;--i)
				rKglob.erase(rKglob.begin()+i*dofs+3*id+2);

			rVector.erase(rVector.begin()+3*id+2);
			--dofs;
		}
		if((constraint0&16)== 16) //   y-dir. is constraint
		{
			rKglob.erase(rKglob.begin()+(3*id+1)*dofs,rKglob.begin()+(3*id+1)*dofs+dofs);
			for(int i=dofs-1;i>=0;--i)
				rKglob.erase(rKglob.begin()+i*dofs+3*id+1);

			rVector.erase(rVector.begin()+3*id+1);
			--dofs;
		}

		if((constraint0&8)== 8) //   x-dir. is constraint
		{
			rKglob.erase(rKglob.begin()+(3*id)*dofs,rKglob.begin()+(3*id)*dofs+dofs);
			for(int i=dofs-1;i>=0;--i)
				rKglob.erase(rKglob.begin()+i*dofs+3*id);

			rVector.erase(rVector.begin()+3*id);
			--dofs;
		}
	}
}

//! @brief get weighting factor for preconditioner
//! @return rWeight ... weighting factor
double  NuTo::OctreeGrid::GetWeightingFactor()const
{
	return mWeightingFactor;
}

//! @brief set weighting factor for preconditioner
//! @param rWeight ... weighting factor
void  NuTo::OctreeGrid::SetWeightingFactor(double rWeight)
{
	mWeightingFactor=rWeight;
}

void NuTo::OctreeGrid::Hessian(std::vector<double>&  rDiagHessian)
{
	if(mUseDiagHessian)
		HessianDiag(rDiagHessian);
}
void NuTo::OctreeGrid::HessianDiag(std::vector<double>& rHessianDiag)
{
	size_t numNodes=mData.size();
	if(mHessianDiag.empty()==true)
	{
	#pragma data present(rHessianDiag,mData,mLocalCoefficientMatrix0)
	#pragma acc parallel loop
		std::vector<double> rProduct(3*numNodes);

		std::vector<size_t> nodeNumbers(8);

		#pragma acc parallel loop
		for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
		{
			double elemStiff=it->second.weight;
			if(elemStiff>0)
			{
				elemStiff*=mLocalCoefficientMatrix0[0][0];
				// get nodes

				uint32_t key=it->first;
				uint32_t x=0,y=0,z=0,
						 level=it->second.level,
						 neighbor,
						 fac=0;
				//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
				if(level==0)
					fac=1;
				else
					fac=(2<<(level-1));

				rProduct[3*(it->second.id)]+=elemStiff;
				rProduct[3*(it->second.id)+1]+=elemStiff;
				rProduct[3*(it->second.id)+2]+=elemStiff;
				for(uint32_t node=1;node<8;++node)
				{
					x=MortonOrder::DecodeMorton3X(node)*fac;
					y=MortonOrder::DecodeMorton3Y(node)*fac;
					z=MortonOrder::DecodeMorton3Z(node)*fac;
					neighbor=MortonOrder::Neighbor3D(key,x,y,z);
					std::map<uint32_t,data>::const_iterator it_neighbor=mData.find(neighbor);
					// take into account the level of the element
					rProduct[3*(it_neighbor->second.id)]+=elemStiff;
					rProduct[3*(it_neighbor->second.id)+1]+=elemStiff;
					rProduct[3*(it_neighbor->second.id)+2]+=elemStiff;
				}
			}
		}
		size_t numDOFs=mData.size()*3;
		for (size_t i=0;i<numDOFs;++i)
		{
			if (rProduct[i]!=0)
			{
				rProduct[i]=1./rProduct[i];
			}
			else
				rProduct[i]=1.;
		}
		CalcScalingFactors(rProduct);
		mHessianDiag=rProduct;

		std::cout << " Hessian: ";
//		for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
//		{
//			if(it->second.constraint==0)
//				std::cout<<rProduct[3*(it->second.id)]<<" "<<rProduct[3*(it->second.id)+1]<<" "<<rProduct[3*(it->second.id)+2]<<" ";
//
//		}
//		std::cout<<"\n";
	}
	for (size_t para=0;para<numNodes*3;++para)
		rHessianDiag[para]*=mHessianDiag[para];
}
void NuTo::OctreeGrid::CalcScalingFactors(std::vector<double> &p)
{
	//diagonal scaling with scaling factor
	double scalefactor=GetWeightingFactor();
	std::cout<<"OctreeGrid: Weighting factor "<<scalefactor<<"\n";
	if(scalefactor==0.)
	{
#ifdef ENABLE_OPTIMIZE
		if (mUseMisesWielandt)
		{
			SetMisesWielandt(false); // if not, get a infinite loop
			NuTo::MisesWielandt myEigenCalculator(GetNumNodes()*3);
			NuTo::MisesWielandt::eObjectiveType type=myEigenCalculator.GetObjectiveType();
			myEigenCalculator.SetObjectiveType("MAX_EIGENVALUE_OF_PRECOND_MATRIX");
			// scale factor is 2/(2-lambda_max-lambda_min) [Meister: Num. lin. GLS] for Jacobi-Relaxation-Method
			myEigenCalculator.SetVerboseLevel(GetVerboseLevel());
			myEigenCalculator.SetCallback(this);

			myEigenCalculator.Optimize();
			double lambda_max=myEigenCalculator.GetObjective();
	//			double lambda_min=lambda_max/2.;
			// Jacobi-Relaxation-weighting
	//		scalefactor=2./(2-lambda_max-lambda_min);
			// my scale factor
	//		scalefactor=2./(lambda_max*lambda_max);
			scalefactor=1./lambda_max;
			// damping Jacobi: lampda of D-1 K Arbenz_2007
	//		scalefactor=4./(3.*lambda_max);
			SetMisesWielandt(true);
			myEigenCalculator.SetObjectiveType(type);

		}
		else
#endif //ENABLE_OPTIMIZE
			// 1 is needed for MisesWielandt (Hessian), make sure it is 1 in standard
			scalefactor=1.;
			//   scalefactor=1./(double) mNumParameters;
		SetWeightingFactor(scalefactor);
//		std::cout<<"[OctreeGrid] weighting factor "<<scalefactor<<"\n";

	}
	for (size_t count=0; count<GetNumNodes()*3; ++count)
        p[count] *=scalefactor;

}

void NuTo::OctreeGrid::GetEngineeringStrain(const std::vector<double> &rDisplacements, std::vector<double> &rEngineeringStrain)const
{
#ifdef PLANESTRESS
	std::cout<<"[OctreeGrid::GetEngineeringStrain] PLANESTRESS defined \n";
#endif

	//construct a help structure with one element
	NuTo::Structure myHelpStruc(3);

	myHelpStruc.SetVerboseLevel(0);
#ifdef SHOW_TIME
	myHelpStruc.SetShowTime(false);
#endif
	// create material law

	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(3);
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(8);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(0)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(1)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(2)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(3)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(4)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(5)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(6)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(7)=myHelpStruc.NodeCreateDOFs("displacements", nodeCoordinates);

	// first element create

   // elementIncidence.Info();
	int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
	int mySection1 = myHelpStruc.SectionCreate("VOLUME");
	myHelpStruc.ElementSetSection(myHelpElement,mySection1);
	Brick8N* myElementPointer=static_cast<Brick8N*> (myHelpStruc.ElementGetElementPtr(myHelpElement));

	size_t rNodesPerElement=8;
	std::vector<size_t> nodes(rNodesPerElement);
	rEngineeringStrain.resize(mNumElements*6);
	std::vector<double>nodeCoord(rNodesPerElement*3);
	double invJacobian[9], detJac;
	std::vector<double> derivativeShapeFunctionsGlobal(mLocalDerivativeShapeFunctions.size());

#ifdef PLANESTRESS
	for(size_t element=0;element<mNumElements;++element)
	{
		//CalculateNodalCoordinatesAtElement(element,nodes,nodeCoord);
		// get nodes

		myElementPointer->CalculateJacobian(mLocalDerivativeShapeFunctions,nodeCoord, invJacobian, detJac);

		myElementPointer->CalculateDerivativeShapeFunctionsGeometryGlobal(mLocalDerivativeShapeFunctions,invJacobian,
													derivativeShapeFunctionsGlobal);

		for(size_t i=0;i<rNodesPerElement;++i)
		{
			// e_ii=B u
			rEngineeringStrain[6*element+0]+=derivativeShapeFunctionsGlobal[3*i+0] *rDisplacements[3*nodes[i]+0];
			rEngineeringStrain[6*element+1]+=derivativeShapeFunctionsGlobal[3*i+1] *rDisplacements[3*nodes[i]+1];
			//e_xy
			rEngineeringStrain[6*element+5]+=derivativeShapeFunctionsGlobal[3*i+0] *rDisplacements[3*nodes[i]+1]
			                                +derivativeShapeFunctionsGlobal[3*i+1] *rDisplacements[3*nodes[i]+0];
		}
		rEngineeringStrain[6*element+2]=-mPoissonsRatio/(1-mPoissonsRatio)*(rEngineeringStrain[6*element+0]+rEngineeringStrain[6*element+1]);
			//e_yz
		rEngineeringStrain[6*element+3]=0.;
			//e_xz
		rEngineeringStrain[6*element+4]=0.;
	}
#else
	for(size_t element=0;element<mNumElements;++element)
	{
//		CalculateNodalCoordinatesAtElement(element,nodes,nodeCoord);

	    // deactivated temporarily
//		myElementPointer->CalculateJacobian(mLocalDerivativeShapeFunctions,nodeCoord, invJacobian, detJac);
	    double detJacUnusedWarningFixMe = detJac;

		myElementPointer->CalculateDerivativeShapeFunctionsGlobal(mLocalDerivativeShapeFunctions,invJacobian,
													derivativeShapeFunctionsGlobal);

		for(size_t i=0;i<rNodesPerElement;++i)
		{
			// e_ii=B u
			rEngineeringStrain[6*element+0]+=derivativeShapeFunctionsGlobal[3*i+0] *rDisplacements[3*nodes[i]+0];
			rEngineeringStrain[6*element+1]+=derivativeShapeFunctionsGlobal[3*i+1] *rDisplacements[3*nodes[i]+1];
			rEngineeringStrain[6*element+2]+=derivativeShapeFunctionsGlobal[3*i+2] *rDisplacements[3*nodes[i]+2];
			//e_yz
			rEngineeringStrain[6*element+3]+=derivativeShapeFunctionsGlobal[3*i+2] *rDisplacements[3*nodes[i]+1]
			                                +derivativeShapeFunctionsGlobal[3*i+1] *rDisplacements[3*nodes[i]+2];
			//e_xz
			rEngineeringStrain[6*element+4]+=derivativeShapeFunctionsGlobal[3*i+2] *rDisplacements[3*nodes[i]+0]
			                                +derivativeShapeFunctionsGlobal[3*i+0] *rDisplacements[3*nodes[i]+2];
			//e_xy
			rEngineeringStrain[6*element+5]+=derivativeShapeFunctionsGlobal[3*i+1] *rDisplacements[3*nodes[i]+0]
			                                +derivativeShapeFunctionsGlobal[3*i+0] *rDisplacements[3*nodes[i]+1];
		}
	}
#endif //PLANESTRESS
}
void NuTo::OctreeGrid::GetEngineeringStress( std::vector<double>& rEngineeringStrain, std::vector<double>& rEngineeringStress)const
{
#ifdef PLANESTRESS
	std::cout<<"[OctreeGrid::GetEngineeringStress] PLANESTRESS defined \n";
#endif
	//! @ToDo check strain caluclated
	size_t rNodesPerElement=8;
	std::vector<int> nodes(rNodesPerElement);
	rEngineeringStress.resize(mNumElements*6);
	double youngsModulus=0.;
#ifdef PLANESTRESS
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		youngsModulus=it->second.weight;
		if(youngsModulus>0)
		{
			uint32_t element=it->second.id;
			rEngineeringStress[6*element+0]=youngsModulus
											  *(mC11 * rEngineeringStrain[6*element+0]
											   +mC12 * rEngineeringStrain[6*element+1]);
			rEngineeringStress[6*element+1]=youngsModulus
										  *(mC12 * rEngineeringStrain[6*element+0]
										   +mC11 * rEngineeringStrain[6*element+1]);
			rEngineeringStress[6*element+2]=0.;
			//Sigma_yz
			rEngineeringStress[6*element+3]=0.;
			//Sigma_xz
			rEngineeringStress[6*element+4]=0.;
			//Sigma_xy
			rEngineeringStress[6*element+5]=youngsModulus*mC44*rEngineeringStrain[6*element+5];
		}
	}
#else
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		youngsModulus=it->second.weight;
		if(youngsModulus>0)
		{
			uint32_t element=it->second.id;
			rEngineeringStress[6*element+0]=youngsModulus
										  *(mC11 * rEngineeringStrain[6*element+0]
										   +mC12 * rEngineeringStrain[6*element+1]
										   +mC12 * rEngineeringStrain[6*element+2]);
			rEngineeringStress[6*element+1]=youngsModulus
										  *(mC12 * rEngineeringStrain[6*element+0]
										   +mC11 * rEngineeringStrain[6*element+1]
										   +mC12 * rEngineeringStrain[6*element+2]);
			rEngineeringStress[6*element+2]=youngsModulus
										  *(mC12 * rEngineeringStrain[6*element+0]
										   +mC12 * rEngineeringStrain[6*element+1]
										   +mC11 * rEngineeringStrain[6*element+2]);
			//Sigma_yz
			rEngineeringStress[6*element+3]=youngsModulus*mC44*rEngineeringStrain[6*element+3];
			//Sigma_xz
			rEngineeringStress[6*element+4]=youngsModulus*mC44*rEngineeringStrain[6*element+4];
			//Sigma_xy
			rEngineeringStress[6*element+5]=youngsModulus*mC44*rEngineeringStrain[6*element+5];
		}
	}
#endif //PLANESTRESS
}

// export to Vtk Datafile
void NuTo::OctreeGrid::ExportVTKUnstructuredGridDataFile(const std::string& rFilename) const
{
	std::ofstream file(rFilename.c_str());
    if (!file.is_open())
    	throw MechanicsException(std::string("[NuTo::OctreeGrid::ExportVTKUnstructuredGridDataFile] Error opening file ")+rFilename.c_str());

    // header /////////////////////////////////////////////////////////////////
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Data file was generated by NuTo" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID " << std::endl;
    ///////////////////////////////////////////////////////////////////////////
//    file << "DIMENSIONS "<<mGridDimension[0]-1<<" "<<mGridDimension[1]-1<<" "<<mGridDimension[2]-1<<"\n";
	size_t numNodes=mData.size();
    file << "POINTS "<<numNodes<<" floats\n";
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
		file<<MortonOrder::DecodeMorton3X(it->first)*mVoxelSpacing[0]<<" "<<MortonOrder::DecodeMorton3Y(it->first)*mVoxelSpacing[1]<<" "<<MortonOrder::DecodeMorton3Z(it->first)*mVoxelSpacing[2]<<"\n";

    file << "CELLS "<<mNumElements<<" "<<mNumElements*9<<"\n";//number of elements and number of total data in this block
    for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
    {
    	if(it->second.weight>0)
    	{
    		uint32_t key=it->first;
    		uint32_t level=it->second.level;
    		uint32_t fac=level-1;
			if(level==0)
				fac=1;
			else
				fac=(2<<fac);

    		file<<"8 "<<it->second.id<<" "<<mData.find(MortonOrder::Neighbor3D(key,fac,0,0))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,0,fac,0))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,fac,fac,0))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,0,0,fac))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,fac,0,fac))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,0,fac,fac))->second.id<<" "<<
    				mData.find(MortonOrder::Neighbor3D(key,fac,fac,fac))->second.id<<"\n";

    	}
    }
    file << "CELL_TYPES "<<mNumElements<<"\n";//number of elements
    for(size_t i=0;i<mNumElements;++i)
    	file <<"11 ";
    file<<"\n";
    file<<"CELL_DATA "<<mNumElements<<"\n";
    file<<"SCALARS mat float 1\n";
    file<<"LOOKUP_TABLE DEFAULT\n";
    for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
    {
   		if(it->second.weight>0)
   			file<<it->second.weight<<"\n";
   	}
	file<<"SCALARS elemId float 1\n";
	file<<"LOOKUP_TABLE DEFAULT\n";
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		if(it->second.weight>0)
			file<<it->second.id<<"\n";
	}

    file<<"POINT_DATA "<<numNodes<<"\n";

    file<<"SCALARS id int 1\n";
    file<<"LOOKUP_TABLE DEFAULT\n";
    for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
    {
   		file<<it->second.id<<"\n";
   	}
//    file<<"SCALARS key int 1\n";// set cell nbr
//    file<<"LOOKUP_TABLE DEFAULT\n";
//    for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
//    {
//   		file<<it->first<<"\n";
//   	}


	assert(3*(numNodes+1)==mDisplacements.size());
	file << "VECTORS displacements double \n";
	for (size_t i=0;i<numNodes;++i) // loop over grid points
		file<<mDisplacements[3*i]<<" " <<mDisplacements[3*i+1]<<" "<<mDisplacements[3*i+2]<<"\n";

//	file << "SCALARS constraints int \n";
//	file<<"LOOKUP_TABLE DEFAULT\n";
//	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
//		file<<(int) it->second.constraint<<"\n";
//
	file << "SCALARS hangingNodes int \n";
	file<<"LOOKUP_TABLE DEFAULT\n";
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		if((it->second.constraint&3)==3 || (it->second.constraint&5)==5 || (it->second.constraint&6)==6)
			file<<4<<" ";
		else if ((it->second.constraint&1)==1 || (it->second.constraint&2)==2 || (it->second.constraint&4)==4)
			file<<2<<" ";
		else
			file<<0<<" ";
	}
	file << "\n";

	file << "VECTORS constraints int \n";
//	file<<"LOOKUP_TABLE DEFAULT\n";
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		int help[3]={0,0,0};
		if((it->second.constraint&8)==8)
			help[0]=1;
		if((it->second.constraint&16)==16)
			help[1]=1;
		if((it->second.constraint&32)==32)
			help[2]=1;
		file<<help[0]<<" "<<help[1]<<" "<<help[2]<<"\n";
	}
//	file<<"\n";

   file.close();
//	std::vector<double> rStrainVector(0);
//	GetEngineeringStrain(mDisplacements,rStrainVector);
////	for(size_t i=0;i<rStrainVector.size();++i)
////	{
////		if(rStrainVector[i]!=0.)
////			std::cout<<" strain["<<i<<"]="<<rStrainVector[i]<<"\n";
////	}
//	std::vector<double> rStressVector(0);
//	GetEngineeringStress(rStrainVector,rStressVector);
//
////	for (size_t i=0;i<rStrainVector.size();++i)
////	{
////		if(abs(rStrainVector[i])<minDigit)
////			rStrainVector[i]=0.0;
////		if(abs(rStressVector[i])<minDigit)
////			rStressVector[i]=0.0;
////	}
//    file << "CELL_DATA "<<mNumVoxel<<"\n";
//	file << "TENSORS strains double \n";
//	size_t j=0;
//	for (size_t i=0;i<mNumVoxel;++i)
//	{
//		if(mVoxelId[j]==i)
//		{
//			file<<rStrainVector[6*j]<<" "<<rStrainVector[6*j+5]<<" "<<rStrainVector[6*j+4]<<"\n"
//			  <<rStrainVector[6*j+5]<<" "<<rStrainVector[6*j+1]<<" "<<rStrainVector[6*j+3]<<"\n"
//			  <<rStrainVector[6*j+4]<<" "<<rStrainVector[6*j+3]<<" "<<rStrainVector[6*j+2]<<"\n";
//			++j;
//		}
//		else
//			file<<"0.0 0.0 0.0 \n0.0 0.0 0.0 \n0.0 0.0 0.0 \n";
//	}
//
//	file << "TENSORS stress double \n";
//	j=0;
//	for (size_t i=0;i<mNumVoxel;++i)
//	{
//		if(mVoxelId[j]==i)
//		{
//			file<<rStressVector[6*j]<<" "<<rStressVector[6*j+5]<<" "<<rStressVector[6*j+4]<<"\n"
//			  <<rStressVector[6*j+5]<<" "<<rStressVector[6*j+1]<<" "<<rStressVector[6*j+3]<<"\n"
//			  <<rStressVector[6*j+4]<<" "<<rStressVector[6*j+3]<<" "<<rStressVector[6*j+2]<<"\n";
//			++j;
//		}
//		else
//			file<<"0.0 0.0 0.0 \n0.0 0.0 0.0 \n0.0 0.0 0.0 \n";
//	}
//	countNodes=0;


}

std::vector<double>&  NuTo::OctreeGrid::GetParameters()
{
	return mDisplacements;
}
std::vector<double>&  NuTo::OctreeGrid::GetResidual()
{
	return mResidual;
}
std::vector<double>&  NuTo::OctreeGrid::GetExtForces()
{
	return mExtForces;
}

void NuTo::OctreeGrid::SetParameters(std::vector<double>& rParameters)
{
	mDisplacements=rParameters;
}

void NuTo::OctreeGrid::SetResidual(std::vector<double>& rResidual)
{
	mResidual=rResidual;
}

void NuTo::OctreeGrid::SetExtForces(std::vector<double>& rExtForces)
{
	mExtForces=rExtForces;
}

//! @brief MultiGrid routines
//! @brief initialize on coarser grid
//! @param OctreeGrid reference
void NuTo::OctreeGrid::SetCoarserGridLevel(NuTo::OctreeGrid *rCoarseGrid)
{
	std::cout<<" todo \n";
}

//! @brief MultiGrid routines
//! @brief set pointer to coarser grid
void NuTo::OctreeGrid::SetCoarseGridPtr(NuTo::OctreeGrid* rCoarseGrid)
{
		mpCoarseGrid=rCoarseGrid;
}

//! @brief MultiGrid routines
//! @brief set pointer to finer grid
void NuTo::OctreeGrid::SetFineGridPtr(NuTo::OctreeGrid* rFineGrid)
{
		mpFineGrid=rFineGrid;
}



void NuTo::OctreeGrid::CalculateMultiGridCorrolations(std::string restrictionType,std::vector<double>& rRestriction,std::vector<double>& rProlongation)
{
	if (restrictionType=="TWENTYSEVENPOINTS")
	{
		// considered neighbor nodes:
		// 7pt: -dim0*dim1   -dim0  -1 0 -1 dim0  dim0*dim1
		// 27pt all neighbors
		const int neighbors=27;
		const std::vector<double> helpRestriction={1,2,1,2,4,2,1,2,1, 2,4,2,4,8,4,2,4,2, 1,2,1,2,4,2,1,2,1};

		rRestriction.resize(neighbors);
		rProlongation.resize(neighbors);
		// 1/8:  ..........  r coarse first step = r coarse calculated r_c=-K_c*u_c
		double unitFactorRestriction=1./8.;
		double unitFactorProlongation=1./8.;
		for (int i=0;i<neighbors;++i)
		{
			rRestriction[i]=helpRestriction[i]*unitFactorRestriction;
			rProlongation[i]=helpRestriction[i]*unitFactorProlongation;
		}
	}
	else
	{
		throw 	MechanicsException(std::string("[MultiGrid] Wrong restriction type."));
	}

}
void NuTo::OctreeGrid::Restriction(std::vector<double> &rRestrictionFactor)
{
	std::cout<<" todo \n";
}
void NuTo::OctreeGrid::Prolongation(std::vector<double> &rProlongationFactor)
{
	std::cout<<" todo \n";
}

void  NuTo::OctreeGrid::AnsysInput(std::vector<double> &rDisplVector) const
{
	// open file
	std::ofstream file;
    file.open("ansysInput");
    file<<"fini \n/clear,nostart \n/prep7 \n";
    file<<"!Ansys Input File: \n !DIM: "<<mGridDimension[0]<<" DOFS: "<<mData.size()*3<<"\n";
    file<<"et,1,solid185 \nkeyopt,1,2,3 \n";

    assert(mNumMaterials);
    int mat=1;
    std::map<int,double> E;
    double weight=0;
    bool newMat=true;
    std::map<uint32_t,data>::const_iterator it=mData.begin();
	E[mat]=it->second.weight;
	file<<"mp,ex,"<<mat<<","<<it->second.weight<<" \n";
	file<<"mp,prxy,"<<mat<<",.2 \n";
	++mat;

 	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
    {
 		weight=it->second.weight;
 		if(weight>0)
     	{
			for(std::map<int,double>::iterator it_e=E.begin(); it_e==E.end();++it_e)
			{
				if(it_e->second==weight)
				{
					newMat=false;
					it_e=E.end();
				}
			}
			if(newMat) //true
			{
				E[mat]=weight;
				file<<"mp,ex,"<<mat<<","<<weight<<" \n";
				file<<"mp,prxy,"<<mat<<",.2 \n";
				++mat;
				newMat=false;
			}

      	}
    }

    // create nodes
 	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
 		file<<"n,"<<(it->second.id)+1<<","<<MortonOrder::DecodeMorton3X(it->first)*mVoxelSpacing[0]<<","<<MortonOrder::DecodeMorton3Y(it->first)*mVoxelSpacing[1]<<"," 		<<MortonOrder::DecodeMorton3Z(it->first)*mVoxelSpacing[2]<<"\n";


    // create elements
    mat=1;
    file<<"type,1 \n mat,"<<mat<<" \n";
   std::vector<size_t> nodes(8);
    int elem=1;

	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
     	weight=it->second.weight;
		if(weight>0)
     	{
     		uint32_t key=it->first;
     		uint32_t level=it->second.level;
     		uint32_t fac=level-1;
 			if(level==0)
 				fac=1;
 			else
 				fac=(2<<fac);

 			// set mat
			for(std::map<int,double>::iterator it_e=E.begin(); it_e==E.end();++it_e)
			{
				if(it_e->second==weight)
				{
					if(it_e->first!=mat)
						file<<"mat,"<<it_e->first<<"\n";
					it_e==E.end();
				}
			}
     		file<<"en,"<<elem<<","<<(it->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,fac,0,0))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,fac,fac,0))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,0,fac,0))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,0,0,fac))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,fac,0,fac))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,fac,fac,fac))->second.id)+1<<","<<
     				(mData.find(MortonOrder::Neighbor3D(key,0,fac,fac))->second.id)+1<<"\n";
     		++elem;
     	}

    }

 	// create displacement constraints
	// by hand
	///////////////////////////////////////////////////////
	// x=0;ux=0
//	file<<"nsel,s,node,loc,x,0\n";
//	file<<"d,all,ux,0.\n";
//
//	// y=0;uy=0
//	file<<"nsel,s,node,loc,y,0\n";
//	file<<"d,all,uy,0.\n";
//
//	// z=0;uz=0
//	file<<"nsel,s,node,loc,z,0\n";
//	file<<"d,all,uz,0.\n";
//
//	// y=max;uy=1
//
//	file<<"alls \n"
//		"*GET,umax,node,,mxloc,y\n";
//	file<<"nsel,s,node,loc,y,umax\n";
//	file<<"d,all,uy,1.\n";
//	file<<"alls\n";
// 	// create displacement constraints
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{
		if((it->second.constraint&8)==8)
		{
			file<<"nsel,s,node,,"<<(it->second.id)+1<<"\n";
			file<<"d,all,ux,"<<rDisplVector[3*(it->second.id)]<<"\n";
		}
		if((it->second.constraint&16)==16)
		{
			file<<"nsel,s,node,,"<<(it->second.id)+1<<"\n";
			file<<"d,all,uy,"<<rDisplVector[3*(it->second.id)+1]<<"\n";
		}
		if((it->second.constraint&32)==32)
		{
			file<<"nsel,s,node,,"<<(it->second.id)+1<<"\n";
			file<<"d,all,uz,"<<rDisplVector[3*(it->second.id)+2]<<"\n";
		}
	}
    file<<"alls\n";

    //create linear constraint equation
    // no node is displacements constraint!!!
    int count=0; //number of constraint equations
	for(std::map<uint32_t,data>::const_iterator it=mData.begin();it!=mData.end();++it)
	{

		//check for hanging node
		// for all six cases of constraint nodes
		// for all nodes not free
		uint32_t constraint=it->second.constraint;
		uint32_t id=it->second.id;
		uint32_t key=it->first;
		std::map<uint32_t,data>::const_iterator it_node;
		uint32_t neighbor,
				 level=it->second.level;
		uint32_t fac=0;
		//level 1: 2^1=2 ~ <<0,level 2: 2^2=4 ~ <<1, ...
		if(level==0)
			fac=1;
		else
			fac=(2<<(level-1));

		int32_t xc=0,yc=0,zc=0;
		int32_t signA =-1;
		// for all three cases of constraint nodes on planes extra
		if((constraint&3)==3 || (constraint&5)==5 ||	(constraint&6)==6)
		{
//			std::cout<<"\nplane - key "<<key<<"(constraint= "<<constraint<<") : ";
			int32_t signB =-1;
			int32_t null=0;
			std::vector<int32_t *> signPlane(3);


			if((constraint&1)==1  ) // x
			{
				signPlane[0]=&signA;
				xc=fac;
				if ((constraint&2)==2 ) // y
				{
					signPlane[1]=&signB;
					signPlane[2]=&null;
					yc=fac;zc=0;
				}
				else //z
				{
					signPlane[1]=&null;
					signPlane[2]=&signB;
					yc=0;zc=fac;

				}
			}
			else
			{
				signPlane[0]=&null;
				signPlane[1]=&signA;
				signPlane[2]=&signB;
				xc=0;yc=fac;zc=fac;
			}
			uint32_t planeNodes[4];
			for(int helpA=0;helpA<2;++helpA) // two times loop
			{
				for(int helpB=0;helpB<2;++helpB) // two times loop
				{
					neighbor=MortonOrder::NegNeighbor3D(key,(*signPlane[0])*xc,*signPlane[1]*yc,*signPlane[2]*zc);
//					std::cout<<neighbor<<" ";
					it_node=mData.find(neighbor); // no check if node exist, have to
					planeNodes[2*helpA+helpB]=it_node->second.id;
					//change sign
					signB*=-1;
				}
				//change sign
				signA*=-1;
			}
			// CE,lhs,0,HNid,dof,factor,node,dof,0.25, ....
			file<<"CE,"<<++count<<",0,"<<id+1<<",ux,1,"<<planeNodes[0]+1<<",ux,-0.25,"
						<<planeNodes[1]+1<<",ux,-0.25 \n";
			file<<"CE,"<<count<<",0,"<<id+1<<",ux,1,"<<planeNodes[2]+1<<",ux,-0.25,"
						<<planeNodes[3]+1<<",ux,-0.25 \n";

			file<<"CE,"<<++count<<",0,"<<id+1<<",uy,1,"<<planeNodes[0]+1<<",uy,-0.25,"
						<<planeNodes[1]+1<<",uy,-0.25 \n";
			file<<"CE,"<<count<<",0,"<<id+1<<",uy,1,"<<planeNodes[2]+1<<",uy,-0.25,"
						<<planeNodes[3]+1<<",uy,-0.25 \n";

			file<<"CE,"<<++count<<",0,"<<id+1<<",uz,1,"<<planeNodes[0]+1<<",uz,-0.25,"
						<<planeNodes[1]+1<<",uz,-0.25 \n";
			file<<"CE,"<<count<<",0,"<<id+1<<",uz,1,"<<planeNodes[2]+1<<",uz,-0.25,"
						<<planeNodes[3]+1<<",uz,-0.25 \n";
		}
		// for all three cases of constraint nodes on edges
		else if((constraint&1)==1 || (constraint&2)==2 || (constraint&4)==4) //  110 | 010 =  110 != 010 -> result only equal value when constraint equal value
		{
//			std::cout<<"\nedge - key "<<key<<"(constraint= "<<constraint<<") : ";
			if((constraint&1)==1)
			{
				xc=fac;yc=0;zc=0;
			}
			else if((constraint&2)==2)
			{
				xc=0;yc=fac;zc=0;
			}
			else if((constraint&4)==4)
			{
				xc=0;yc=0;zc=fac;
			}

			uint32_t edgeNode[2];
			for(int help=0;help<2;++help) // two times loop
			{
				// only one value differs from zero
				neighbor=MortonOrder::NegNeighbor3D(key,signA*xc,signA*yc,signA*zc);
//				std::cout<<neighbor<<" ";
				it_node=mData.find(neighbor);
				edgeNode[help]=it_node->second.id;
				//change sign
				signA*=-1;
			}
			// CE,lhs,0,HNid,dof,factor,node,dof,0.25, ....
			file<<"CE,"<<++count<<",0,"<<id+1<<",ux,1,"<<edgeNode[0]+1<<",ux,-0.5,"
						<<edgeNode[1]+1<<",ux,-0.5 \n";
			file<<"CE,"<<++count<<",0,"<<id+1<<",uy,1,"<<edgeNode[0]+1<<",uy,-0.5,"
						<<edgeNode[1]+1<<",uy,-0.5 \n";
			file<<"CE,"<<++count<<",0,"<<id+1<<",uz,1,"<<edgeNode[0]+1<<",uz,-0.5,"
						<<edgeNode[1]+1<<",uz,-0.5 \n";

		}
	}

	// solution
	file<<"alls\n";
	file<<"!nummeriert Knoten neu durch\n"
			"!numcmp,node \n"
			"/solu \n"
			"eqslv,sparse	! solve with sparse direct sover \n"
			"! to get .full file \n"
			"outr,all,all	! alle Zwischenergebnisse speichern  \n"
			"wrfull,on		! stops analysis after writing file \n"
			"solve	           	! Lösen des Gleichungssystems \n"
			"finish	      ! Verlassen des Lösungsprozessors \n"
			"/AUX2 \n"
			" FILE,,full,		!DEFINE THE FILE \n"
			"HBMAT,,,,ASCII,stiff,yes,		!DUMP MATRIX FILE IN HARWELL-BOEING FORMAT \n"
			"FINISH \n";




//	/post1
//	!only for direct solver
//	!*DMAT,mat_gl,D,import,FULL,file.full,STIFF
//	!export,mat_gl,MMF,mat_gl.mmf,


//	!plnsol,u,x
    file.close();
}
