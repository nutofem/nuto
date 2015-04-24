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
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif
# ifdef _OPENMP
    #include <omp.h>
# endif


#include "nuto/mechanics/structures/grid/StructureGrid.h"
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
NuTo::StructureGrid::StructureGrid(int rDimension):  CallbackHandlerGrid ()
#else
NuTo::StructureGrid::StructureGrid(int rDimension)
#endif
{
   if (rDimension!=3)
	{
		throw MechanicsException("[StructureGrid::StructureGrid] The dimension of the grid structure is either so far 3.");
	}
    mDimension = rDimension;
    mNumVoxel=0;  //number of voxels
    mGridDimension.resize(rDimension);
    mVoxelSpacing.resize(rDimension);	//spacing between center of neighbor voxels / dimension of each voxel
    mGridOrigin.resize(rDimension);		// origin of the model , in the center of the first voxel
    mMatrixFreeMethod=0;
    mUseDiagHessian =true;
   	mUseMisesWielandt =true;
    mImageDataFile="InputFile";
    mNumMaterials=0;
    mNumBasisMaterials=1;
    mNumConstraintDofs=0;
    mC11=0,mC12=0,mC44=0;
    mCurrentGridNumber=0;
    mWeightingFactor=0.;
    mpFineGrid=0;
    mpCoarseGrid=0;
}

NuTo::StructureGrid::~StructureGrid()
{	// coarse and fine grid pointer have not to be deleted, will be deleted in MultiGrid
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureGrid::Info()const
{
	std::cout<<"Structured grid Info\n "
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
		std::cout<<"-   Youngs modulus .............................. "<<mYoungsModulus[0]<<"\n";
//	}
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
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of grid structure" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP (mDimension)
       & BOOST_SERIALIZATION_NVP (mNumVoxel)
       & BOOST_SERIALIZATION_NVP (mGridDimension)
       & BOOST_SERIALIZATION_NVP (mVoxelSpacing)
       & BOOST_SERIALIZATION_NVP (mGridOrigin)
       & BOOST_SERIALIZATION_NVP (mImageDataFile)
       & BOOST_SERIALIZATION_NVP (mNumMaterials)
       & BOOST_SERIALIZATION_NVP (mNumBasisMaterials)
       & BOOST_SERIALIZATION_NVP (mNumConstraintDofs)
       & BOOST_SERIALIZATION_NVP (mLocalCoefficientMatrix0)
       & BOOST_SERIALIZATION_NVP (mBasisEdgeCoefficientMatrix0)
       & BOOST_SERIALIZATION_NVP (mLocalDerivativeShapeFunctions)
       & BOOST_SERIALIZATION_NVP (mNeighborNodesNE)
       & BOOST_SERIALIZATION_NVP (mEdgeId)
       & BOOST_SERIALIZATION_NVP (mNodeId)
       & BOOST_SERIALIZATION_NVP (mVoxelId)
//	serialization does not work
//       & BOOST_SERIALIZATION_NVP (mDofIsConstraint)
       & BOOST_SERIALIZATION_NVP (mDisplacements)
       & BOOST_SERIALIZATION_NVP (mLinearElasticEngineeringStrains)
       & BOOST_SERIALIZATION_NVP (mLinearElasticEngineeringStresses)
       & BOOST_SERIALIZATION_NVP (mYoungsModulus)
       & BOOST_SERIALIZATION_NVP (mMatrixFreeMethod)
       & BOOST_SERIALIZATION_NVP(mUseDiagHessian)
       & BOOST_SERIALIZATION_NVP(mUseMisesWielandt)
       & BOOST_SERIALIZATION_NVP (mFineEdgeId)
       & BOOST_SERIALIZATION_NVP (mWeightingFactor)
       & BOOST_SERIALIZATION_NVP (mpFineGrid)
       & BOOST_SERIALIZATION_NVP (mpCoarseGrid);
#ifdef DEBUG_SERIALIZATION
     std::cout << "finish serialization of grid structure" << std::endl;
#endif
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
    catch ( boost::archive::archive_exception &e )
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
    catch ( boost::archive::archive_exception &e )
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
void NuTo::StructureGrid::ImportFromVtkASCIIFileHeader(std::string rFileName,size_t *rGridDimension,double *rVoxelSpacing,double *rGridOrigin, size_t rNumVoxel)
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
       if (parse(line.c_str(),("DIMENSIONS " >> uint_p[assign_a(rGridDimension[0])] >> ' '
                                           >> uint_p[assign_a(rGridDimension[1])] >> ' '
                                           >> uint_p[assign_a(rGridDimension[2])] >> *space_p)).full == false)
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
        if (parse(line.c_str(),("POINT_DATA ">> uint_p[assign_a(rNumVoxel)] >>  *space_p)).full == false)
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
void NuTo::StructureGrid::ImportFromVtkASCIIFileHeader(std::string rFileName)
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

void NuTo::StructureGrid::ImportFromVtkASCIIFile(const std::string rFileName,std::vector<int> &rData)
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
const std::vector<double> NuTo::StructureGrid::GetVoxelSpacing() const
{
    return mVoxelSpacing;
}

//! @brief returns GridOrigin
 //! @return GridOrigin
const std::vector<double> NuTo::StructureGrid::GetGridOrigin() const
{
     return mGridOrigin;
}

 //! @brief returns GridDimension
 //! @return GridDimension
const std::vector<size_t> NuTo::StructureGrid::GetGridDimension() const
{
     return mGridDimension;
}
//! @brief returns GridDimension
 //! @return GridDimension
void  NuTo::StructureGrid::SetGridDimension(std::vector<size_t> &rGridDimension)
{
	mGridDimension=rGridDimension;
}

//! @brief Get NumMaterials
//! @return NumMaterial
const int NuTo::StructureGrid::GetNumMaterials() const
{
    return mNumMaterials;
}
//! @brief Get NumBasisMaterials
//! @return NumBasisMaterial
const int NuTo::StructureGrid::GetNumBasisMaterials() const
{
    return mNumBasisMaterials;
}

//! @brief Set NumBasisMaterials
void NuTo::StructureGrid::SetNumBasisMaterials(int rNumBasisMaterials)
{
	mNumBasisMaterials=rNumBasisMaterials;
}
//! @brief Get CurrentGridNumber
//! @return rCurrentGridNumber
const int NuTo::StructureGrid::GetCurrentGridNumber() const
{
	return mCurrentGridNumber;
}

//! @brief Set CurrentGridNumber
void NuTo::StructureGrid::SetCurrentGridNumber(int rCurrentGridNumber)
{
	mCurrentGridNumber=rCurrentGridNumber;
}

//! @brief Get number of Constraints
 //! @return NumConstraints
const size_t NuTo::StructureGrid::GetNumConstraints() const
{
	 return mNumConstraintDofs;
}

void NuTo::StructureGrid::SetMatrixFreeMethod(bool rMatrixFreeMethod)
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

bool NuTo::StructureGrid::GetMatrixFreeMethod()
{
	return mMatrixFreeMethod;
}

//! @brief set basis element stiffness
//! @param rVoxelSpacing ... element length,
//! @param rPoissonsRatio
//! @return stiffnesMatrix ... basis element stiffness matrix
void NuTo::StructureGrid::SetBasisElementStiffnessMatrix(double rPoissonsRatio,int rBasisMaterialNum)
{
	try
	{
		if (mDimension!=3)
		   throw MechanicsException("[StructureGrid::SetBasisElementStiffnessMatrix] Only 3D is implemented.");
		if(rBasisMaterialNum>=mNumBasisMaterials)
			   throw MechanicsException("[StructureGrid::SetBasisElementStiffnessMatrix] Basis Material number too large.");

		if(rBasisMaterialNum<(int) mLocalCoefficientMatrix0.size())
		   std::cout<<"[StructureGrid::SetBasisElementStiffnessMatrix] This basis element stiffness matrix for material "<<rBasisMaterialNum<<" will be replaces.\n";

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
		elementIncidence(0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(1)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(2)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
		elementIncidence(3)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(4)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(5)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(6)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
		nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
		nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
		elementIncidence(7)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

		// first element create

	   // elementIncidence.Info();
		int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
		myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);
		int mySection1 = myHelpStruc.SectionCreate("VOLUME");
		myHelpStruc.ElementSetSection(myHelpElement,mySection1);

		// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
		NuTo::FullVector<int,Eigen::Dynamic> rows;
		NuTo::FullVector<int,Eigen::Dynamic> coluums;
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrix;
		myHelpStruc.ElementStiffness(0,stiffnessMatrix,rows,coluums );
//		if(mVerboseLevel>3)
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
			NuTo::Brick8N* myElementPointer;
			myElementPointer=static_cast<NuTo::Brick8N*> (myHelpStruc.ElementGetElementPtr(myHelpElement));
			myElementPointer->CalculateDerivativeShapeFunctionsGeometryNatural( rLocalCoordinates, mLocalDerivativeShapeFunctions);
//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" mLocalDerivativeShapeFunctions "<<mLocalDerivativeShapeFunctions[0]<<" "<<mLocalDerivativeShapeFunctions[1]<<" "<<mLocalDerivativeShapeFunctions[2]<<"\n";
//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" mLocalDerivativeShapeFunctions "<<mLocalDerivativeShapeFunctions[21]<<" "<<mLocalDerivativeShapeFunctions[22]<<" "<<mLocalDerivativeShapeFunctions[23]<<"\n";
#ifdef PLANESTRESS
			static_cast<NuTo::LinearElasticEngineeringStress*> (myElementPointer->GetConstitutiveLaw(0))->CalculateCoefficients2DPlainStress(mC11,mC12,mC44);
#else
			static_cast<NuTo::LinearElasticEngineeringStress*> (myElementPointer->GetConstitutiveLaw(0))->CalculateCoefficients3D(mC11,mC12,mC44);
#endif //PLANESTRESS
//			std::cout<<"[StructureGrid] (line "<<__LINE__<<" stiffnessTensorCoefficients "<<mC11<<" "<<mC12<<" "<<" "<<mC44<<"\n";
		}
		else
			throw MechanicsException("[StructureGrid::SetBasisElementStiffnessMatrix] Basis material number is wrong.");
		// for strain and stress computations
		// get shape function derivatives for element center
	}
	catch( ... )
	{
		throw MechanicsException("[StructureGrid::SetBasisElementStiffnessMatrix] Error.");
	}


}

//! @brief set basis edge stiffnesses
//! @param rBasisMaterialNum ... number of material,
void NuTo::StructureGrid::SetBasisEdgeStiffnessMatrices(int rBasisMaterialNum)
{
	const std::vector<double> *baseStiffness=GetLocalCoefficientMatrix0(rBasisMaterialNum);
	if(baseStiffness->size()==0)
        throw MechanicsException("[StructureGrid::SetBasisEdgeStiffnessMatrices] Calculate basis element stiffness matrix first.");
	if((int) mBasisEdgeCoefficientMatrix0.size()>rBasisMaterialNum)
        throw MechanicsException("[StructureGrid::SetBasisEdgeStiffnessMatrices] Material exist already.");
	else if ((int) mBasisEdgeCoefficientMatrix0.size()<rBasisMaterialNum)
        throw MechanicsException("[StructureGrid::SetBasisEdgeStiffnessMatrices] Material number to high.");

	std::vector<double> edgeStiffness(9*64);
	// init edgeStiffness in order,
	// description: node nbr, elem nbr, row nbr (nodes), col nbr (nodes) of base
	int count=0;
	for(int row=0;row<3;++row)
	{
		for(int col=0;col<3;++col)
		{
			count=-1;
			// 0 - 0 - 6  0
			edgeStiffness[(++count)*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+0*3+col];
			// 1 - 0 - 6  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+1*3+col];
			// 1 - 1 - 7  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+0*3+col];
			// 2 - 1 - 7  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+1*3+col];
			// 3 - 0 - 6  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+3*3+col];
			// 3 - 2 - 5  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+0*3+col];
			// 4 - 0 - 6  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+2*3+col];
			// 4 - 1 - 7  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+3*3+col];
			// 4 - 2 - 5  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+1*3+col];
			// 4 - 3 - 4  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+0*3+col];
			// 5 - 1 - 7  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+2*3+col];
			// 5 - 3 - 4  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+1*3+col];
			// 6 - 2 - 5  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+3*3+col];
			// 7 - 2 - 5  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+2*3+col];
			// 7 - 3 - 4  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+3*3+col];
			// 8 - 3 - 4  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+2*3+col];
			// 9 - 0 - 6  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+4*3+col];
			// 9 - 4 - 2  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+0*3+col];
			//10- 0 - 6  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+5*3+col];
			//10 - 1 - 7  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+4*3+col];
			//10 - 4 - 2  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+1*3+col];
			//10 - 5 - 3  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+0*3+col];
			//11 - 1 - 7  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+5*3+col];
			//11 - 5 - 3  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+1*3+col];
			//12 - 0 - 6  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+7*3+col];
			//12 - 2 - 5  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+4*3+col];
			//12 - 4 - 2  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+3*3+col];
			//12 - 6 - 1  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+0*3+col];
			//13 - 0 - 6  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(6*3+row)*24+6*3+col];
			//13 - 1 - 7  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+7*3+col];
			//13 - 2 - 5  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+5*3+col];
			//13 - 3 - 4  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+4*3+col];
			//13 - 4 - 2  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+2*3+col];
			//13 - 5 - 3  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+3*3+col];
			//13 - 6 - 1  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+1*3+col];
			//13 - 7 - 0  0
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+0*3+col];
			//14 - 1 - 7  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(7*3+row)*24+6*3+col];
			//14 - 3 - 4  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+5*3+col];
			//14 - 5 - 3  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+2*3+col];
			//14 - 7 - 0  1
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+1*3+col];
			//15 - 2 - 5  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+7*3+col];
			//15 - 6 - 1  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+3*3+col];
			//16 - 2 - 5  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(5*3+row)*24+6*3+col];
			//16 - 3 - 4  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+7*3+col];
			//16 - 6 - 1  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+2*3+col];
			//16 - 7 - 0  3
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+3*3+col];
			//17 - 3 - 4  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(4*3+row)*24+6*3+col];
			//17 - 7 - 0  2
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+2*3+col];
			//18 - 4 - 2  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+4*3+col];
			//19 - 4 - 2  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+5*3+col];
			//19 - 5 - 3  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+4*3+col];
			//20 - 5 - 3  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+5*3+col];
			//21 - 4 - 2  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+7*3+col];
			//21 - 6 - 1  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+4*3+col];
			//22 - 4 - 2  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(2*3+row)*24+6*3+col];
			//22 - 5 - 3  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+7*3+col];
			//22 - 6 - 1  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+5*3+col];
			//22 - 7 - 0  4
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+4*3+col];
			//23 - 5 - 3  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(3*3+row)*24+6*3+col];
			//23 - 7 - 0  5
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+5*3+col];
			//24 - 6 - 1  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+7*3+col];
			//25 - 6 - 1  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(1*3+row)*24+6*3+col];
			//25 - 7 - 0  7
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+7*3+col];
			//26 - 7 - 0  6
			edgeStiffness[++count*9+(3*row+col)]=(*baseStiffness)[(0*3+row)*24+6*3+col];

		}
	}
	mBasisEdgeCoefficientMatrix0.push_back(edgeStiffness);
}

//! @brief set general neighbor nodes for node-edge based routines
//! @brief in function of grid dimension
void NuTo::StructureGrid::SetNeighborNodesNE()
{
	mNeighborNodesNE.resize(27);
	mNeighborNodesNE[0]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1)-1;
	mNeighborNodesNE[1]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1);
	mNeighborNodesNE[2]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1)+1;
	mNeighborNodesNE[3]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)-1;
	mNeighborNodesNE[4]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1);
	mNeighborNodesNE[5]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)+1;
	mNeighborNodesNE[6]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1)-1;
	mNeighborNodesNE[7]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1);
	mNeighborNodesNE[8]=(int) -(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1)+1;
	mNeighborNodesNE[9]=(int) -(mGridDimension[0]+1)-1;
	mNeighborNodesNE[10]=(int) -(mGridDimension[0]+1);
	mNeighborNodesNE[11]=(int) -(mGridDimension[0]+1)+1;
	mNeighborNodesNE[12]=(int) -1;
	mNeighborNodesNE[13]=(int) 0;
	mNeighborNodesNE[14]=(int) +1;
	mNeighborNodesNE[15]=(int) +(mGridDimension[0]+1)-1;
	mNeighborNodesNE[16]=(int) +(mGridDimension[0]+1);
	mNeighborNodesNE[17]=(int) +(mGridDimension[0]+1)+1;
	mNeighborNodesNE[18]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1)-1;
	mNeighborNodesNE[19]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1);
	mNeighborNodesNE[20]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)-(mGridDimension[0]+1)+1;
	mNeighborNodesNE[21]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)-1;
	mNeighborNodesNE[22]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1);
	mNeighborNodesNE[23]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)+1;
	mNeighborNodesNE[24]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1)-1;
	mNeighborNodesNE[25]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1);
	mNeighborNodesNE[26]=(int) +(mGridDimension[0]+1)*(mGridDimension[1]+1)+(mGridDimension[0]+1)+1;
}

//! @brief ..calculate node numbers at one element
//! @param elementNumber ... number of the element
//! @param nodeNumbers ... nodes at element
void NuTo::StructureGrid::CalculateNodesAtElement(size_t elementNumber,std::vector<size_t>& nodeNumbers)const
{
	size_t rVoxelLocation[3];
	size_t residual1=mVoxelId[elementNumber]%((mGridDimension[0])*(mGridDimension[1]));
	rVoxelLocation[0]=residual1%(mGridDimension[0]);
	rVoxelLocation[1]=residual1/(mGridDimension[0]);
	rVoxelLocation[2]=mVoxelId[elementNumber]/((mGridDimension[0])*(mGridDimension[1]));

	nodeNumbers[0] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[1] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[2] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[3] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[4] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[5] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[6] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[7] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]];

//	std::cout<<" nodesAtElem "<<elementNumber<<": "<<nodeNumbers[0]<<" "<<nodeNumbers[1]<<" "<<nodeNumbers[2]<<" "<<nodeNumbers[3]<<" "<<nodeNumbers[4]<<" "<<nodeNumbers[5]<<" "<<nodeNumbers[6]<<" "<<nodeNumbers[7]<<"\n";
}
//! @brief ..calculate node numbers and node coordinates at one element
//! @param elementNumber ... number of the element
//! @param nodeNumbers ... nodes at element
//! @param nodalCoords ... nodal coordinates at element
void NuTo::StructureGrid::CalculateNodalCoordinatesAtElement(size_t elementNumber,std::vector<size_t>& nodeNumbers,std::vector<double>& nodalCoord)const
{
	size_t rVoxelLocation[3];
	size_t residual1=mVoxelId[elementNumber]%((mGridDimension[0])*(mGridDimension[1]));
	rVoxelLocation[0]=residual1%(mGridDimension[0]);
	rVoxelLocation[1]=residual1/(mGridDimension[0]);
	rVoxelLocation[2]=mVoxelId[elementNumber]/((mGridDimension[0])*(mGridDimension[1]));

	nodeNumbers[0] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[1] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]* (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[2] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[3] = mNodeId[rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[4] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[0]+1) + rVoxelLocation[0]];
	nodeNumbers[5] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[6] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]+1];
	nodeNumbers[7] = mNodeId[(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[0]+1) + rVoxelLocation[0]];


	nodalCoord[0] = rVoxelLocation[0]* mVoxelSpacing[0];
	nodalCoord[1] = rVoxelLocation[1]* mVoxelSpacing[1] ;
	nodalCoord[2] = rVoxelLocation[2]* mVoxelSpacing[2];

	nodalCoord[3] = (rVoxelLocation[0]+1)* mVoxelSpacing[0];
	nodalCoord[4] = rVoxelLocation[1]* mVoxelSpacing[1];
	nodalCoord[5] = rVoxelLocation[2]* mVoxelSpacing[2];

	nodalCoord[6] = (rVoxelLocation[0]+1)* mVoxelSpacing[0];
	nodalCoord[7] = (rVoxelLocation[1]+1)* mVoxelSpacing[1];
	nodalCoord[8] = rVoxelLocation[2]* mVoxelSpacing[2];

	nodalCoord[9] = rVoxelLocation[0]* mVoxelSpacing[0];
	nodalCoord[10] = (rVoxelLocation[1]+1) * mVoxelSpacing[1];
	nodalCoord[11] = rVoxelLocation[2]* mVoxelSpacing[2];

	nodalCoord[12] = rVoxelLocation[0]* mVoxelSpacing[0];
	nodalCoord[13] = (rVoxelLocation[1]+1)* mVoxelSpacing[1];
	nodalCoord[14] = (rVoxelLocation[2]+1)* mVoxelSpacing[2];

	nodalCoord[15] = (rVoxelLocation[0]+1)* mVoxelSpacing[0];
	nodalCoord[16] = rVoxelLocation[1]* mVoxelSpacing[1];
	nodalCoord[17] = (rVoxelLocation[2]+1)* mVoxelSpacing[2];

	nodalCoord[18] = (rVoxelLocation[0]+1)* mVoxelSpacing[0];
	nodalCoord[19] = (rVoxelLocation[1]+1)* mVoxelSpacing[1];
	nodalCoord[20] = (rVoxelLocation[2]+1)* mVoxelSpacing[2];

	nodalCoord[21] = rVoxelLocation[0]* mVoxelSpacing[0];
	nodalCoord[22] = (rVoxelLocation[1]+1) * mVoxelSpacing[1];
	nodalCoord[23] = (rVoxelLocation[2]+1)* mVoxelSpacing[2];

//	std::cout<<" nodesAtElem "<<elementNumber<<": "<<nodeNumbers[0]<<" "<<nodeNumbers[1]<<" "<<nodeNumbers[2]<<" "<<nodeNumbers[3]<<" "<<nodeNumbers[4]<<" "<<nodeNumbers[5]<<" "<<nodeNumbers[6]<<" "<<nodeNumbers[7]<<"\n";
//	std::cout<<"coordsAtElem "<<elementNumber<<": "<<nodalCoord[0]<<" "<<nodalCoord[1]<<" "<<nodalCoord[2]<<" "<<nodalCoord[3]<<" "<<nodalCoord[4]<<" "<<nodalCoord[5]<<" "<<nodalCoord[6]<<" "<<nodalCoord[7]<<"\n";
}

//! @brief set material number at edges for node-edge based routines
void NuTo::StructureGrid::SetMaterialNumberForEdges()
{
	std::vector<double> matOfEdge(0);
	matOfEdge.resize(64*mEdgeId.size(),0);
	int orderEdgesForElem[8][8]={
		{35,39,47,45,57,59,63,62},// 7.
		{27,34,44,41,53,56,61,60},// 6.
		{17,20,32,26,48,49,54,52},// 4.
		{21,23,38,33,50,51,58,55},//5.
		{9,11,15,14,31,37,46,43},//3.
		{5,8,13,12,25,30,42,40}, //2.
		{0,1,6,4,16,18,28,24},//0. elem
		{2,3,10,7,19,22,36,29},//1. elem
	};
	std::vector<int> numElemPerNeigh={1,2,1,2,4,2,1,2,1,2,4,2,4,8,4,2,4,2,1,2,1,2,4,2,1,2,1};
	std::vector<int> edgesNumOfNeigh={0,1,3,4,6,10,12,13,15,17,21,23,27,35,39,41,45,47,48,50,51,53,57,59,60,62,63};
	std::vector<size_t> nodeNumbers(8);
	for(size_t element=0;element<mVoxelId.size();++element)
	{
		CalculateNodesAtElement(element,nodeNumbers);
		// over local 8 nodes to get correct neighbor element
		for(int locNode=0;locNode<8;++locNode)
		{
			// save edge matrix at all edges of element
			for(int node=0;node<8;++node)
			{
				matOfEdge[64*nodeNumbers[locNode]+orderEdgesForElem[locNode][node]]=mYoungsModulus[element];

			}
		}
	}
	mYoungsModulus.assign(matOfEdge.begin(),matOfEdge.end());
}

//! @brief create grid data
//! @param rThresholdMaterialValue ... threshold between material one and two
//! @param imageValues ... vector of image data
//! @param rColorToMaterialData ... vector of material data (Young's Modulus) mapped to color points
void NuTo::StructureGrid::CreateGrid(int rThresholdMaterialValue, std::string fileName,
		 std::vector<double>& rColorToMaterialData)
{
	mImageDataFile=fileName;
	std::vector<int> imageValues (mNumVoxel);
	NuTo::StructureGrid::ImportFromVtkASCIIFile( fileName,imageValues);
	size_t numElems=0;
	size_t numNodes=0;
	size_t numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
	boost::dynamic_bitset<> nodeExist(numGridNodes); //0 = false, all 0 here
	mNodeId.resize(numGridNodes);		//saves the id of the edges of this node

	std::vector<int> allEdgesAtVoxel(mNumVoxel*8);
	std::vector<size_t>globVoxelLocation(mNumVoxel*3);
	for (size_t element=0;element<mNumVoxel;++element)
	{
			size_t numDimxy=element/((mGridDimension[0])*(mGridDimension[1]));
			size_t numDimx=0;
			size_t residual1=element%((mGridDimension[0])*(mGridDimension[1]));
			size_t residual2=0;
			numDimx=residual1/(mGridDimension[0]);
			residual2=residual1%(mGridDimension[0]);
			size_t rVoxelLocation[3];
			rVoxelLocation[0]=residual2;
			rVoxelLocation[1]=numDimx;
			rVoxelLocation[2]=numDimxy;
			globVoxelLocation[element+0]=residual2;
			globVoxelLocation[element+1]=numDimx;
			globVoxelLocation[element+2]=numDimxy;
			//! @brief for grid creation calculate and save temporaly all nodes at all voxels
			allEdgesAtVoxel[8*element+0] = (int) rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0];
			allEdgesAtVoxel[8*element+1] =  (int)rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+2] = (int) rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0] +1;
			allEdgesAtVoxel[8*element+3] = (int) rVoxelLocation[2]*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0];

			allEdgesAtVoxel[8*element+4] = (int) (rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0];
			allEdgesAtVoxel[8*element+5] =  (int)(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + rVoxelLocation[1]     * (mGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+6] = (int) (rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+7] =  (int)(rVoxelLocation[2]+1)*(mGridDimension[0]+1)*(mGridDimension[1]+1) + (rVoxelLocation[1]+1) * (mGridDimension[1]+1) + rVoxelLocation[0];
	}

	for (size_t countVoxels=0;countVoxels<mNumVoxel;++countVoxels)
	{
		if (imageValues[countVoxels]<rThresholdMaterialValue)
		{
			mVoxelId.push_back(countVoxels);

//			std::cout<<" voxel "<<countVoxels<<" elem "<<numElems;
			++numElems;
			mYoungsModulus.push_back(rColorToMaterialData[imageValues[countVoxels]]);

			for (int node=0;node<8;++node)
			{
				nodeExist.set(allEdgesAtVoxel[8*countVoxels+node],true);
//				std::cout<<" nodeExist "<<allEdgesAtVoxel[8*countVoxels+node];
			}
  		}
	}
	mNumBasisMaterials=1;
	mNumMaterials=mYoungsModulus.size();
 	for (size_t node=0;node<numGridNodes;++node)
	{
		if (nodeExist[node])
		{
			mEdgeId.push_back(node);
			mNodeId[node]=numNodes++;
		}
		else
		{
			mNodeId[node]=numGridNodes;
		}
	}
 	for (size_t node=0;node<numGridNodes;++node)
 	{
		if (mNodeId[node]==numGridNodes)
		{
			mNodeId[node]=mEdgeId.size();
		}

 	}

 	// add nodes without element to mData
 	std::vector<size_t> nodeNumbers(8);
 	std::vector<double> nodalCoord(8*3);

 	mRightHandSide.resize(3*(mEdgeId.size()+1),0.0);// initialized with zero

 	if(mVerboseLevel>1)
	{
		std::cout<<"[NuTo::StructureGrid] numElems = "<<numElems<<"\n";
		std::cout<<"[NuTo::StructureGrid] numNodes = "<<mEdgeId.size()<<"\n";
		std::cout<<"[NuTo::StructureGrid] numNodes = "<<numNodes<<"\n";
		std::cout<<"[NuTo::StructureGrid] Young's Modulus elem 0= "<<mYoungsModulus[0]<<"\n";
	}

	if(mVerboseLevel>5)
	{
		std::cout<<"\n[NuTo::StructureGrid] : voxelId ";
		for(size_t i=0;i<numElems;++i)
			std::cout<<mVoxelId[i]<<" ";
		std::cout<<"\n";

		std::cout<<"[NuTo::StructureGrid] : nodeId ";
		for(size_t i=0;i<numGridNodes;++i)
			std::cout<<mNodeId[i]<<" ";
		std::cout<<"\n";

		std::cout<<"[NuTo::StructureGrid] : edgeId ";
		for(size_t i=0;i<numNodes;++i)
			std::cout<<mEdgeId[i]<<" ";
		std::cout<<"\n";
	}

#pragma acc data copyin(mVoxelId,mYoungsModulus,mLocalCoefficientMatrix0,mRightHandSide)
}

//! @brief set displacement boundary conditions
//! @brief all nodes at one plane perpendicular to direction are constrained
//! @param rDirection ... 0 = x,1 = y,2 = z
//! @param rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
//! @param rValue ... value of boundary displacement
//! @param rDisplVector ... output initial displacement vector
void NuTo::StructureGrid::SetDisplacementConstraints(size_t rDirection,size_t *rGridLocation,double rValue,std::vector<double> &rDisplVector)
{
	// rGridLocation is with considering frame nodes
	if(mDofIsConstraint.size()==0)
	{
		mDofIsConstraint.resize(3*mEdgeId.size(),false); //0 = false, all 0 here
		mDisplacements.resize(3*(mEdgeId.size()+1),0.0);// initialized with zero
		mNumConstraintDofs=0;
	}
	size_t gridNode=0;
//	if (mVerboseLevel>1)
//	{
//		std::cout<<"[StructureGrid::SetDisplacementConstraints] constraint regions x: "<<rGridLocation[0]<<" - "<<rGridLocation[1];
//		std::cout<<" y: "<<rGridLocation[2]<<" - "<<rGridLocation[3];
//		std::cout<<" z: "<<rGridLocation[4]<<" - "<<rGridLocation[5];
//	}
	for(size_t dim2=rGridLocation[4];dim2<=rGridLocation[5];++dim2)
	{
		for(size_t dim1=rGridLocation[2];dim1<=rGridLocation[3];++dim1)
		{
			gridNode= (dim2)*(mGridDimension[0]+1)*(mGridDimension[1]+1)
					 +(dim1)*(mGridDimension[0]+1)
					 +     rGridLocation[0]  ;
			for(size_t dim0=rGridLocation[0];dim0<=rGridLocation[1];++dim0)
			{
				if (mNodeId[gridNode]<mEdgeId.size()) //check if node exist and increment counter
				{
//					std::cout<<" grid node "<<gridNode<<" mNodeId[gridNode] "<<mNodeId[gridNode]<<" rDirection "<<rDirection<<"\n ";

					mDofIsConstraint.set(mNodeId[gridNode]*3+rDirection,true);
					mDisplacements[mNodeId[gridNode]*3+rDirection]=rValue;
					++mNumConstraintDofs;
				}
					++gridNode;
			}
		}
	}
	if(mVerboseLevel>3)
		std::cout<<" numConstraints "<<mNumConstraintDofs<<"\n";
	rDisplVector=mDisplacements;
}

//! @brief get DisplacementConstaints
//! @return dynamic_bitset of constraints
const boost::dynamic_bitset<> NuTo::StructureGrid::GetDisplacementConstaints()
{
	return mDofIsConstraint;
}


//! @brief Get LocalCoefficientMatrix0
//! @param NumLocalCoefficientMatrix0 number of stiffness matrix
const std::vector<double>* NuTo::StructureGrid::GetLocalCoefficientMatrix0(int rNumLocalCoefficientMatrix0) const
{
 	if (rNumLocalCoefficientMatrix0<0 || rNumLocalCoefficientMatrix0>=GetNumBasisMaterials())
        throw MechanicsException("[NuTo::StructureGrid::GetLocalCoefficientMatrix0] No valid material number.");
    return &mLocalCoefficientMatrix0[rNumLocalCoefficientMatrix0];
}
//! @brief calculate gradient

void NuTo::StructureGrid::Gradient (std::vector<double>& rValues,std::vector<double>& rGradient)
{
	if(!mMatrixFreeMethod)
		CalculateMatrixVectorProductEBE(rValues,rGradient);
	else
		CalculateMatrixVectorProductNBN(rValues,rGradient);
}
//! @brief ... calculate matrix-vector-product in element-by-element way
//! @brief ... with local vectors
//! @param ... u - parameters input, r - gradient output
void NuTo::StructureGrid::CalculateMatrixVectorProductEBE(std::vector<double>& u,std::vector<double>& r)const
{
	size_t numParas=r.size();
	size_t numElems=mVoxelId.size();
	int dofsElem=24;


	// make sure residual vector is reset with zero //
#pragma acc parallel loop present(r[0:numParas])
{
	for(size_t i=0;i<numParas;++i)
		r[i]=0;
}

	//move local data in loop for acc ?
	#pragma acc firstprivate(residual[0,24],displacement[0,24],nodeNumbers[0,numNodes],youngsModulus)
	std::vector<double> residual (24);
	std::vector<double> displacement(24);
	int numNodes=8;
	std::vector<size_t> nodeNumbers(numNodes);
	double youngsModulus=0;

	//loop over all elements
	#pragma acc parallel loop present(u[0:numParas],r[0:numParas],mVoxelId[0:numElems],mYoungsModulus[0:numElems],	mLocalCoefficientMatrix0[0:576],mDofIsConstraint[0:numParas])
	for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
	{
		//calculate local return vector with all dofs: r=Ku
		#ifndef NODESATELEM
		CalculateNodesAtElement(elementNumber,nodeNumbers);
		#endif
		youngsModulus=mYoungsModulus[elementNumber];

		#pragma acc parallel loop
		for(int node=0;node<numNodes;++node)
		{
			#ifdef NODESATELEM
			nodeNumbers[node]=mNodesAtElem[8*elementNumber+node];
			#endif
//			displacement[3*node]=u.at(3*nodeNumbers[node]);
			displacement[3*node]=u[3*nodeNumbers[node]];
			displacement[3*node+1]=u[3*nodeNumbers[node]+1];
			displacement[3*node+2]=u[3*nodeNumbers[node]+2];
		}
		#pragma acc parallel loop //needed or for each thread initialized already with zero
		for(int i=0;i<dofsElem;++i)
		{
			residual[i]=0;
			#pragma acc parallel loop
			for(int j=0;j<dofsElem;++j)
				residual[i]+=youngsModulus*mLocalCoefficientMatrix0[0][dofsElem*i+j]*displacement[j];
		}

		for(int node=0;node<8;++node)
		{
			if(!mDofIsConstraint[3*nodeNumbers[node]])
				r[3*nodeNumbers[node]]+=residual[3*node];
			if (!mDofIsConstraint[3*nodeNumbers[node]+1])
				r[3*nodeNumbers[node]+1]+=residual[3*node+1];
			if (!mDofIsConstraint[3*nodeNumbers[node]+2])
			      r[3*nodeNumbers[node]+2]+=residual[3*node+2];
//				r.at(3*nodeNumbers[node]+2)+=residual[3*node+2];

		}
	}
}

//! @brief ... calculate matrix-vector-product in node-by-node way
//! @brief ... with local vectors
//! @param ... u - prarmeters input, r - gradient output
void NuTo::StructureGrid::CalculateMatrixVectorProductNBN(std::vector<double> &u,std::vector<double> &r)const
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME
//	if (mVerboseLevel>3)
//		std::cout<<"[StructureGrid::CalculateMatrixVectorProductNBN] "<<std::endl;

	// make sure RightHandSide vector is reset with zero //
	for(size_t i=0;i<r.size();++i)
		r[i]=0;
	//----------------------------------------------//

	size_t numNodes=mEdgeId.size();
//	int  nodeNumbers[64];
	const int nodeNbrOfEdge[64]={0,1,1,2,3,3,4,4,4,4,5,5,6,7,7,8,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,13,13,14,14,14,14,15,15,16,16,16,16,17,17,18,19,19,20,21,21,22,22,22,22,23,23,24,25,25,26};
//	const int nbrOfEdges[27]={1,2,1,2,4,2,1,2,1,2,4,2,4,8,4,2,4,2,1,2,1,2,4,2,1,2,1};

	for (size_t node=0;node<numNodes;++node)
    {
		size_t neighbors[64] ;
		double youngsModulus[64];
		for (size_t i=0;i<64;++i)
		{
			neighbors[i]=mNodeId[mEdgeId[node]+mNeighborNodesNE[nodeNbrOfEdge[i]]];

			youngsModulus[i] = mYoungsModulus[64*node+i];
		}


//		std::cout<<" NBN : grid "<<nodeNumber<<" nodeid "<<node<<" grid neigh: ";
		r[3*node+0]=0;
		r[3*node+1]=0;
		r[3*node+2]=0;
		//27 neighbor nodes
		//calculate local return vector: r=Ku
		for(int i=0;i<64;++i)
		{
//			std::cout<<" dof "<<3*neighbors[i]<<" node "<<node<<" i "<<i<<" edge neigh "<<mEdgeId[node]+mNeighborNodesNE[nodeNbrOfEdge[i]]<<"="<<mEdgeId[node]<<"+"<<mNeighborNodesNE[nodeNbrOfEdge[i]]<<" nodeId " <<mNodeId[mEdgeId[node]+mNeighborNodesNE[nodeNbrOfEdge[i]]]<<"\n";
			r[3*node+0]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i]*u[3*neighbors[i]];
			r[3*node+0]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+1]*u[3*neighbors[i]+1];
			r[3*node+0]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+2]*u[3*neighbors[i]+2];

			r[3*node+1]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+3]*u[3*neighbors[i]];
			r[3*node+1]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+4]*u[3*neighbors[i]+1];
			r[3*node+1]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+5]*u[3*neighbors[i]+2];

			r[3*node+2]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+6]*u[3*neighbors[i]];

			r[3*node+2]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+7]*u[3*neighbors[i]+1];
			r[3*node+2]+=youngsModulus[i]*mBasisEdgeCoefficientMatrix0[0][9*i+8]*u[3*neighbors[i]+2];
		}
		r[3*node+0]*=!mDofIsConstraint[3*node];
		r[3*node+1]*=!mDofIsConstraint[3*node+1];
		r[3*node+2]*=!mDofIsConstraint[3*node+2];

//		std::cout<<"\n";
	}


//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::StructureGrid::CalculateMatrixVectorProductNBN] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief get weighting factor for preconditioner
//! @return rWeight ... weighting factor
double  NuTo::StructureGrid::GetWeightingFactor()const
{
	return mWeightingFactor;
}

//! @brief set weighting factor for preconditioner
//! @param rWeight ... weighting factor
void  NuTo::StructureGrid::SetWeightingFactor(double rWeight)
{
	mWeightingFactor=rWeight;
}

void NuTo::StructureGrid::Hessian(std::vector<double>&  rDiagHessian)
{
	if(mUseDiagHessian)
	{
		HessianDiag(rDiagHessian);
	}
	else
	{
		std::cout<<"Implement\n";
		//set RightHandSide in structure as input
//		MultiGridSolve();
	}

}
//! @brief calculate preconditioned vector z with diagonal preconditioner
//! @param rHessianDiag ... input residual, output z=p*r
void NuTo::StructureGrid::HessianDiag(std::vector<double>& rHessianDiag)
{
#pragma data present(rHessianDiag,mYoungsModulus,mLocalCoefficientMatrix0)
#pragma acc parallel loop
	// input vector is residual which has to preconditioned
	size_t numNodes=mEdgeId.size();
	std::vector<double> rProduct(3*numNodes);
	if(mHessianDiag.empty()==true)
	{
		if (!mMatrixFreeMethod)
		{
			size_t numElems=mVoxelId.size();
	//		std::cout<<"[StructureGrid::HessianDiag] EBE"<<std::endl;
			std::vector<size_t> nodeNumbers(8);

	#pragma acc parallel loop
			for (size_t elementNumber=0;elementNumber<numElems;++elementNumber)
			{
				double elemStiff=0;
				elemStiff=mYoungsModulus[elementNumber]*mLocalCoefficientMatrix0[0][0];
	#ifdef NODESATELEM
				nodeNumbers[node]=mNodesAtElem[8*elementNumber+node];
	#else
				CalculateNodesAtElement(elementNumber,nodeNumbers);
	#endif
				for(int node=0;node<8;++node)
				{
					rProduct[3*nodeNumbers[node]]+=elemStiff;
					rProduct[3*nodeNumbers[node]+1]+=elemStiff;
					rProduct[3*nodeNumbers[node]+2]+=elemStiff;
				}
			}
		}
		else
		{
	//		std::cout<<"[StructureGrid::HessianDiag] NBN"<<std::endl;
			for (size_t nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
			{
				double elemStiff[3]={0,0,0};
				for (size_t edge=28;edge<36;++edge)
				{
						elemStiff[0]+=mYoungsModulus[64*nodeNumber+edge]*mBasisEdgeCoefficientMatrix0[0][9*edge];
						elemStiff[1]+=mYoungsModulus[64*nodeNumber+edge]*mBasisEdgeCoefficientMatrix0[0][9*edge+4];
						elemStiff[2]+=mYoungsModulus[64*nodeNumber+edge]*mBasisEdgeCoefficientMatrix0[0][9*edge+8];
				}
				rProduct[3*nodeNumber]=elemStiff[0];
				rProduct[3*nodeNumber+1]=elemStiff[1];
				rProduct[3*nodeNumber+2]=elemStiff[2];
			}

		}
		size_t numDOFs=mEdgeId.size()*3;
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
//		std::cout << " Hessian: ";
//		for(size_t i=0;i<numNodes*3;++i)
//			std::cout<< rProduct[i]<<" ";
//		std::cout<<"\n";

	}
	for (size_t para=0;para<numNodes*3;++para)
		rHessianDiag[para]*=mHessianDiag[para];
}
void NuTo::StructureGrid::CalcScalingFactors(std::vector<double> &p)
{
	//diagonal scaling with scaling factor
	double scalefactor=GetWeightingFactor();
	if(scalefactor==0.)
	{
#ifdef ENABLE_OPTIMIZE
		if (mUseMisesWielandt)
		{
			SetMisesWielandt(false); // if not, get a infinite loop, for Hessian
			// scale factor is 2/(2-lambda_max-lambda_min) [Meister: Num. lin. GLS] for Jacobi-Relaxation-Method
			NuTo::MisesWielandt myEigenCalculator(GetNumNodes()*3);
			myEigenCalculator.SetVerboseLevel(GetVerboseLevel()+1);
			myEigenCalculator.SetAccuracyGradient(1e-4);
			myEigenCalculator.SetCallback(this);
			myEigenCalculator.Optimize();
			double lambda_max=myEigenCalculator.GetObjective();

	//			double lambda_min=lambda_max/2.;
			// Jacobi-Relaxation-weighting
	//		scalefactor=2./(2-lambda_max-lambda_min);
			// my scale factor
	//		scalefactor=2./(lambda_max*lambda_max);
//			scalefactor=1e-2/lambda_max;
			scalefactor=1./lambda_max;
			// damping Jacobi: lampda of D-1 K Arbenz_2007
	//		scalefactor=4./(3.*lambda_max);
			SetMisesWielandt(true);
		}
		else
#endif //ENABLE_OPTIMIZE
			// 1 is needed for MisesWielandt (Hessian), make sure it is 1 in standard
			scalefactor=1.;
			//   scalefactor=1./(double) mNumParameters;
		SetWeightingFactor(scalefactor);
//		std::cout<<"[StructureGrid] weighting factor "<<scalefactor<<"\n";

	}
	for (size_t count=0; count<GetNumNodes()*3; ++count)
        p[count] *=scalefactor;

}

void NuTo::StructureGrid::GetEngineeringStrain(const std::vector<double> &rDisplacements, std::vector<double> &rEngineeringStrain)const
{
#ifdef PLANESTRESS
	std::cout<<"[StructureGrid::GetEngineeringStrain] PLANESTRESS defined \n";
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
	elementIncidence(0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(1)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(2)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = -mVoxelSpacing[2] * 0.5;
	elementIncidence(3)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(4)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = -mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(5)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(6)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0) = -mVoxelSpacing[0] * 0.5;
	nodeCoordinates(1) = mVoxelSpacing[1] * 0.5;
	nodeCoordinates(2) = mVoxelSpacing[2] * 0.5;
	elementIncidence(7)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	// first element create

   // elementIncidence.Info();
	int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
	int mySection1 = myHelpStruc.SectionCreate("VOLUME");
	myHelpStruc.ElementSetSection(myHelpElement,mySection1);
	Brick8N* myElementPointer=static_cast<Brick8N*> (myHelpStruc.ElementGetElementPtr(myHelpElement));

	size_t rNodesPerElement=8;
	std::vector<size_t> nodes(rNodesPerElement);
	rEngineeringStrain.resize(mVoxelId.size()*6);
	std::vector<double>nodeCoord(rNodesPerElement*3);
	double invJacobian[9], detJac;
	std::vector<double> derivativeShapeFunctionsGlobal(mLocalDerivativeShapeFunctions.size());

#ifdef PLANESTRESS
	for(size_t element=0;element<mVoxelId.size();++element)
	{
		CalculateNodalCoordinatesAtElement(element,nodes,nodeCoord);

		myElementPointer->CalculateJacobian(mLocalDerivativeShapeFunctions,nodeCoord, invJacobian, detJac);

		myElementPointer->CalculateDerivativeShapeFunctionsGlobal(mLocalDerivativeShapeFunctions,invJacobian,
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
	for(size_t element=0;element<mVoxelId.size();++element)
	{
		CalculateNodalCoordinatesAtElement(element,nodes,nodeCoord);

		myElementPointer->CalculateJacobian(mLocalDerivativeShapeFunctions,nodeCoord, invJacobian, detJac);

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
void NuTo::StructureGrid::GetEngineeringStress( std::vector<double>& rEngineeringStrain, std::vector<double>& rEngineeringStress)const
{
#ifdef PLANESTRESS
	std::cout<<"[StructureGrid::GetEngineeringStress] PLANESTRESS defined \n";
#endif
	//! @ToDo check strain caluclated
	size_t rNodesPerElement=8;
	std::vector<int> nodes(rNodesPerElement);
	rEngineeringStress.resize(mVoxelId.size()*6);
#ifdef PLANESTRESS
	for(size_t element=0;element<mVoxelId.size();++element)
	{
			rEngineeringStress[6*element+0]=mYoungsModulus[element]
			                              *(mC11 * rEngineeringStrain[6*element+0]
			                               +mC12 * rEngineeringStrain[6*element+1]);
			rEngineeringStress[6*element+1]=mYoungsModulus[element]
			                              *(mC12 * rEngineeringStrain[6*element+0]
			                               +mC11 * rEngineeringStrain[6*element+1]);
			rEngineeringStress[6*element+2]=0.;
			//Sigma_yz
			rEngineeringStress[6*element+3]=0.;
			//Sigma_xz
			rEngineeringStress[6*element+4]=0.;
			//Sigma_xy
			rEngineeringStress[6*element+5]=mYoungsModulus[element]*mC44*rEngineeringStrain[6*element+5];
	}
#else
	for(size_t element=0;element<mVoxelId.size();++element)
	{
			rEngineeringStress[6*element+0]=mYoungsModulus[element]
			                              *(mC11 * rEngineeringStrain[6*element+0]
			                               +mC12 * rEngineeringStrain[6*element+1]
			                               +mC12 * rEngineeringStrain[6*element+2]);
			rEngineeringStress[6*element+1]=mYoungsModulus[element]
			                              *(mC12 * rEngineeringStrain[6*element+0]
			                               +mC11 * rEngineeringStrain[6*element+1]
			                               +mC12 * rEngineeringStrain[6*element+2]);
			rEngineeringStress[6*element+2]=mYoungsModulus[element]
			                              *(mC12 * rEngineeringStrain[6*element+0]
			                               +mC12 * rEngineeringStrain[6*element+1]
			                               +mC11 * rEngineeringStrain[6*element+2]);
			//Sigma_yz
			rEngineeringStress[6*element+3]=mYoungsModulus[element]*mC44*rEngineeringStrain[6*element+3];
			//Sigma_xz
			rEngineeringStress[6*element+4]=mYoungsModulus[element]*mC44*rEngineeringStrain[6*element+4];
			//Sigma_xy
			rEngineeringStress[6*element+5]=mYoungsModulus[element]*mC44*rEngineeringStrain[6*element+5];
	}
#endif //PLANESTRESS
}
// export to Vtk Datafile
void NuTo::StructureGrid::ExportVTKStructuredDataFile(const std::string& rFilename) const
{
	std::ofstream file(rFilename.c_str());
    if (!file.is_open())
    	throw MechanicsException(std::string("[NuTo::StructureGrid::ExportVTKStructuredDataFile] Error opening file ")+rFilename.c_str());

    // header /////////////////////////////////////////////////////////////////
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "Data file was generated by NuTo" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_POINTS" << std::endl;
    ///////////////////////////////////////////////////////////////////////////
    file << "DIMENSIONS "<<mGridDimension[0]+1<<" "<<mGridDimension[1]+1<<" "<<mGridDimension[2]+1<<"\n";
    file << "SPACING "<<mVoxelSpacing[0]<<" "<<mVoxelSpacing[1]<<" "<<mVoxelSpacing[2]<<"\n";
    file << "ORIGIN "<<mGridOrigin[0]<<" "<<mGridOrigin[1]<<" "<<mGridOrigin[2]<<"\n";
    file << "CELL_DATA "<<mNumVoxel<<"\n";
    file << "SCALARS imageData int 1\n";
    file << "LOOKUP_TABLE default\n";

    using namespace boost::spirit::classic;
   std::ifstream input(mImageDataFile, std::ios::in);
    if (input.is_open() == false)
        throw MechanicsException("[StructureGrid::ExportVTKStructuredDataFile] error opening input file.");
    std::string line;
    for (int count=0;count<10;count++)
        getline (input, line);

    int value;
    size_t count=0;

    while(getline(input,line))
    {
    	std::istringstream iss(line);
    	while(iss >> value && count<mNumVoxel)
    	{
    		file << value<<" \n" ;
    		++count;
    	}
    }
     // close file
	input.close();

	size_t numNodes=mEdgeId.size();
	size_t numEdges=mNodeId.size();
	size_t countNodes=0;
	assert(3*(numNodes+1)==mDisplacements.size());
	file << "POINT_DATA "<<numEdges<<"\n";
	file << "VECTORS displacements double \n";
//			   std::cout<<" VTK mDisplacements "<<mDisplacements.size()<<"\n";
	for (size_t i=0;i<mNodeId.size();++i) // loop over grid points
	{
		if(mNodeId[i]<numNodes)
		{
			file<<mDisplacements[3*mNodeId[i]]<<" " <<mDisplacements[3*mNodeId[i]+1]<<" "<<mDisplacements[3*mNodeId[i]+2]<<"\n";
//			file<<mRightHandSide[3*mNodeId[i]]<<" " <<mRightHandSide[3*mNodeId[i]+1]<<" "<<mRightHandSide[3*mNodeId[i]+2]<<"\n";
			++countNodes;
		}
		else
			file<<"0.0 0.0 0.0\n";
	}
	assert(countNodes==numNodes);

//
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
//	file << "POINT_DATA "<<numEdges<<"\n";
//	file << "VECTORS constraints int \n";
//	for (size_t i=0;i<numEdges;++i) // loop over grid points
//	{
//		if(mNodeId[i]<numNodes)
//		{
//			file<<mDofIsConstraint[3*mNodeId[i]]<<" " <<mDofIsConstraint[3*mNodeId[i]+1]<<" "<<mDofIsConstraint[3*mNodeId[i]+2]<<"\n";
//		}
//		else
//			file<<"0 0 0\n";
//	}
//
//	file << "POINT_DATA "<<numEdges<<"\n";
//	file << "SCALARS id int \n";
//    file << "LOOKUP_TABLE default\n";
//	for (size_t i=0;i<numEdges;++i) // loop over grid points
//	{
//		file<<mNodeId[i]<<"\n";
//	}


    file.close();
}

std::vector<double>&  NuTo::StructureGrid::GetParameters()
{
	return mDisplacements;
}
std::vector<double>&  NuTo::StructureGrid::GetRightHandSide()
{
	return mRightHandSide;
}
void NuTo::StructureGrid::SetParameters(std::vector<double>& rParameters)
{
	mDisplacements=rParameters;
}

void NuTo::StructureGrid::SetRightHandSide(std::vector<double>& rRightHandSide)
{
	mRightHandSide=rRightHandSide;
}
//! @brief MultiGrid routines
//! @brief initialize on coarser grid
//! @param StructureGrid reference
void NuTo::StructureGrid::SetCoarserGridLevel(NuTo::StructureGrid *rCoarseGrid)
{
	std::vector<size_t> rGridDimension=GetGridDimension();
	std::vector<size_t> rOrigDimension=GetGridDimension();

	size_t numElems=GetNumElements();
	std::vector<size_t> rElemId(GetNumVoxels(),numElems);

	for(size_t elementNumber=0;elementNumber<numElems;++elementNumber)
		rElemId[mVoxelId[elementNumber]]=elementNumber;


	if(!mCurrentGridNumber)
		mRightHandSide.resize((GetNumNodes()+1)*3,0.);

	size_t voxel=0;
	size_t edge=0;

	std::vector<size_t> neighborVoxel{0,1,rGridDimension[0],rGridDimension[0]+1,
		(rGridDimension[0])*(rGridDimension[1]),(rGridDimension[0])*(rGridDimension[1])+1,
		(rGridDimension[0])*(rGridDimension[1])+rGridDimension[0],
		(rGridDimension[0])*(rGridDimension[1])+rGridDimension[0]+1};
	rGridDimension[0]*=0.5;
	rGridDimension[1]*=0.5;
	rGridDimension[2]*=0.5;

	// add "right" frame element
	rGridDimension[0]+=1;
	rGridDimension[1]+=1;
	rGridDimension[2]+=1;
	// new structure settings
	rCoarseGrid->SetGridDimension(rGridDimension);
	rCoarseGrid->mNumVoxel=rGridDimension[0]*rGridDimension[1]*rGridDimension[2];
	rCoarseGrid->mGridOrigin[0]=mGridOrigin[0]-mVoxelSpacing[0];
	rCoarseGrid->mGridOrigin[1]=mGridOrigin[1]-mVoxelSpacing[1];
	rCoarseGrid->mGridOrigin[2]=mGridOrigin[2]-mVoxelSpacing[2];
	rCoarseGrid->mVoxelSpacing[0]=2*mVoxelSpacing[0];
	rCoarseGrid->mVoxelSpacing[1]=2*mVoxelSpacing[1];
	rCoarseGrid->mVoxelSpacing[2]=2*mVoxelSpacing[2];
	rCoarseGrid->SetNumBasisMaterials(GetNumBasisMaterials());
	rCoarseGrid->SetWeightingFactor(GetWeightingFactor());
	rCoarseGrid->SetMisesWielandt(GetMisesWielandt());
	rCoarseGrid->SetMatrixFreeMethod(GetMatrixFreeMethod());
	rCoarseGrid->SetCurrentGridNumber(GetCurrentGridNumber()+1);
	rCoarseGrid->SetFineGridPtr(this);
	SetCoarseGridPtr(rCoarseGrid);

	// same matrix for simplicity copy to this grid
	rCoarseGrid->mLocalCoefficientMatrix0=mLocalCoefficientMatrix0;

	size_t coarseGridNumElements=0;
	size_t numGridNodes=(rGridDimension[0]+1)*(rGridDimension[1]+1)*(rGridDimension[2]+1);//all nodes of the grid
	boost::dynamic_bitset<> nodeExist(numGridNodes); //0 = false, all 0 here

	if (!mNumBasisMaterials)
		  throw MechanicsException("[StructureGrid::SetMultiGridStructure] only one Poisson ration possible.");

	for(size_t dim2=1;dim2<rGridDimension[2]-1;++dim2)
	{
		for(size_t dim1=1;dim1<rGridDimension[1]-1;++dim1)
		{
			for(size_t dim0=1;dim0<rGridDimension[0]-1;++dim0)
			{
				voxel= (2*dim2-1)*(rOrigDimension[0])*(rOrigDimension[1])
					 +(2*dim1-1)*(rOrigDimension[0]) + 2*dim0-1  ; // frame +1, fine grid voxel number 2*
				double rYoungsModulus=0.;
				for (size_t i=0;i<8;++i)
				{
					if(rElemId[voxel+neighborVoxel[i]]<numElems) //element exist in fine grid
					{
						//create coarse element
//							rYoungsModulus+=GetElementYoungsModulus(rElemId[voxel+neighborVoxel[i]]); // coarse grid material
							rYoungsModulus+=mYoungsModulus[rElemId[voxel+neighborVoxel[i]]]; // coarse grid material
					}
				}

				rCoarseGrid->mYoungsModulus.push_back(rYoungsModulus*0.125*2); // average and length considered
				++coarseGridNumElements; // one more already
				rCoarseGrid->mVoxelId.push_back((dim2)*(rGridDimension[0])*(rGridDimension[1])
						 +(dim1)*(rGridDimension[0]) + dim0);


				edge=(dim2)*(rGridDimension[0]+1)*(rGridDimension[1]+1)
					 +(dim1)*(rGridDimension[0]+1) + dim0  ;

				nodeExist.set(edge,true);
				nodeExist.set(edge+1,true);
				nodeExist.set(edge+(rGridDimension[0]+1),true);
				nodeExist.set(edge+(rGridDimension[0]+1)+1,true);
				nodeExist.set(edge+(rGridDimension[0]+1)*(rGridDimension[1]+1),true);
				nodeExist.set(edge+(rGridDimension[0]+1)*(rGridDimension[1]+1)+1,true);
				nodeExist.set(edge+(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1),true);
				nodeExist.set(edge+(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1)+1,true);
				//nodes later


			}
		}
	}
	//create nodes
	size_t numNodes=0;
//		mNodeId. neue größe
	for (size_t edge=0;edge<numGridNodes;++edge)
	{
		if (nodeExist[edge])
		{
			rCoarseGrid->mEdgeId.push_back(edge);
			rCoarseGrid->mNodeId.push_back(numNodes++);
		}
		else
			rCoarseGrid->mNodeId.push_back(numGridNodes);
	}
	for (size_t edge=0;edge<numGridNodes;++edge)
	{
		if (rCoarseGrid->mNodeId[edge]==numGridNodes)
			rCoarseGrid->mNodeId[edge]=mEdgeId.size();
	}
	size_t fineEdge=0;
	rCoarseGrid->mDofIsConstraint.resize(numNodes*3,false);
	rCoarseGrid->mDisplacements.resize((numNodes+1)*3,0.);
	rCoarseGrid->mRightHandSide.resize((numNodes+1)*3,0.);
	rCoarseGrid->mFineEdgeId.resize(numNodes);
	rCoarseGrid->mpFineGrid=this;
	for(size_t dim2=1;dim2<rGridDimension[2];++dim2)
	{
		for(size_t dim1=1;dim1<rGridDimension[1];++dim1)
		{
			for(size_t dim0=1;dim0<rGridDimension[0];++dim0)
			{

				fineEdge= (2*dim2-1)*(rOrigDimension[0]+1)*(rOrigDimension[1]+1)
					 +(2*dim1-1)*(rOrigDimension[0]+1) + 2*dim0-1  ; // frame +1, fine grid voxel number 2*
				edge=(dim2)*(rGridDimension[0]+1)*(rGridDimension[1]+1)
							 +(dim1)*(rGridDimension[0]+1) + dim0  ;
				if(rCoarseGrid->mNodeId[edge]<numNodes)
				{
					rCoarseGrid->mFineEdgeId[rCoarseGrid->mNodeId[edge]]=fineEdge;
					rCoarseGrid->mDofIsConstraint.set(3*rCoarseGrid->mNodeId[edge],mDofIsConstraint[3*mNodeId[fineEdge]]);
					rCoarseGrid->mDofIsConstraint.set(3*rCoarseGrid->mNodeId[edge]+1,mDofIsConstraint[3*mNodeId[fineEdge]+1]);
					rCoarseGrid->mDofIsConstraint.set(3*rCoarseGrid->mNodeId[edge]+2,mDofIsConstraint[3*mNodeId[fineEdge]+2]);
					if(rCoarseGrid->mDofIsConstraint[3*rCoarseGrid->mNodeId[edge]])
					{
//						rCoarseGrid->mDisplacements[3*rCoarseGrid->mNodeId[edge]]=mDisplacements[3*mNodeId[fineEdge]];
						rCoarseGrid->mNumConstraintDofs++;
					}
					if(rCoarseGrid->mDofIsConstraint[3*rCoarseGrid->mNodeId[edge]+1])
					{
//						rCoarseGrid->mDisplacements[3*rCoarseGrid->mNodeId[edge]+1]=mDisplacements[3*mNodeId[fineEdge]+1];
						rCoarseGrid->mNumConstraintDofs++;
					}
					if(rCoarseGrid->mDofIsConstraint[3*rCoarseGrid->mNodeId[edge]+2])
					{
//						rCoarseGrid->mDisplacements[3*rCoarseGrid->mNodeId[edge]+2]=mDisplacements[3*mNodeId[fineEdge]+2];
						rCoarseGrid->mNumConstraintDofs++;
					}
				}
			}
		}
	}
	rCoarseGrid->Info();
	mpCoarseGrid=rCoarseGrid;
}

//! @brief MultiGrid routines
//! @brief set pointer to coarser grid
void NuTo::StructureGrid::SetCoarseGridPtr(NuTo::StructureGrid* rCoarseGrid)
{
		mpCoarseGrid=rCoarseGrid;
}

//! @brief MultiGrid routines
//! @brief set pointer to finer grid
void NuTo::StructureGrid::SetFineGridPtr(NuTo::StructureGrid* rFineGrid)
{
		mpFineGrid=rFineGrid;
}



void NuTo::StructureGrid::CalculateMultiGridCorrolations(std::string restrictionType,std::vector<double>& rRestriction,std::vector<double>& rProlongation)
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
void NuTo::StructureGrid::Restriction(std::vector<double> &rRestrictionFactor,
		std::vector<double> & rResidual,std::vector<double> &rResidualCoarser)
{
//	if(CallbackHandlerGrid::mVerboseLevel>2)
//		std::cout<<"[StructureGrid::Restriction] \n";
	std::vector<size_t> rCoarseDimension=mpCoarseGrid->GetGridDimension();
	std::vector<size_t> rFineDimension=GetGridDimension();

	// 27 neighbor nodes
	if (mNeighborNodesNE.empty())
		SetNeighborNodesNE();

	int fineGridEdge=0;
	int rNumNodes=mpCoarseGrid->GetNumNodes();
	rResidualCoarser.resize(rNumNodes*3+3);

	for (int node = 0; node < rNumNodes; ++node)
	{
		std::vector<double> help={0,0,0};
		//corresponding fine grid node
		fineGridEdge=(int) mpCoarseGrid->mFineEdgeId[node];
		if(!mpCoarseGrid->mDofIsConstraint[3*node])
		{
			for (int i = 0; i < 27; ++i)
			{
				if(mNodeId[fineGridEdge+mNeighborNodesNE[i]]<GetNumNodes())
					help[0]+=rResidual[3*mNodeId[fineGridEdge+mNeighborNodesNE[i]]]*rRestrictionFactor[i];
			}
			rResidualCoarser[3*node]=help[0];
		}

		if(!mpCoarseGrid->mDofIsConstraint[3*node+1])
		{
			for (int i = 0; i < 27; ++i)
			{
				if(mNodeId[fineGridEdge+mNeighborNodesNE[i]]<GetNumNodes())
					help[1]+=rResidual[3*mNodeId[fineGridEdge+mNeighborNodesNE[i]]+1]*rRestrictionFactor[i];
			}
			rResidualCoarser[3*node+1]=help[1];

		}
		if(!mpCoarseGrid->mDofIsConstraint[3*node+2])
		{
			for (int i = 0; i < 27; ++i)
			{
				if(mNodeId[fineGridEdge+mNeighborNodesNE[i]]<GetNumNodes())
					help[2]+=rResidual[3*mNodeId[fineGridEdge+mNeighborNodesNE[i]]+2]*rRestrictionFactor[i];
			}
			rResidualCoarser[3*node+2]=help[2];
		}
	}
}

void NuTo::StructureGrid::Prolongation(std::vector<double> &rProlongationFactor,std::vector<double>& rParameters,std::vector<double>& rParametersFine)
{
//	std::cout<<"[StructureGrid::Prolongation] \n";
	size_t numFineNodes=GetNumNodes();
	std::vector<size_t> rCoarseDimension=mpCoarseGrid->GetGridDimension();
	std::vector<size_t> rFineDimension=GetGridDimension();
	// 27 neighbor nodes
	if (mNeighborNodesNE.empty())
		SetNeighborNodesNE();
	size_t numNeighbors=mNeighborNodesNE.size();
	size_t fineGridEdge=0;
	int numNodes=mpCoarseGrid->GetNumNodes();
	for (int node = 0; node < numNodes; ++node)
	{
		//corresponding fine grid node
		fineGridEdge=mpCoarseGrid->mFineEdgeId[node];
		for (size_t i = 0; i < numNeighbors; ++i)
		{
			size_t nodeId=mNodeId[fineGridEdge+mNeighborNodesNE[i]];
			if (nodeId<numFineNodes)
			{
				if(!mDofIsConstraint[3*nodeId])
					rParametersFine[3*nodeId]+=rParameters[3*node]*rProlongationFactor[i];
				if(!mDofIsConstraint[3*nodeId+1])
					rParametersFine[3*nodeId+1]+=rParameters[3*node+1]*rProlongationFactor[i];
				if(!mDofIsConstraint[3*nodeId+2])
					rParametersFine[3*nodeId+2]+=rParameters[3*node+2]*rProlongationFactor[i];
			}
		}
	}
}

void  NuTo::StructureGrid::AnsysInput(std::vector<double> &rDisplVector) const
{
	// open file
	std::ofstream file;
    file.open("ansysInput");
    file<<"!Ansys Input File: \n !DIM: "<<mGridDimension[0]<<" DOFS: "<<mEdgeId.size()*3<<"\n";
    file<<"/prep7 \net,1,solid185 \nkeyopt,1,2,3 \n";

    assert(mNumMaterials);
    assert((int) mYoungsModulus.size()==mNumMaterials);
    int countMat=mNumMaterials;
    for(int i=0;i<countMat;++i)
    {
        file<<"mp,ex,"<<i+1<<","<<mYoungsModulus[i]<<" \n";
        file<<"mp,prxy,"<<i+1<<",.2 \n";

    }
    int node=0;
    for (size_t z=0;z<mGridDimension[2]+1;++z)
    {
    	for (size_t y=0;y<mGridDimension[1]+1;++y)
    	{
    		for (size_t x=0;x<mGridDimension[0]+1;++x)
			{
    			if(mNodeId[node]<mEdgeId.size())
    			{
    				file<<"n,"<<mNodeId[node]+1<<","<<x*mVoxelSpacing[0]<<","<<y*mVoxelSpacing[1]<<","<<z*mVoxelSpacing[2]<<"\n";
   				if (mDofIsConstraint[3*mNodeId[node]])
    				{
    					file<<"nsel,s,node,,"<<mNodeId[node]+1<<"\n";
    				   	file<<"d,all,ux,"<<rDisplVector[3*mNodeId[node]]<<"\n";
     				}
   				if (mDofIsConstraint[3*mNodeId[node]+1])
    				{
    					file<<"nsel,s,node,,"<<mNodeId[node]+1<<"\n";
    				   	file<<"d,all,uy,"<<rDisplVector[3*mNodeId[node]+1]<<"\n";
     				}
   				if (mDofIsConstraint[3*mNodeId[node]+2])
    				{
    					file<<"nsel,s,node,,"<<mNodeId[node]+1<<"\n";
    				   	file<<"d,all,uz,"<<rDisplVector[3*mNodeId[node]+2]<<"\n";
     				}
   			}
    			node++;
			}
    	}
    }
    file<<"alls\n";
    int mat=1;
    file<<"type,1 \n mat,"<<mat<<" \n";
    std::vector<size_t> nodes(8);
    // attention for more materials
    for(int i=0;i<countMat;++i)
    {
		file<<"mat,"<<i+1<<"\n";
		CalculateNodesAtElement(i,nodes);
 		file << "en," << i + 1 << "," << nodes[0] + 1
				<< "," <<nodes[1] + 1 << ","
				<< nodes[2] + 1 << ","
				<< nodes[3] + 1 << ","
				<< nodes[4] + 1 << ","
				<< nodes[5] + 1 << ","
				<< nodes[6] + 1 << ","
				<< nodes[7] + 1 << "\n";

    }
    file.close();
}

void  NuTo::StructureGrid::LSDynaInput() const
{
	// open file
	std::ofstream file;

    file.open("lsdynaInput");
      assert(mNumMaterials);
    assert((int) mYoungsModulus.size()==mNumMaterials);
    int countMat=mNumMaterials;
    // change so, that materials arent doubled
    if(file)
    {
			//material id, youngsmodulus, nu
		for(int i=0;i<countMat;++i)
		{

			file<<i+1<<","<<mYoungsModulus[i]<<" , 0.2\n";

		}
		file<<"\n";
		//node id, x,y,z
		int node=0;
		for (size_t z=0;z<mGridDimension[2]+1;++z)
		{
			for (size_t y=0;y<mGridDimension[1]+1;++y)
			{
				for (size_t x=0;x<mGridDimension[0]+1;++x)
				{
					if(mNodeId[node]<mEdgeId.size())
					{
						file<<mNodeId[node]+1<<","<<x*mVoxelSpacing[0]<<","<<y*mVoxelSpacing[1]<<","<<z*mVoxelSpacing[2]<<"\n";
					}
					node++;
				}
			}
		}
		file<<"\n";
		for(size_t elems=0;elems<mVoxelId.size();++elems)
		{
			std::vector<size_t> nodes(8);
			CalculateNodesAtElement(elems,nodes);
			file <<  elems + 1 << "," <<elems+1 << ","<< nodes[0] + 1
					<< "," <<nodes[1] + 1 << ","
					<< nodes[2] + 1 << ","
					<< nodes[3] + 1 << ","
					<< nodes[4] + 1 << ","
					<< nodes[5] + 1 << ","
					<< nodes[6] + 1 << ","
					<< nodes[7] + 1 << "\n";

		}
			// element id, material id, node 1, 2, 3, ..

    }
    file.close();
}
