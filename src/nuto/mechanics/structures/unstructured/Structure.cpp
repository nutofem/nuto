#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#include <boost/serialization/map.hpp>
#include <cstdint>
#endif // ENABLE_SERIALIZATION

#include <boost/assign/ptr_map_inserter.hpp>

# ifdef _OPENMP
#include <omp.h>
# endif

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/base/Timer.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputDummy.h"
#include "nuto/mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputBlockVectorInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputFullVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"

#include <ANN/ANN.h>
#include <set>

//! @brief constructor
//! @param mDimension number of nodes
NuTo::Structure::Structure(int rDimension) :
        StructureBase(rDimension)
{
}

//! @brief destructor
NuTo::Structure::~Structure()
{
	mElementMap.clear();
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Structure::Info() const
{
    StructureBase::Info();
    NodeInfo(mVerboseLevel);
    ElementInfo(mVerboseLevel);
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Structure::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Structure::save(Archive & ar, const unsigned int version) const {}

template<class Archive>
void NuTo::Structure::saveImplement(Archive & ar) const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start save of structure" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase);
    ar & boost::serialization::make_nvp ("elementMap", mElementMap);
    ar & boost::serialization::make_nvp ("nodeMap", mNodeMap);

    /***************************** Pointer update *****************************/
    // cast the 'mNodeMap' to a map containing pairs (int,uintptr_t)
    std::map<int, std::uintptr_t> mNodeMapCast;
    for (boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        mNodeMapCast.insert(std::pair<int, std::uintptr_t>(it->first, reinterpret_cast<std::uintptr_t>(it->second)));
    }
    ar & boost::serialization::make_nvp("mNodeMapCast", mNodeMapCast);

    // cast the 'mElementMap' to a map containing pairs (int,uintptr_t)
    std::map<int, std::uintptr_t> mElementMapCast;
    for (boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.begin(); it!= mElementMap.end(); it++)
    {
        mElementMapCast.insert(std::pair<int, std::uintptr_t>(it->first, reinterpret_cast<std::uintptr_t>(it->second)));
    }
    ar & boost::serialization::make_nvp("mElementMapCast", mElementMapCast);

#ifdef _OPENMP
    int size = mMIS.size();
    ar & boost::serialization::make_nvp("mMIS_size", size);
    int i = 0;
    for (std::vector<std::vector<ElementBase*>>::iterator it =  mMIS.begin(); it!=mMIS.end(); it++, i++)
    {
        int size = it->size();
        const std::uintptr_t* mMISAddress = reinterpret_cast<const std::uintptr_t*>(it->data());
        std::string name = "size_" + std::to_string(i);
        ar & boost::serialization::make_nvp(name.c_str(), size);
        name = "mMIS_" + std::to_string(i);
        ar & boost::serialization::make_nvp(name.c_str(), boost::serialization::make_array(mMISAddress, size));
    }
#endif

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish save of structure" << std::endl;
#endif
}


//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Save (const std::string &filename, std::string rType ) const
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

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

            saveImplement(oba);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);

            saveImplement(oxa);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);

            saveImplement(ota);
        }
        else
        {
            throw MechanicsException ( "[NuTo::Structure::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception &e )
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

void NuTo::Structure::SaveUpdate (const std::string &filename, std::string rType ) const
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::Structure::SaveUpdate] Error opening file.");
        }
        // write data to file
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            // TODO!
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            // TODO!
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            // TODO!
        }
        else
        {
            throw MechanicsException ( "[NuTo::Structure::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception &e )
    {
        std::string s ( std::string ( "[NuTo::Structure::SaveUpdate]File save exception in boost - " ) + std::string ( e.what() ) );
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

// serializes the class
template void NuTo::Structure::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Structure::load(Archive & ar, const unsigned int version){}

template<class Archive>
void NuTo::Structure::loadImplement(Archive & ar)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start load of structure" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase);

    ar & boost::serialization::make_nvp ("elementMap", mElementMap);
    ar & boost::serialization::make_nvp ("nodeMap", mNodeMap);

    /***************************** Pointer update *****************************/
    // node
    std::map<int, std::uintptr_t> mNodeMapCast;
    ar & boost::serialization::make_nvp("mNodeMapCast", mNodeMapCast);
    std::map<std::uintptr_t, std::uintptr_t> mNodeMapOldNewPtr;
    std::map<int, std::uintptr_t>::iterator itCastNode = mNodeMapCast.begin();
    for (boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++, itCastNode++)
    {
        mNodeMapOldNewPtr.insert(std::pair<std::uintptr_t, std::uintptr_t>(itCastNode->second, reinterpret_cast<std::uintptr_t>(it->second)));
    }


    // element
    std::map<int, std::uintptr_t> mElementMapCast;
    ar & boost::serialization::make_nvp("mElementMapCast", mElementMapCast);
    std::map<std::uintptr_t, std::uintptr_t> mElementMapOldNewPtr;
    std::map<int, std::uintptr_t>::iterator itCastElement = mElementMapCast.begin();
    for (boost::ptr_map<int,ElementBase>::iterator it = mElementMap.begin(); it!= mElementMap.end(); it++, itCastElement++)
    {
        mElementMapOldNewPtr.insert(std::pair<std::uintptr_t, std::uintptr_t>(itCastElement->second, reinterpret_cast<std::uintptr_t>(it->second)));
    }

#ifdef _OPENMP
    int size = 0;
    ar & boost::serialization::make_nvp("mMIS_size", size);
    mMIS.resize(size);
    int i = 0;
    for (std::vector<std::vector<ElementBase*>>::iterator it =  mMIS.begin(); it!=mMIS.end(); it++, i++)
    {
        int size = 0;
        std::string name = "size_" + std::to_string(i);
        ar & boost::serialization::make_nvp(name.c_str(), size);
        std::uintptr_t* mMISAddress = new std::uintptr_t[size];

        name = "mMIS_" + std::to_string(i);
        ar & boost::serialization::make_nvp(name.c_str(), boost::serialization::make_array(mMISAddress, size));

        for(int i = 0; i < size; i++)
        {
            std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mElementMapOldNewPtr.find(mMISAddress[i]);
            if (itCast!=mElementMapOldNewPtr.end())
            {
                mMISAddress[i] = itCast->second;
            }
            else
                throw MechanicsException("[NuTo::Structure::loadImplement] The ElementBase-Pointer could not be updated.");
        }

        it->assign(reinterpret_cast<ElementBase**>(&mMISAddress[0]), reinterpret_cast<ElementBase**>(&mMISAddress[size]));
    }
#endif

    // update the node ptr in elements
    for(boost::ptr_map<int, ElementBase>::iterator itElements = mElementMap.begin(); itElements!=mElementMap.end(); itElements++)
    {
        itElements->second->SetNodePtrAfterSerialization(mNodeMapOldNewPtr);
    }

    // update the nonlocal elements in elementData
    for(boost::ptr_map<int, ElementBase>::iterator itElements = mElementMap.begin(); itElements!=mElementMap.end(); itElements++)
    {
        itElements->second->GetDataPtr()->SetElementPtrAfterSerialization(mElementMapOldNewPtr);
    }

    // exchange node pointer in constraints
    for (boost::ptr_map<int,ConstraintBase>::iterator itConstraints=mConstraintMap.begin(); itConstraints!=mConstraintMap.end(); itConstraints++)
    {
        // cast the Address to a NodeBase-Pointer
        itConstraints->second->SetNodePtrAfterSerialization(mNodeMapOldNewPtr);
    }

    // exchange node pointer in loads
    for (boost::ptr_map<int,LoadBase>::iterator itLoads=mLoadMap.begin(); itLoads!=mLoadMap.end(); itLoads++)
    {
        // cast the Address to a NodeBase-Pointer
        itLoads->second->SetNodePtrAfterSerialization(mNodeMapOldNewPtr);
    }

    // exchange element pointer in loads
    for (boost::ptr_map<int,LoadBase>::iterator itLoads=mLoadMap.begin(); itLoads!=mLoadMap.end(); itLoads++)
    {
        // cast the Address to a NodeBase-Pointer
        itLoads->second->SetElementPtrAfterSerialization(mElementMapOldNewPtr);
    }

    // exchange node AND element pointer in groups
    std::map<std::uintptr_t, std::uintptr_t> mNodeAndElementMapOldNewPtr(mNodeMapOldNewPtr);
    mNodeAndElementMapOldNewPtr.insert(mElementMapOldNewPtr.begin(), mElementMapOldNewPtr.end());
    for (boost::ptr_map<int,GroupBase>::iterator itGroups=mGroupMap.begin(); itGroups!=mGroupMap.end(); itGroups++)
    {
        // cast the Address to a NodeBase-Pointer
        itGroups->second->SetNodePtrAfterSerialization(mNodeAndElementMapOldNewPtr);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish load of structure" << std::endl;
#endif
}

//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Restore (const std::string &filename, std::string rType )
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
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

            loadImplement(oba);
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

            loadImplement(oxa);
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

            loadImplement(ota);
        }
        else
        {
            throw MechanicsException ( "[NuTo::Structure::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception &e )
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

void NuTo::Structure::RestoreUpdate (const std::string &filename, std::string rType )
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::Structure::RestoreUpdate] Error opening file.");
        }
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            // TODO!
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            // TODO!
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            // TODO!
        }
        else
        {
            throw MechanicsException ( "[NuTo::Structure::RestoreUpdate]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception &e )
    {
        std::string s ( std::string ( "[NuTo::Structure::RestoreUpdate] File save exception in boost - " ) + std::string ( e.what() ) );
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


void NuTo::Structure::Evaluate(const NuTo::ConstitutiveInputMap& rInput, std::map<eStructureOutput, StructureOutputBase*> &rStructureOutput)
{
    std::string outputs = " ";
    for (auto it : rStructureOutput)
        outputs += StructureOutputToString(it.first) + " ";

    Timer timer(std::string(__FUNCTION__) + outputs, GetShowTime(), GetLogger());
    try
    {
        if (rStructureOutput.empty())
            return;     // ! ---> may occur if matrices have been identified as constant

        if (this->mNodeNumberingRequired)
            throw MechanicsException(__PRETTY_FUNCTION__, "Node numbering required! Call NodeBuildGlobalDofs() first.");

        // build global tmp static data
        if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
            throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");

        for (auto iteratorOutput : rStructureOutput)
        {
            iteratorOutput.second->SetZero();
        }

        eError errorGlobal = eError::SUCCESSFUL;
#ifdef _OPENMP
        if (mNumProcessors!=0)
        {
            omp_set_num_threads(mNumProcessors);
        }

        if (mMIS.size()==0)
        {
            CalculateMaximumIndependentSets();
        }
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
#pragma omp parallel shared(rStructureOutput) //firstprivate(elementOutputMap)
            {
#endif // _OPENMP
        // The allocation of the elementOutputMap is inside the openmp block
        // since the every thread needs a copy of the map.
        // This special case cannot (to my knowledge) be handled with the
        // omp firstprivate directive, since a copy of a shared_ptr is
        // not a deep copy of the underlying data - which makes perfectly sense.

        // BEWARE (!!!) Do not perform a SetZero on the rStructureOutput here
        // since it will remove the allocation of other MIS.

        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;

        // allocate element outputs and resize the structure outputs
        for (auto iteratorOutput : rStructureOutput)
        {
            switch (iteratorOutput.first)
            {
            case NuTo::eStructureOutput::HESSIAN0:
            {
                elementOutputMap[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(mDofStatus);
                break;
            }
            case NuTo::eStructureOutput::HESSIAN1:
            {
                elementOutputMap[Element::eOutput::HESSIAN_1_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(mDofStatus);
                break;
            }
            case NuTo::eStructureOutput::HESSIAN2:
            {
                elementOutputMap[Element::eOutput::HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(mDofStatus);
                break;
            }
            case NuTo::eStructureOutput::HESSIAN2_LUMPED:
            {
                elementOutputMap[Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockVectorDouble>(mDofStatus);
                break;
            }
            case NuTo::eStructureOutput::INTERNAL_GRADIENT:
            {
                elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(mDofStatus);
                break;
            }
            case NuTo::eStructureOutput::UPDATE_STATIC_DATA:
            {
                elementOutputMap[Element::eOutput::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();
                break;
            }
            default:
            {
                throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + std::string("] Output request not implemented."));
            }
            }
        }
        elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(mDofStatus);
        elementOutputMap[Element::eOutput::GLOBAL_COLUMN_DOF] = std::make_shared<ElementOutputBlockVectorInt>(mDofStatus);
#ifdef _OPENMP
        for (auto elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = *elementIter;

#else
        for (auto elementIter : this->mElementMap)
        {
            ElementBase* elementPtr = elementIter->second;
#endif

            // calculate element contribution
            //bool symmetryFlag = false;
            eError error = elementPtr->Evaluate(rInput, elementOutputMap);

            if (error != eError::SUCCESSFUL)
            {
                if (errorGlobal == eError::SUCCESSFUL)
                {
                    errorGlobal = error;
                } else if (errorGlobal != error)
                {
                    throw MechanicsException(__PRETTY_FUNCTION__, "elements have returned multiple different error codes, can't handle that.");
                }
            }

            const auto & elementVectorGlobalDofsRow = elementOutputMap.at(Element::eOutput::GLOBAL_ROW_DOF)->GetBlockFullVectorInt();
            const auto & elementVectorGlobalDofsColumn = elementOutputMap.at(Element::eOutput::GLOBAL_COLUMN_DOF)->GetBlockFullVectorInt();

            for (auto& iteratorOutput : rStructureOutput)
            {
                StructureOutputBase* structureOutput = iteratorOutput.second;

                switch (iteratorOutput.first)
                {
                case NuTo::eStructureOutput::HESSIAN0:
                {
                    const auto& elementMatrix = elementOutputMap.at(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE)->GetBlockFullMatrixDouble();
                    structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(elementPtr, elementMatrix, elementVectorGlobalDofsRow, elementVectorGlobalDofsColumn, mToleranceStiffnessEntries, GetDofStatus().HasInteractingConstraints());
                    break;
                }
                case NuTo::eStructureOutput::HESSIAN1:
                {
                    const auto& elementMatrix = elementOutputMap.at(Element::eOutput::HESSIAN_1_TIME_DERIVATIVE)->GetBlockFullMatrixDouble();
                    structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(elementPtr, elementMatrix, elementVectorGlobalDofsRow, elementVectorGlobalDofsColumn, mToleranceStiffnessEntries, GetDofStatus().HasInteractingConstraints());
                    break;
                }

                case NuTo::eStructureOutput::HESSIAN2:
                {
                    const auto& elementMatrix = elementOutputMap.at(Element::eOutput::HESSIAN_2_TIME_DERIVATIVE)->GetBlockFullMatrixDouble();
                    structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(elementPtr, elementMatrix, elementVectorGlobalDofsRow, elementVectorGlobalDofsColumn, mToleranceStiffnessEntries, true); // always calculate the KJ and KK
                                                                                                                                                                                                              // since its most likely only needed once,
                                                                                                                                                                                                              // and causes troubles in the test files.
                    break;
                }

                case NuTo::eStructureOutput::HESSIAN2_LUMPED:
                {
                    const auto& elementVector = elementOutputMap.at(Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->GetBlockFullVectorDouble();

                    structureOutput->AsStructureOutputBlockMatrix().AddElementVectorDiagonal(elementVector, elementVectorGlobalDofsRow, mToleranceStiffnessEntries);
                    break;
                }

                case NuTo::eStructureOutput::INTERNAL_GRADIENT:
                {
                    const auto& elementVector = elementOutputMap.at(Element::eOutput::INTERNAL_GRADIENT)->GetBlockFullVectorDouble();

                    structureOutput->AsStructureOutputBlockVector().AddElementVector(elementVector, elementVectorGlobalDofsRow);
                    break;
                }

                case NuTo::eStructureOutput::UPDATE_STATIC_DATA:
                    break;

                default:
                {
                    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, StructureOutputToString(iteratorOutput.first) + " requested but not implemented.");
                }
                }
            }

#ifdef _OPENMP
        }
    }   // end loop over elements
}   // end parallel region
}   // end loop over independent sets
#else
        }   // end loop over elements
#endif

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error evaluating the structure.");
        throw e;
    }

}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(int rConstitutiveId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented");
    /*
    boost::ptr_map<int, ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveId);
    if (itConstitutive == mConstitutiveLawMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    try
    {
        BuildNonlocalData(itConstitutive->second);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error calculating nonlocal data.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error calculating nonlocal data.");
    }
     */
}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(const ConstitutiveBase* rConstitutive)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented");
    /*
    double R(rConstitutive->GetParameterDouble(Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS));
    double R2(R * R);
    std::vector<ElementBase*> indexElement;
    std::vector<int> indexIp;
    std::vector<double> indexIpVolume;

    // build up search tree with all integration points
    boost::ptr_map<int, ElementBase>::iterator elementIter;
    for (elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        if (elementPtr == 0)
            continue;

//TODO check element type
//TODO check constitutive type
        Eigen::VectorXd ipVolume = elementPtr->GetIntegrationPointVolume();

        //calculate element contribution and iterate over all integration points
        for (int theIp = 0; theIp < elementPtr->GetNumIntegrationPoints(); theIp++)
        {
            //theWeight = elementIter->second->GetGlobalIntegrationPointWeight(theIp);
            if (&elementPtr->GetConstitutiveLaw(theIp) == rConstitutive)
            {
                elementPtr->DeleteNonlocalElements();
                indexElement.push_back(elementPtr);
                indexIp.push_back(theIp);
                indexIpVolume.push_back(ipVolume(theIp,0));
            }
        }
    }

    // build kd_tree
    ANNpointArray dataPoints;
    ANNkd_tree* kdTree;

    if (mDimension == 2)
    {
        dataPoints = annAllocPts(indexIp.size(), 2);
        for (unsigned int count = 0; count < indexIp.size(); count++)
        {
            ANNpoint thePoint = dataPoints[count];
            Eigen::Vector3d coordinates = indexElement[count]->GetGlobalIntegrationPointCoordinates(indexIp[count]);
            //the third parameter is probably zero, but in order to avoid writing another routine ...
            thePoint[0] = coordinates[0];
            thePoint[1] = coordinates[1];
        }
        kdTree = new ANNkd_tree(dataPoints, indexIp.size(), 2);
    }
    else
    {
        dataPoints = annAllocPts(indexIp.size(), 3);
        for (unsigned int count = 0; count < indexIp.size(); count++)
        {
            ANNpoint thePoint = dataPoints[count];
            Eigen::Vector3d coordinates = indexElement[count]->GetGlobalIntegrationPointCoordinates(indexIp[count]);
            thePoint[0] = coordinates[0];
            thePoint[1] = coordinates[1];
            thePoint[2] = coordinates[2];
        }
        kdTree = new ANNkd_tree(dataPoints, indexIp.size(), 2);
    }

    // find the neighbors in radius R
    unsigned int allocatedResultPoints = 100;
    ANNidxArray nnIdx = new ANNidx[allocatedResultPoints];
    ANNdistArray dists = new ANNdist[allocatedResultPoints];

    for (unsigned int theIp = 0; theIp < indexIp.size(); theIp++)
    {
        ElementBase* elementPtr = indexElement[theIp];
        int localIpNumber = indexIp[theIp];
        unsigned int numNeighborPoints = 0;
        do
        {
            numNeighborPoints = kdTree->annkFRSearch(dataPoints[theIp], R2, allocatedResultPoints, nnIdx, dists, 0);
            if (numNeighborPoints > allocatedResultPoints)
            {
                allocatedResultPoints = numNeighborPoints;
                delete[] nnIdx;
                delete[] dists;
                nnIdx = new ANNidx[allocatedResultPoints];
                dists = new ANNdist[allocatedResultPoints];
                numNeighborPoints = 0;
            }
        } while (numNeighborPoints == 0);

        //calculate total weight to scale the weights
        double totalVolume(0);
        for (unsigned int theNeighborPoint = 0; theNeighborPoint < numNeighborPoints; theNeighborPoint++)
        {
            if (dists[theNeighborPoint] < R2)
            {
                dists[theNeighborPoint] = (1 - dists[theNeighborPoint] / R2);
                dists[theNeighborPoint] *= dists[theNeighborPoint];
                totalVolume += indexIpVolume[nnIdx[theNeighborPoint]] * dists[theNeighborPoint];
            }
            else
                dists[theNeighborPoint] = 0.;

        }
        totalVolume = 1. / totalVolume;

        //add all the nonlocal integration points to the nonlocal data of theIp
        for (unsigned int theNeighborPoint = 0; theNeighborPoint < numNeighborPoints; theNeighborPoint++)
        {
            int theNeighborIndex(nnIdx[theNeighborPoint]);
            elementPtr->SetNonlocalWeight(localIpNumber, indexElement[theNeighborIndex], indexIp[theNeighborIndex], indexIpVolume[theNeighborIndex] * dists[theNeighborPoint] * totalVolume);
        }
    }

    //just check the sum of the weights to be one
    for (elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        if (elementPtr == 0)
            continue;

        for (int theIp = 0; theIp < elementPtr->GetNumIntegrationPoints(); theIp++)
        {
            const std::vector<const NuTo::ElementBase*>& nonlocalElements(elementPtr->GetNonlocalElements());

            double sumW(0.);
            for (int countNonlocalElement = 0; countNonlocalElement < (int) nonlocalElements.size(); countNonlocalElement++)
            {
                const std::vector<double>& weights(elementPtr->GetNonlocalWeights(theIp, countNonlocalElement));

                assert((int )weights.size() == nonlocalElements[countNonlocalElement]->GetNumIntegrationPoints());

                //Go through all the integration points
                for (int theIP = 0; theIP < (int) weights.size(); theIP++)
                {
                    //and add up to nonlocal equivalent plastic strain
                    sumW += weights[theIP];
                }
            }
        }
    }

    delete kdTree;
    if (dataPoints != 0)
        annDeallocPts(dataPoints);

    delete[] nnIdx;
    delete[] dists;

    annClose();
     if(EXIT_SUCCESS != 0)
     {
     INTERPRET_INTERN error_mess("ELEMENT_BUILD_NL_ELEMENTS_ANN: Error using ANN library.");
     return -1;
     }
     */
}

void NuTo::Structure::CalculateInitialValueRates(NuTo::TimeIntegrationBase& rTimeIntegrationScheme)
{
    assert(mNumTimeDerivatives==1 && "Using this function for 0 time derivatives seems to make no sense. More than one time derivative is not implemented so far!");


    constexpr const unsigned int maxIterations = 20;


    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    rTimeIntegrationScheme.CalculateStaticAndTimeDependentExternalLoad();





    // declare necessary variables
    StructureOutputBlockMatrix Hessian_1(mDofStatus,true);

    StructureOutputBlockVector delta_dof_dt1 (mDofStatus,true);

    StructureOutputBlockVector dof_dt1(mDofStatus,true);
    StructureOutputBlockVector residual(mDofStatus, true);
    StructureOutputBlockVector trialResidual(mDofStatus, true);
    StructureOutputBlockVector intForce(mDofStatus, true);
    StructureOutputBlockVector extForce(mDofStatus, true);


    // declare and fill output map for structure
    std::map<eStructureOutput, StructureOutputBase*> StructureOutputsTrial;
    StructureOutputsTrial[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    StructureOutputsTrial[eStructureOutput::HESSIAN1] = &Hessian_1;


    dof_dt1 = NodeExtractDofValues(1);

    ConstitutiveInputMap StructureInputs;
    StructureInputs[Constitutive::eInput::CALCULATE_INITIALIZE_VALUE_RATES] = nullptr;

    Evaluate(StructureInputs,StructureOutputsTrial);

    extForce = rTimeIntegrationScheme.CalculateCurrentExternalLoad(0);

    residual = intForce - extForce;


    unsigned int iteration = 0;


    while(residual.J.CalculateInfNorm()>rTimeIntegrationScheme.GetToleranceResidual())
    {
        ++iteration;
        if(iteration > maxIterations)
            throw MechanicsException(__PRETTY_FUNCTION__,"No convergence while solving for initial value rates!");
        trialResidual = intForce - extForce;


        trialResidual.J -= Hessian_1.JJ * delta_dof_dt1.J;

        delta_dof_dt1.J = SolveBlockSystem(Hessian_1.JJ,residual.J);

        dof_dt1.J += delta_dof_dt1.J;

        NodeMergeDofValues(1,dof_dt1.J,dof_dt1.K);

        Evaluate(StructureInputs,StructureOutputsTrial);

        residual = intForce -extForce;

    }
}

NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ImportFromGmsh(const std::string& rFileName)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> ids;

    try
    {
        ids = ImportFromGmshAux(rFileName);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error importing from Gmsh.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error importing from Gmsh.");
    }
    return ids;
}

class gmsh_node
{
public:
    gmsh_node() :
            id(0)
    {
        this->Coordinates[0] = 0.;
        this->Coordinates[1] = 0.;
        this->Coordinates[2] = 0.;
    }
    ;
    unsigned int id;
    double Coordinates[3];
};

class gmsh_element
{
public:
    gmsh_element() :
            id(0), type(0)
    {
    }
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

NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ImportFromGmshAux(const std::string& rFileName)
{
    const unsigned int num_elm_nodes[24] =
    { 0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 0, 10, 0, 15 };

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
    std::ifstream file(rFileName.c_str(), std::ios::in);
    if (file.is_open() == false)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error opening input file for read access.");
    }

    std::vector<gmsh_node> nodes;
    std::vector<gmsh_element> elements;

    // read first line
    std::string line;
    getline(file, line);

    // get file format
    std::getline(file, line);
    unsigned int mayor_version, minor_version, binary, double_size;

    std::string::iterator itFirst(line.begin());
    bool match(phrase_parse(itFirst, line.end(), uint_[ref(mayor_version) = qi::_1] >> '.' >> uint_[ref(minor_version) = qi::_1] >> uint_[ref(binary) = qi::_1] >> uint_[ref(double_size) = qi::_1], ascii::space));

    if (!match || itFirst != line.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error reading file information.");
    }

    //std::cout << "version " << mayor_version <<"," << minor_version  << " binary " << binary << " double_size " << double_size << std::endl;

    if (mayor_version != 2)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Invalid file format.");
    }

    if (binary == 0) // read ASCII - file
    {
        getline(file, line);
        if (line != "$EndMeshFormat")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndMeshFormat not found.");
        }

        // begin node section
        getline(file, line);
        if (line != "$Nodes")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$Nodes not found.");
        }

        // read number of nodes
        getline(file, line);
        itFirst = line.begin();
        match = phrase_parse(itFirst, line.end(), uint_[ref(num_nodes) = qi::_1], ascii::space);

        if (!match || itFirst != line.end())
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "error reading number of nodes.");
        }

        // read node data
        nodes.resize(num_nodes);
        for (unsigned int node_count = 0; node_count < num_nodes; node_count++)
        {
            getline(file, line);
            itFirst = line.begin();
            match = phrase_parse(itFirst, line.end(), uint_[ref(nodes[node_count].id) = qi::_1] >> double_[ref(nodes[node_count].Coordinates[0]) = qi::_1] >> double_[ref(nodes[node_count].Coordinates[1]) = qi::_1] >> double_[ref(nodes[node_count].Coordinates[2]) = qi::_1], ascii::space);

            if (!match || itFirst != line.end())
            {
                throw MechanicsException(__PRETTY_FUNCTION__, "error reading node data.");
            }
        }

        // end node section
        getline(file, line);
        if (line != "$EndNodes")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndNodes not found.");
        }

        // begin element section
        getline(file, line);
        if (line != "$Elements")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$Elements not found.");
        }

        // read number of elements
        getline(file, line);
        itFirst = line.begin();
        match = phrase_parse(itFirst, line.end(), uint_[ref(num_elements) = qi::_1], ascii::space);

        if (!match || itFirst != line.end())
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "error reading number of elements.");
        }

        // read element data
        elements.resize(num_elements);
        std::vector<unsigned int> tmp_elem_data;
        for (unsigned int element_count = 0; element_count < num_elements; element_count++)
        {
            // read data to vector
            getline(file, line);
            itFirst = line.begin();
            match = phrase_parse(itFirst, line.end(), *(uint_[push_back(boost::phoenix::ref(tmp_elem_data), qi::_1)]), ascii::space);

            if (!match || itFirst != line.end())
            {
                throw MechanicsException(__PRETTY_FUNCTION__, "error reading element data.");
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
            for (unsigned int tag_count = 0; tag_count < num_tags; tag_count++)
            {
                elements[element_count].tags[tag_count] = *iter;
                iter++;
            }

            // get number of nodes and nodes
            unsigned int cur_num_elm_nodes = num_elm_nodes[elements[element_count].type];
            elements[element_count].nodes.resize(cur_num_elm_nodes);
            for (unsigned int elem_node_count = 0; elem_node_count < cur_num_elm_nodes; elem_node_count++)
            {
                elements[element_count].nodes[elem_node_count] = *iter;
                iter++;
            }

            // check size
            if (iter != tmp_elem_data.end())
            {
                throw MechanicsException(__PRETTY_FUNCTION__, "invalid number of element data.");
            }
            tmp_elem_data.clear();
        }

        // end element section
        getline(file, line);
        if (line != "$EndElements")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndElements not found.");
        }
    }

    else // binary format
    {
        if (double_size != sizeof(double))
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Invalid size of double.");
        }

        // close file and open as binary
        file.close();
        file.open(rFileName.c_str(), std::ios::in | std::ios::binary);
        if (file.is_open() == false)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Error opening input file for read access.");
        }

        // read the first two lines again
        getline(file, line);
        getline(file, line);

        // check size of integer
        int one;
        file.read((char *) &one, sizeof(int));
        if (one != 1)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Invalid binary format.");
        }
        file.seekg(1, std::ios::cur);

        getline(file, line);
        if (line != "$EndMeshFormat")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndMeshFormat not found.");
        }

        // begin node section
        getline(file, line);
        if (line != "$Nodes")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$Nodes not found.");
        }

        // read number of nodes
        file >> num_nodes;
        getline(file, line); //endl

        // read node data
        nodes.resize(num_nodes);

        for (unsigned int node_count = 0; node_count < num_nodes; node_count++)
        {
            file.read((char *) &nodes[node_count].id, sizeof(int));
            file.read((char *) nodes[node_count].Coordinates, 3 * sizeof(double));
        }

        //endl
        getline(file, line);

        //$EndNodes
        getline(file, line);

        if (line != "$EndNodes")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndNodes not found.");
        }

        // begin element section
        getline(file, line);
        if (line != "$Elements")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$Elements not found.");
        }

        // read number of elements
        file >> num_elements;
        getline(file, line); //endl

        elements.resize(num_elements);

        // read element data
        unsigned int element_type;
        unsigned int num_elm_follow;
        unsigned int num_tags;
        unsigned int cur_num_elm_nodes;

        for (unsigned int elemCount = 0; elemCount < num_elements; elemCount++)
        {

            // Read element type
            file.read((char *) &element_type, sizeof(int));

            // Read num of Elem with the same header
            file.read((char *) &num_elm_follow, sizeof(int));

            // set num_elemt_node
            cur_num_elm_nodes = num_elm_nodes[element_type];

            // Read numOfTags
            file.read((char *) &num_tags, sizeof(int));

            for (unsigned int indexH = 0; indexH < num_elm_follow; indexH++)
            {

                // set element type
                elements[elemCount].type = element_type;

                // read element number
                file.read((char *) &elements[elemCount].id, sizeof(int));

                elements[elemCount].tags.resize(num_tags);
                elements[elemCount].nodes.resize(cur_num_elm_nodes);

                //read tags
                for (unsigned int tagCount = 0; tagCount < num_tags; tagCount++)
                    file.read((char *) &elements[elemCount].tags[tagCount], sizeof(int));

                //read nodes
                for (unsigned int nodeCount = 0; nodeCount < cur_num_elm_nodes; nodeCount++)
                    file.read((char *) &elements[elemCount].nodes[nodeCount], sizeof(int));

                elemCount += indexH;
            }

        }
        getline(file, line); //endl

        // end element section
        getline(file, line);
        if (line != "$EndElements")
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "$EndElements not found.");
        }
    }
///////////////end binary
    //create the nodes
    NuTo::FullVector<double, Eigen::Dynamic> coordinates;
    switch (mDimension)
    {
    case 2:
        coordinates.Resize(2);
        break;
    case 3:
        coordinates.Resize(3);
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Only implemented for 2D and 3D.");
    }
    std::map<int, int> newNodeNumber;
    for (unsigned int nodeCount = 0; nodeCount < nodes.size(); nodeCount++)
    {
        coordinates(0) = nodes[nodeCount].Coordinates[0];
        coordinates(1) = nodes[nodeCount].Coordinates[1];
        if (mDimension == 3)
            coordinates(2) = nodes[nodeCount].Coordinates[2];
        newNodeNumber[nodes[nodeCount].id] = NodeCreate(coordinates);
    }

    // allocate data structure for group id and interpolation type id
    std::map<int, std::set<int>> groupInterpolationIds;

    NuTo::FullVector<int, Eigen::Dynamic> nodeNumbers;
    for (unsigned int elementCount = 0; elementCount < elements.size(); elementCount++)
    {
        nodeNumbers.Resize(elements[elementCount].nodes.size());
        for (unsigned int countNode = 0; countNode < elements[elementCount].nodes.size(); countNode++)
            nodeNumbers(countNode) = newNodeNumber[elements[elementCount].nodes[countNode]];
        //std::cout << "element " << elementCount << " with nodes " << nodeNumbers.Trans() << std::endl;
        Interpolation::eShapeType shapeType;
        Interpolation::eTypeOrder typeOrder;

        switch (elements[elementCount].type)
        {
        case 1: // 	2-node line in 2d/3d
            shapeType = Interpolation::eShapeType::TRUSSXD;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;
        case 2: // 3-node triangle.
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 3: // 4-node quadrangle.
            shapeType = Interpolation::eShapeType::QUAD2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 4: // 4-node tetrahedron.
            shapeType = Interpolation::eShapeType::TETRAHEDRON3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

        case 5: // 8-node hexahedron.
            shapeType = Interpolation::eShapeType::BRICK3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT1;
            break;

//    	case 6: // 6-node prism.

//    	case 7: // 5-node pyramid.

    	case 8: // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
            shapeType = Interpolation::eShapeType::TRUSSXD;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            //ordering is different than in gmsh, fix this first
            {
            NuTo::FullVector<int, Eigen::Dynamic> nodeNumbersCopy(nodeNumbers);
            nodeNumbers(0)  = nodeNumbersCopy(0);
            nodeNumbers(1)  = nodeNumbersCopy(2);
            nodeNumbers(2)  = nodeNumbersCopy(1);
            }
            break;
        case 9: // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

        case 10: // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
            shapeType = Interpolation::eShapeType::QUAD2D;
            typeOrder = Interpolation::eTypeOrder::LOBATTO2;
            break;

        case 11: // 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
            shapeType = Interpolation::eShapeType::TETRAHEDRON3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

    	case 12: // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
    	{
            shapeType = Interpolation::eShapeType::BRICK3D;
            typeOrder = Interpolation::eTypeOrder::LOBATTO2;
            //ordering is different than in gmsh, fix this first
            NuTo::FullVector<int, Eigen::Dynamic> nodeNumbersGmsh(nodeNumbers);
            nodeNumbers(0)  = nodeNumbersGmsh(4);
            nodeNumbers(1)  = nodeNumbersGmsh(16);
            nodeNumbers(2)  = nodeNumbersGmsh(5);
            nodeNumbers(3)  = nodeNumbersGmsh(10);
            nodeNumbers(4)  = nodeNumbersGmsh(21);
            nodeNumbers(5)  = nodeNumbersGmsh(12);
            nodeNumbers(6)  = nodeNumbersGmsh(0);
            nodeNumbers(7)  = nodeNumbersGmsh(8);
            nodeNumbers(8)  = nodeNumbersGmsh(1);
            nodeNumbers(9)  = nodeNumbersGmsh(17);
            nodeNumbers(10) = nodeNumbersGmsh(25);
            nodeNumbers(11) = nodeNumbersGmsh(18);
            nodeNumbers(12) = nodeNumbersGmsh(22);
            nodeNumbers(13) = nodeNumbersGmsh(26);
            nodeNumbers(14) = nodeNumbersGmsh(23);
            nodeNumbers(15) = nodeNumbersGmsh(9);
            nodeNumbers(16) = nodeNumbersGmsh(20);
            nodeNumbers(17) = nodeNumbersGmsh(11);
            nodeNumbers(18) = nodeNumbersGmsh(7);
            nodeNumbers(19) = nodeNumbersGmsh(19);
            nodeNumbers(20) = nodeNumbersGmsh(6);
            nodeNumbers(21) = nodeNumbersGmsh(15);
            nodeNumbers(22) = nodeNumbersGmsh(24);
            nodeNumbers(23) = nodeNumbersGmsh(14);
            nodeNumbers(24) = nodeNumbersGmsh(3);
            nodeNumbers(25) = nodeNumbersGmsh(13);
            nodeNumbers(26) = nodeNumbersGmsh(2);
            break;
    	}
//    	case 13: // 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).

//    	case 14: // 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).

//    	case 15: // 1-node point.

        case 16: // 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
            shapeType = Interpolation::eShapeType::QUAD2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

        case 17: // 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
            shapeType = Interpolation::eShapeType::BRICK3D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT2;
            break;

//    	case 18: // 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).

//    	case 19: // 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

//    	case 20: // 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)

        case 21: // 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT3;
            break;

//    	case 22: // 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)

        case 23: // 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
            shapeType = Interpolation::eShapeType::TRIANGLE2D;
            typeOrder = Interpolation::eTypeOrder::EQUIDISTANT4;
            break;

//    	case 24: // 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)

//    	case 25: // 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)

//    	case 26: // 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)

//    	case 27: // 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)

//    	case 28: // 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)

//    	case 29: // 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)

//    	case 30: // 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)

//    	case 31: // 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)

//    	case 92: // 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)

//    	case 93: //	125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)

        default:
            std::cout << "element type in gmsh " << elements[elementCount].type << std::endl;
            throw MechanicsException(__PRETTY_FUNCTION__, "Element type not implemented in the import routine.");
        }

        // get gmsh group id and create a corresponding nuto group if needed
        int groupId = elements[elementCount].tags[0]; // NuTo groupId == gmsh groupId. // This might cause errors if groups exist before the gmsh import.
        boost::ptr_map<int, GroupBase>::iterator itGroup(mGroupMap.find(groupId));
        if (itGroup == mGroupMap.end())
            GroupCreate(groupId, NuTo::eGroupId::Elements);

        int interpolationTypeId;

        // find the group in the id map, add it if needed
        auto itGroupInterpolation = groupInterpolationIds.find(groupId);
        if (itGroupInterpolation == groupInterpolationIds.end())
        {
            // create a new interpolation type
            interpolationTypeId = InterpolationTypeCreate(Interpolation::ShapeTypeToString(shapeType));
            InterpolationTypeAdd(interpolationTypeId, Node::eDof::COORDINATES, typeOrder);

            // add it to the map
            groupInterpolationIds[groupId].insert(interpolationTypeId);
        }
        else
        {
            bool exists = false;
            for (int id : itGroupInterpolation->second)
            {
                // check if the corresponding interpolation is already in the group
                const InterpolationType& interpolationType = mInterpolationTypeMap.at(id);
                Interpolation::eShapeType groupShapeType = interpolationType.GetShapeType();
                Interpolation::eTypeOrder groupTypeOrder = interpolationType.Get(Node::eDof::COORDINATES).GetTypeOrder();
                if (groupShapeType == shapeType and groupTypeOrder == typeOrder)
                {
                    exists = true;
                    interpolationTypeId = id;
                    break;
                }
            }

            if (not exists)
            {
                // create a new interpolation type
                interpolationTypeId = InterpolationTypeCreate(Interpolation::ShapeTypeToString(shapeType));
                InterpolationTypeAdd(interpolationTypeId, Node::eDof::COORDINATES, typeOrder);

                // add it to the map
                groupInterpolationIds[groupId].insert(interpolationTypeId);
            }



        } // now, a valid group id and a valid interpolation type exists

        int elementId = ElementCreate(interpolationTypeId, nodeNumbers);
        mGroupMap.find(groupId)->second->AddMember(elementId, ElementGetElementPtr(elementId));
    }

    // translate groupInterpolationIds to NuToIntMatrix
    int numGroups = groupInterpolationIds.size();
    int numMaxInterpolationTypesPerGroup = 0;
    for (auto groupInterpolationId : groupInterpolationIds)
    {
        int size = groupInterpolationId.second.size();
        numMaxInterpolationTypesPerGroup = std::max(numMaxInterpolationTypesPerGroup, size);
    }

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> ids(numGroups, numMaxInterpolationTypesPerGroup+1); // +1 since first column is group id
    ids.fill(-1); // invalid value
    int iGroup = 0;
    for (auto groupInterpolationId : groupInterpolationIds)
    {
        ids(iGroup, 0) = groupInterpolationId.first; // group id
        std::set<int> groupInterpolationTypes = groupInterpolationId.second;
        int iIP = 1;
        for (int groupInterpolationType : groupInterpolationTypes)
        {
            ids(iGroup, iIP) = groupInterpolationType;
            iIP++;
        }
        iGroup++;
    }

    return ids;
}


//! @brief copy and move the structure
//! most of the data is kept, but e.g. nonlocal data and
//! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
void NuTo::Structure::CopyAndTranslate(NuTo::FullVector<double, Eigen::Dynamic>& rOffset)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    try
    {
        std::map<NodeBase*, NodeBase*> old2NewNodePointer;
        std::map<ElementBase*, ElementBase*> old2NewElementPointer;
        CopyAndTranslate(rOffset, old2NewNodePointer, old2NewElementPointer);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error translating and copying structure.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error translating and copying structure.");
    }
}

void NuTo::Structure::CopyAndTranslate(NuTo::FullVector<double, Eigen::Dynamic>& rOffset, std::map<NodeBase*, NodeBase*>& rOld2NewNodePointer, std::map<ElementBase*, ElementBase*>& rOld2NewElementPointer)
{
    if (rOffset.GetNumRows() != mDimension)
        throw MechanicsException(__PRETTY_FUNCTION__, "offset has to have the same dimension as the structure.");
    if (rOffset.GetNumColumns() != 1)
        throw MechanicsException(__PRETTY_FUNCTION__, "offset has to have a single column.");

    std::vector<NodeBase*> nodeVector;
    GetNodesTotal(nodeVector);
    for (unsigned int countNode = 0; countNode < nodeVector.size(); countNode++)
    {
        NodeBase* newNode = nodeVector[countNode]->Clone();
        rOld2NewNodePointer[nodeVector[countNode]] = newNode;

        //find unused integer id
        int id(mNodeMap.size());
        boost::ptr_map<int, NodeBase>::iterator it = mNodeMap.find(id);
        while (it != mNodeMap.end())
        {
            id++;
            it = mNodeMap.find(id);
        }

        // add node to map
        this->mNodeMap.insert(id, newNode);

        newNode->Set(Node::eDof::COORDINATES, nodeVector[countNode]->Get(Node::eDof::COORDINATES) + rOffset);

    }
    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired = true;

    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    std::set<ConstitutiveBase*> constitutiveWithNonlocalData;
    for (unsigned int countElement = 0; countElement < elements.size(); countElement++)
    {
        ElementBase* oldElementPtr = elements[countElement];
        int numNodes = oldElementPtr->GetNumNodes();
        std::vector<NodeBase*> nodeVector(numNodes);
        for (int countNode = 0; countNode < numNodes; countNode++)
        {
            nodeVector[countNode] = rOld2NewNodePointer[oldElementPtr->GetNode(countNode)];
        }
        // find interpolation type
        int interpolationTypeId = 0;
        const InterpolationType& interpolationTypeOld = oldElementPtr->GetInterpolationType();
        for (auto it = mInterpolationTypeMap.begin(); it != mInterpolationTypeMap.end(); it++)
            if ((it->second) == &interpolationTypeOld)
            {
                interpolationTypeId = it->first;
                break;
            }

        int newElementId = ElementCreate(interpolationTypeId, nodeVector);
        ElementBase* newElementPtr = ElementGetElementPtr(newElementId);
        rOld2NewElementPointer[oldElementPtr] = newElementPtr;

        //set integration type
        const IntegrationTypeBase& integrationType = oldElementPtr->GetIntegrationType();
        newElementPtr->SetIntegrationType(integrationType);

        //set section
        const SectionBase& section = oldElementPtr->GetSection();
        newElementPtr->SetSection(section);

        //set constitutive model
        ConstitutiveBase& constitutive = oldElementPtr->GetConstitutiveLaw(0);
        newElementPtr->SetConstitutiveLaw(constitutive);
    }

#ifdef _OPENMP
    //there seems to be a problem with the nearest neighbor search library
#pragma omp critical
#endif
    {
        //rebuild nonlocal data
//        for (auto it = constitutiveWithNonlocalData.begin(); it != constitutiveWithNonlocalData.end(); it++)
//            BuildNonlocalData(*it);
    }
}

void NuTo::Structure::NuToSerializeSave(SerializeStreamOut& rStream)
{
    // be super carefull to symmetrically implement the same stuff to NuToSerializeLoad(...)

    // serialize nodes
    for (int i = 0; i < GetNumTimeDerivatives(); ++i)
    {
        auto nodalValues = NodeExtractDofValues(i);
        rStream << nodalValues.J;
        rStream.Separator();
        rStream << nodalValues.K;
        rStream.Separator();
    }

    // serialize element static data
    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    for (ElementBase* element : elements)
    {
        rStream << element->GetIPData();
        rStream.Separator();
    }
}

void NuTo::Structure::NuToSerializeLoad(SerializeStreamIn& rStream)
{
    // serialize nodes
    for (int i = 0; i < GetNumTimeDerivatives(); ++i)
    {
        auto nodalValues = NodeExtractDofValues(i);
        rStream >> nodalValues.J;
        rStream.Separator();
        rStream >> nodalValues.K;
        rStream.Separator();
        NodeMergeDofValues(i, nodalValues.J, nodalValues.K);
    }

    // serialize element static data
    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    for (ElementBase* element : elements)
    {
        rStream >> element->GetIPData();
        rStream.Separator();
    }
}



#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Structure)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
