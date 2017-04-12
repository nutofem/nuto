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

#include "mechanics/structures/unstructured/Structure.h"

#include "base/Timer.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputDummy.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"
#include "mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "mechanics/elements/ElementOutputFullVectorDouble.h"
#include "mechanics/elements/ElementOutputVectorInt.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationBase.h"

#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/Assembler.h"

#include <ANN/ANN.h>
#include <set>

NuTo::Structure::Structure(int rDimension) :
        StructureBase(rDimension)
{
}

NuTo::Structure::~Structure()
{
	mElementMap.clear();
}

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
        throw;
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
        throw;
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
        throw;
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
        throw;
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

    if (rStructureOutput.empty())
        return;     // ! ---> may occur if matrices have been identified as constant

    NodeBuildGlobalDofs();

    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    for (auto iteratorOutput : rStructureOutput)
    {
        iteratorOutput.second->SetZero();
    }

#ifdef _OPENMP
    std::string exceptionMessage = "";
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
            elementOutputMap[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
            break;
        }
        case NuTo::eStructureOutput::HESSIAN1:
        {
            elementOutputMap[Element::eOutput::HESSIAN_1_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
            break;
        }
        case NuTo::eStructureOutput::HESSIAN2:
        {
            elementOutputMap[Element::eOutput::HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
            break;
        }
        case NuTo::eStructureOutput::HESSIAN2_LUMPED:
        {
            elementOutputMap[Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
            break;
        }
        case NuTo::eStructureOutput::INTERNAL_GRADIENT:
        {
            elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
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
    // calculate element contribution
    elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());
    elementOutputMap[Element::eOutput::GLOBAL_COLUMN_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());
#ifdef _OPENMP
    for (auto elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
    {
#pragma omp single nowait
        {
            ElementBase* elementPtr = *elementIter;
            // in OpenMP, exceptions may not leave the parallel region
            try
            {
                elementPtr->Evaluate(rInput, elementOutputMap);
            }
            catch (std::exception& e)
            {
                exceptionMessage = e.what();
            }

#else
    for (auto elementIter : this->mElementMap)
    {
        ElementBase* elementPtr = elementIter->second;
        elementPtr->Evaluate(rInput, elementOutputMap);
#endif

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

    if (exceptionMessage != "")
        throw MechanicsException(exceptionMessage);
#else
    }   // end loop over elements
#endif

}


void NuTo::Structure::CalculateInitialValueRates(NuTo::TimeIntegrationBase& rTimeIntegrationScheme)
{
    assert(mNumTimeDerivatives==1 && "Using this function for 0 time derivatives seems to make no sense. More than one time derivative is not implemented so far!");


    constexpr const unsigned int maxIterations = 20;


    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    rTimeIntegrationScheme.CalculateStaticAndTimeDependentExternalLoad();





    // declare necessary variables
    StructureOutputBlockMatrix Hessian_1(GetDofStatus(),true);

    StructureOutputBlockVector delta_dof_dt1 (GetDofStatus(),true);

    StructureOutputBlockVector dof_dt1(GetDofStatus(),true);
    StructureOutputBlockVector residual(GetDofStatus(), true);
    StructureOutputBlockVector trialResidual(GetDofStatus(), true);
    StructureOutputBlockVector intForce(GetDofStatus(), true);
    StructureOutputBlockVector extForce(GetDofStatus(), true);


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


std::vector<std::pair<int, int>> NuTo::Structure::ImportFromGmsh(const std::string& rFileName)
{
    return MeshCompanion::ImportFromGmsh(*this, rFileName);
}


void NuTo::Structure::CopyAndTranslate(Eigen::VectorXd& rOffset)
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
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error translating and copying structure.");
    }
}

void NuTo::Structure::CopyAndTranslate(Eigen::VectorXd& rOffset, std::map<NodeBase*, NodeBase*>& rOld2NewNodePointer, std::map<ElementBase*, ElementBase*>& rOld2NewElementPointer)
{
    if (rOffset.rows() != mDimension)
        throw MechanicsException(__PRETTY_FUNCTION__, "offset has to have the same dimension as the structure.");
    if (rOffset.cols() != 1)
        throw MechanicsException(__PRETTY_FUNCTION__, "offset has to have a single column.");

    std::vector<NodeBase*> nodeVector;
    GetNodesTotal(nodeVector);
    for (auto& node : nodeVector)
    {
        NodeBase* newNode = node->Clone();
        rOld2NewNodePointer[node] = newNode;

        // find unused integer id
        int id(mNodeMap.size());
        boost::ptr_map<int, NodeBase>::iterator it = mNodeMap.find(id);
        while (it != mNodeMap.end())
        {
            id++;
            it = mNodeMap.find(id);
        }

        // add node to map
        this->mNodeMap.insert(id, newNode);

        newNode->Set(Node::eDof::COORDINATES, node->Get(Node::eDof::COORDINATES) + rOffset);
    }
    //renumbering of dofs for global matrices required
    GetAssembler().SetNodeVectorChanged();

    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    std::set<ConstitutiveBase*> constitutiveWithNonlocalData;
    for (auto oldElementPtr : elements)
    {
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
        std::shared_ptr<const Section> section = oldElementPtr->GetSection();
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
