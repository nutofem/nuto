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
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>


#ifdef SHOW_TIME
#include <ctime>
#endif

# ifdef _OPENMP
#include <omp.h>
# endif

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>

extern "C" {
#include <dSFMT.h>
}


#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NLobatto.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureOutputBase.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include "nuto/visualize/VisualizeComponent.h"
#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#endif // ENABLE_VISUALIZE


NuTo::StructureBase::StructureBase(int rDimension)  : NuTo::NuToObject::NuToObject()
{
    if (rDimension!=1 && rDimension!=2 && rDimension!=3)
    {
        throw MechanicsException("[StructureBase::StructureBase] The dimension of a structure is either 1, 2 or 3.");
    }
    mDimension = rDimension;
    mPrevTime = 0.;
    mTime = 0.;
    mNumActiveDofs = 0;
    mNumDofs   = 0;
    mNodeNumberingRequired = true;
    mNumExtrapolatedCycles = 0;

    mMappingIntEnum2String.resize(NuTo::IntegrationType::NumIntegrationTypes);
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss1Ip]=
        NuTo::IntegrationType1D2NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss2Ip]=
        NuTo::IntegrationType1D2NGauss2Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss3Ip]=
        NuTo::IntegrationType1D2NGauss3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss4Ip]=
        NuTo::IntegrationType1D2NGauss4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss5Ip]=
        NuTo::IntegrationType1D2NGauss5Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NLobatto3Ip]=
        NuTo::IntegrationType1D2NLobatto3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NLobatto4Ip]=
        NuTo::IntegrationType1D2NLobatto4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NLobatto5Ip]=
        NuTo::IntegrationType1D2NLobatto5Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss13Ip]=
        NuTo::IntegrationType2D3NGauss13Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss16Ip]=
        NuTo::IntegrationType2D3NGauss16Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss1Ip]=
        NuTo::IntegrationType2D3NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss3Ip]=
        NuTo::IntegrationType2D3NGauss3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss4Ip]=
            NuTo::IntegrationType2D3NGauss4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss6Ip]=
            NuTo::IntegrationType2D3NGauss6Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss12Ip]=
                NuTo::IntegrationType2D3NGauss12Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss12IpDetail]=
                NuTo::IntegrationType2D3NGauss12IpDetail::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NGauss1Ip]=
        NuTo::IntegrationType2D4NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NGauss4Ip]=
        NuTo::IntegrationType2D4NGauss4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NGauss9Ip]=
        NuTo::IntegrationType2D4NGauss9Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NLobatto9Ip]=
        NuTo::IntegrationType2D4NLobatto9Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NLobatto16Ip]=
        NuTo::IntegrationType2D4NLobatto16Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NLobatto25Ip]=
        NuTo::IntegrationType2D4NLobatto25Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D4NGauss1Ip]=
        NuTo::IntegrationType3D4NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D4NGauss4Ip]=
        NuTo::IntegrationType3D4NGauss4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NGauss1Ip]=
        NuTo::IntegrationType3D8NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip]=
        NuTo::IntegrationType3D8NGauss2x2x2Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NLobatto3x3x3Ip]=
        NuTo::IntegrationType3D8NLobatto<3>::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NLobatto4x4x4Ip]=
        NuTo::IntegrationType3D8NLobatto<4>::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NLobatto5x5x5Ip]=
        NuTo::IntegrationType3D8NLobatto<5>::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NBoundaryGauss3Ip] =
        NuTo::IntegrationType1D2NBoundaryGauss3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType0DBoundary] =
        NuTo::IntegrationType0DBoundary::GetStrIdentifierStatic();

    mNumLoadCases = 1;

    mNumTimeDerivatives = 0;

    mHaveTmpStaticData = false;
    mUpdateTmpStaticDataRequired = true;
    mToleranceStiffnessEntries = 0.;
    mHessianConstant[0] = false;  // consider only the stiffness matrix to be variable, damping and mass are constant
    mHessianConstant[1] = true;   // as a consequence, the inertia term is not considered in the internal gradient routine
    mHessianConstant[2] = true;   // similar the damping term

#ifdef _OPENMP
    mUseMIS = true;
    // then the environment variable is used
    mNumProcessors = 1;
#endif // _OPENMP

}

int NuTo::StructureBase::GetDimension()const
{
    return mDimension;
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::StructureBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::StructureBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::StructureBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::StructureBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::StructureBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::StructureBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::StructureBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    mLogger << "start serialization of structure base" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject);
    ar & BOOST_SERIALIZATION_NVP(mNumTimeDerivatives);
    ar & BOOST_SERIALIZATION_NVP(mDimension);
    ar & BOOST_SERIALIZATION_NVP(mConstitutiveLawMap);
    ar & BOOST_SERIALIZATION_NVP(mConstraintMap);
    ar & BOOST_SERIALIZATION_NVP(mNumLoadCases);
    ar & BOOST_SERIALIZATION_NVP(mLoadMap);
    ar & BOOST_SERIALIZATION_NVP(mGroupMap);
    ar & BOOST_SERIALIZATION_NVP(mIntegrationTypeMap);
    ar & BOOST_SERIALIZATION_NVP(mSectionMap);
    ar & BOOST_SERIALIZATION_NVP(mMappingIntEnum2String);

//    & BOOST_SERIALIZATION_NVP(mVisualizeComponents)

    ar & BOOST_SERIALIZATION_NVP(mNumDofs);
    ar & BOOST_SERIALIZATION_NVP(mNumActiveDofs);
    ar & BOOST_SERIALIZATION_NVP(mNodeNumberingRequired);
    ar & BOOST_SERIALIZATION_NVP(mConstraintMatrix);
    ar & BOOST_SERIALIZATION_NVP(mConstraintRHS);
    ar & BOOST_SERIALIZATION_NVP(mHaveTmpStaticData);
    ar & BOOST_SERIALIZATION_NVP(mUpdateTmpStaticDataRequired);
    ar & BOOST_SERIALIZATION_NVP(mToleranceStiffnessEntries);
#ifdef _OPENMP
    ar & BOOST_SERIALIZATION_NVP(mUseMIS);
    ar & BOOST_SERIALIZATION_NVP(mMIS);
    ar & BOOST_SERIALIZATION_NVP(mNumProcessors);
#endif // _OPENMP
    ar & BOOST_SERIALIZATION_NVP(mLogger);
    ar & boost::serialization::make_nvp("mLogger", mLogger);
    ar & BOOST_SERIALIZATION_NVP(mHessianConstant);
#ifdef DEBUG_SERIALIZATION
    mLogger << "finish serialization of structure base" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureBase::Info()const
{
    mLogger << "dimension : " << mDimension << "\n";

    mLogger << "number of time derivatives : " << mNumTimeDerivatives << "\n";

    mLogger << "num dofs : " << mNumDofs << "\n";
    mLogger << "num active dofs : " << mNumActiveDofs << "\n";
    // print info for sections
    SectionInfo(mVerboseLevel);

    // print info for groups
    GroupInfo(mVerboseLevel);
}

//! @brief ... number of time derivatives for the nodes (0 : static, 1: velocities, 2: accelerations)
void NuTo::StructureBase::SetNumTimeDerivatives(int rNumTimeDerivatives)
{
	if (rNumTimeDerivatives<0 || rNumTimeDerivatives>2)
        throw NuTo::MechanicsException("[NuTo::StructureBase::SetNumTimeDerivatives] number of time derivatives is either 0, 1 or 2.");

	mNumTimeDerivatives = rNumTimeDerivatives;
}

//! @brief ... return number of time derivatives (0 : static, 1: velocities, 2: accelerations)
int NuTo::StructureBase::GetNumTimeDerivatives()const
{
	return mNumTimeDerivatives;
}

//! @brief set the beginning of the time increment to the structure
void NuTo::StructureBase::SetPrevTime(double rPrevTime)
{
	mPrevTime = rPrevTime;
}

//! @brief get the beginning of the time increment of the structure
double NuTo::StructureBase::GetPrevTime() const
{
	return mPrevTime;
}

//! @brief set the end of the time increment to the structure (current time)
void NuTo::StructureBase::SetTime(double rTime)
{
	mTime = rTime;
}

//! @brief get the end of the time increment of the structure (current time)
double NuTo::StructureBase::GetTime() const
{
	return mTime;
}

//! @brief set number of cycles to be extrapolated
void NuTo::StructureBase::SetNumExtrapolatedCycles(int rNumber)
{
	mNumExtrapolatedCycles = rNumber;
}

//! @brief get the number of cycles to be extrapolated
int NuTo::StructureBase::GetNumExtrapolatedCycles() const
{
	return mNumExtrapolatedCycles;
}

// store all elements of a group in a vector
void NuTo::StructureBase::GetElementsByGroup(const Group<ElementBase>* rElementGroup, std::vector<const ElementBase*>& rElements) const
{
    Group<ElementBase>::const_iterator ElementIter = rElementGroup->begin();
    while (ElementIter != rElementGroup->end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a group in a vector
void NuTo::StructureBase::GetElementsByGroup(Group<ElementBase>* rElementGroup, std::vector< ElementBase*>& rElements)
{
    Group<ElementBase>::iterator ElementIter = rElementGroup->begin();
    while (ElementIter != rElementGroup->end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}


// add visualization components for an element group
void NuTo::StructureBase::AddVisualizationComponent(int rElementGroup, VisualizeBase::eVisualizeWhat rVisualizeComponent)
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Element group does not exist.");

    // create a new visualization list for an element group or add components to an already existing list
    if (mGroupVisualizeComponentsMap.find(rElementGroup) == mGroupVisualizeComponentsMap.end())
    {
        std::list<std::shared_ptr<VisualizeComponent>> visualizationPtrList;
        visualizationPtrList.push_back(std::make_shared<VisualizeComponent>(VisualizeComponent(rVisualizeComponent)));

        mGroupVisualizeComponentsMap.insert(std::pair<int,std::list<std::shared_ptr<VisualizeComponent>>>(rElementGroup, visualizationPtrList));
        // mGroupVisualizeComponentsMap.emplace(rElementGroup, visualizationPtrList);       //<- use this for gcc version 4.9 or higher!

        mGroupVisualizationType.insert(std::pair<int, VisualizeBase::eVisualizationType>(rElementGroup, VisualizeBase::VORONOI_CELL));
        // mGroupVisualizationType.emplace(rElementGroup, VisualizeBase::VORONOI_CELL);     //<- use this for gcc version 4.9 or higher!
    } else
    {
        mGroupVisualizeComponentsMap.at(rElementGroup).push_back(std::make_shared<VisualizeComponent>(VisualizeComponent(rVisualizeComponent)));
    }

    for (auto const &iPair : mGroupVisualizeComponentsMap)
    {
        std::cout << "ele group: \t" << iPair.first << std::endl;
        for (auto const &iComponentPtr : iPair.second)
        {
            std::cout << "components: \t " << iComponentPtr->GetComponentName() << std::endl;
        }
    }

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<< __PRETTY_FUNCTION__ << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}


// add visualization components for an element group
void NuTo::StructureBase::AddVisualizationComponent(int rElementGroup, const std::string& rVisualizeComponent)
{

    if (rVisualizeComponent == "Accelerations")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ACCELERATION);
    else if (rVisualizeComponent.compare("AngularAccelerations"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ANGULAR_ACCELERATION);
    else if (rVisualizeComponent.compare("AngularVelocities"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ANGULAR_VELOCITY);
    else if (rVisualizeComponent.compare("BondStress"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::BOND_STRESS);
    else if (rVisualizeComponent.compare("BondStress"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::CONSTITUTIVE);
    else if (rVisualizeComponent.compare("Crack"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::CRACK);
    else if (rVisualizeComponent.compare("Damage"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::DAMAGE);
    else if (rVisualizeComponent.compare("Displacements"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::DISPLACEMENTS);
    else if (rVisualizeComponent.compare("Element"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ELEMENT);
    else if (rVisualizeComponent.compare("EngineeringPlasticStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_PLASTIC_STRAIN);
    else if (rVisualizeComponent.compare("EngineeringStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_STRAIN);
    else if (rVisualizeComponent.compare("EngineeringStress"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_STRESS);
    else if (rVisualizeComponent.compare("HeatFlux"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::HEAT_FLUX);
    else if (rVisualizeComponent.compare("LatticeStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::LATTICE_STRAIN);
    else if (rVisualizeComponent.compare("LatticeStress"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::LATTICE_STRESS);
    else if (rVisualizeComponent.compare("LocalEqStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::LOCAL_EQ_STRAIN);
    else if (rVisualizeComponent.compare("NonlocalEqStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::NONLOCAL_EQ_STRAIN);
    else if (rVisualizeComponent.compare("ParticleRadius"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::PARTICLE_RADIUS);
    else if (rVisualizeComponent.compare("PrincipalEngineeringStress"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);
    else if (rVisualizeComponent.compare("RelativeHumidity"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::RELATIVE_HUMIDITY);
    else if (rVisualizeComponent.compare("Rotations"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::ROTATION);
    else if (rVisualizeComponent.compare("Section"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::SECTION);
    else if (rVisualizeComponent.compare("Slip"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::SLIP);
    else if (rVisualizeComponent.compare("Temperature"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::TEMPERATURE);
    else if (rVisualizeComponent.compare("TotalInelasticEqStrain"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN);
    else if (rVisualizeComponent.compare("Velocities"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::VELOCITY);
    else if (rVisualizeComponent.compare("WaterVolumeFraction"))
        AddVisualizationComponent(rElementGroup, VisualizeBase::WATER_VOLUME_FRACTION);
    else
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Visualization component not implemented or misspelled.");



}
void NuTo::StructureBase::AddVisualizationComponentNonlocalWeights(int rElementGroup, int rElementId, int rIp)
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Element group does not exist.");

    const ElementBase *elementBase = ElementGetElementPtr(rElementId);
    int numIp = elementBase->GetNumIntegrationPoints();

    if (rIp < 0 or rIp >= numIp)
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Integration point number is out of range.");

    try
    {
        // create a new visualization list for an element group or add components to an already existing list
        if (mGroupVisualizeComponentsMap.find(rElementGroup) == mGroupVisualizeComponentsMap.end())
        {
            std::list<std::shared_ptr<VisualizeComponent>> visualizationPtrList;
            visualizationPtrList.push_back(std::make_shared<VisualizeComponentNonlocalWeight>(VisualizeComponentNonlocalWeight(elementBase, rElementId, rIp)));

            mGroupVisualizeComponentsMap.insert(std::pair<int,std::list<std::shared_ptr<VisualizeComponent>>>(rElementGroup, visualizationPtrList));
            // mGroupVisualizeComponentsMap.emplace(rElementGroup, visualizationPtrList);       //<- use this for gcc version 4.9 or higher!

            mGroupVisualizationType.insert(std::pair<int, VisualizeBase::eVisualizationType>(rElementGroup, VisualizeBase::VORONOI_CELL));
            // mGroupVisualizationType.emplace(rElementGroup, VisualizeBase::VORONOI_CELL);     //<- use this for gcc version 4.9 or higher!


        } else
        {
            mGroupVisualizeComponentsMap.at(rElementGroup).push_back(std::make_shared<VisualizeComponentNonlocalWeight>(VisualizeComponentNonlocalWeight(elementBase, rElementId, rIp)));
        }

    }
    catch (NuTo::MechanicsException &e)
     {
        e.AddMessage(std::string(__PRETTY_FUNCTION__) + "\t: error setting element and local ip number.");
        throw e;
     }
    catch(...)
     {
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: error setting element and local ip number.");
     }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<< __PRETTY_FUNCTION__ << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}


void NuTo::StructureBase::SetVisualizationType(const int rElementGroup, const VisualizeBase::eVisualizationType rVisualizationType)
{
    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Element group does not exist.");

    // check if the element group exists
    if (mGroupVisualizationType.find(rElementGroup) == mGroupVisualizationType.end())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Please add a visualization component first before setting the visualization type.");

    mGroupVisualizationType.at(rElementGroup) = rVisualizationType;
}


void NuTo::StructureBase::ClearVisualizationComponents()
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
//    mVisualizeComponents.clear();
    mGroupVisualizeComponentsMap.clear();
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ClearVisualizationComponents] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ExportVtkDataFileNodes(const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif

    for (auto const &it : mGroupVisualizeComponentsMap)
    {
        VisualizeUnstructuredGrid visualize;
        this->DefineVisualizeNodeData(visualize, it.second);
        this->NodeTotalAddToVisualize(visualize, it.second);

        if (rXML)
            visualize.ExportVtuDataFile(rResultFileName);
        else
            visualize.ExportVtkDataFile(rResultFileName);
    }

#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        mLogger << "[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}



void NuTo::StructureBase::ExportVtkDataFileElements(const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif

    for (auto const &it : mGroupVisualizeComponentsMap)
    {
        VisualizeUnstructuredGrid visualize;

        this->DefineVisualizeElementData(visualize, it.second);
        ElementGroupAddToVisualize(it.first, visualize, it.second);

        if (rXML)
            visualize.ExportVtuDataFile(rResultFileName);
        else
            visualize.ExportVtkDataFile(rResultFileName);
    }




    #ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif


    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeElementData(visualize,mGroupVisualizeComponentsMap.at(rGroupIdent));
    this->ElementGroupAddToVisualize(rGroupIdent,visualize,mGroupVisualizeComponentsMap.at(rGroupIdent));

    if (rXML)
        visualize.ExportVtuDataFile(rResultFileName);
    else
        visualize.ExportVtkDataFile(rResultFileName);


#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementGroupExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif // ENABLE_VISUALIZE
}

std::map<int, std::list<std::shared_ptr<NuTo::VisualizeComponent>>>& NuTo::StructureBase::GetGroupVisualizeComponentsMap(void)
{
#ifdef ENABLE_VISUALIZE
    return mGroupVisualizeComponentsMap;
#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::DefineVisualizeElementData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)const
{
#ifdef ENABLE_VISUALIZE

    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {

        case NuTo::VisualizeBase::SECTION:
        case NuTo::VisualizeBase::CONSTITUTIVE:
        case NuTo::VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN:
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
        case NuTo::VisualizeBase::LOCAL_EQ_STRAIN:
        case NuTo::VisualizeBase::ELEMENT:
        case NuTo::VisualizeBase::DAMAGE:
            rVisualize.DefineCellDataScalar(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::SLIP:
        case NuTo::VisualizeBase::CRACK:
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
        case NuTo::VisualizeBase::LATTICE_STRESS:
        case NuTo::VisualizeBase::LATTICE_STRAIN:
        case NuTo::VisualizeBase::LATTICE_PLASTIC_STRAIN:
            rVisualize.DefineCellDataVector(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
        case NuTo::VisualizeBase::BOND_STRESS:
            rVisualize.DefineCellDataTensor(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
        case NuTo::VisualizeBase::TEMPERATURE:
            rVisualize.DefinePointDataScalar(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::DISPLACEMENTS:
//        case NuTo::VisualizeBase::ENGINEERING_STRAIN: // this is a test
        case NuTo::VisualizeBase::VELOCITY:
        case NuTo::VisualizeBase::ACCELERATION:
            rVisualize.DefinePointDataVector(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::PARTICLE_RADIUS:
        case NuTo::VisualizeBase::ROTATION:
        case NuTo::VisualizeBase::ANGULAR_VELOCITY:
        case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
            //do nothing;
            break;

        default:
        	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: undefined visualize components.");
        }

    }
#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::DefineVisualizeNodeData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)const
{
#ifdef ENABLE_VISUALIZE

    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
        case NuTo::VisualizeBase::ROTATION:
        case NuTo::VisualizeBase::VELOCITY:
        case NuTo::VisualizeBase::ACCELERATION:
        case NuTo::VisualizeBase::ANGULAR_VELOCITY:
        case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
            rVisualize.DefinePointDataVector(it.get()->GetComponentName());
            break;
        case NuTo::VisualizeBase::PARTICLE_RADIUS:
        case NuTo::VisualizeBase::TEMPERATURE:
        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
            rVisualize.DefinePointDataScalar(it.get()->GetComponentName());
            break;
        default:
            // do nothing for integration point data in the visualization list. However, the visualization of new dofs needs to be added here!
        	break;
         }

    }
#endif // ENABLE_VISUALIZE
}





//! @brief ... evaluates the structur
void NuTo::StructureBase::Evaluate(std::map<StructureEnum::eOutput, StructureOutputBase *> &rStructureOutput)
{
    throw MechanicsException("[NuTo::StructureBase::Evaluate] Not implemented.");
}

void NuTo::StructureBase::BuildGlobalCoefficientMatrixCheck()
{
    // build global dof numbering if required
    if (this->mNodeNumberingRequired)
    {
        try
        {
            this->NodeBuildGlobalDofs();
        }
        catch (MechanicsException& e)
        {
            e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrixCheck] error building global dof numbering.");
            throw e;
        }
    }

    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientMatrixCheck] First update of tmp static data required.");
    }
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType rType, SparseMatrixCSRGeneral<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullVector<double,Eigen::Dynamic> activeDofValues;
    FullVector<double,Eigen::Dynamic> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(0,activeDofValues, dependentDofValues);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error extracting dof values from node.");
        throw e;
    }

//    rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    // resize output objects
    if (rMatrix.GetNumColumns()!=this->mNumActiveDofs || rMatrix.GetNumRows()!=this->mNumActiveDofs)
    {
        rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    }
    else
    {
        rMatrix.SetZeroEntries();
    }

    rVector.Resize(this->mNumActiveDofs);
    if (this->mConstraintMatrix.GetNumEntries() == 0)
    {
        //mLogger << "non-symmetric, zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesGeneral(rType, rMatrix, coefficientMatrixJK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build equivalent load vector
        rVector = coefficientMatrixJK * NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(dependentDofValues - this->mConstraintRHS);
    }
    else
    {
        //mLogger << "non-symmetric, non-zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRGeneral<double> coefficientMatrixKJ(this->mNumDofs - this->mNumActiveDofs, this->mNumActiveDofs);
        SparseMatrixCSRGeneral<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesGeneral(rType, rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build global matrix
        SparseMatrixCSRGeneral<double> transConstraintMatrix = this->mConstraintMatrix.Transpose();
        rMatrix -= transConstraintMatrix * coefficientMatrixKJ + coefficientMatrixJK * this->mConstraintMatrix;
        rMatrix += transConstraintMatrix * coefficientMatrixKK * this->mConstraintMatrix;

        // build equivalent load vector
        rVector = (transConstraintMatrix * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues);
    }
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType rType, SparseMatrixCSRSymmetric<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullVector<double,Eigen::Dynamic> activeDofValues;
    FullVector<double,Eigen::Dynamic> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(0,activeDofValues, dependentDofValues);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error extracting dof values from node.");
        throw e;
    }

    // resize output objects
    // resize output objects
    if (rMatrix.GetNumColumns()!=this->mNumActiveDofs || rMatrix.GetNumRows()!=this->mNumActiveDofs)
    {
        rMatrix.Resize(this->mNumActiveDofs);
    }
    else
    {
        rMatrix.SetZeroEntries();
    }
    rVector.Resize(this->mNumActiveDofs);
    if (this->mConstraintMatrix.GetNumEntries() == 0)
    {
        //mLogger << "symmetric, zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesSymmetric(rType, rMatrix, coefficientMatrixJK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
    }
    else
    {
        //mLogger << "symmetric, non-zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRSymmetric<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesSymmetric(rType, rMatrix, coefficientMatrixJK, coefficientMatrixKK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build global matrix
        rMatrix.Sub_TransA_Mult_TransB_Plus_B_Mult_A(this->mConstraintMatrix, coefficientMatrixJK);
        rMatrix.Add_TransA_Mult_B_Mult_A(this->mConstraintMatrix, coefficientMatrixKK);

        // build equivalent load vector
        FullVector<double,Eigen::Dynamic> deltaRHS = this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues;
        FullVector<double,Eigen::Dynamic> Kdd_Mult_DeltaRHS = coefficientMatrixKK * deltaRHS;
        rVector = this->mConstraintMatrix.TransMult(Kdd_Mult_DeltaRHS) - coefficientMatrixJK * deltaRHS;
    }
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType rType, SparseMatrixCSRVector2General<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullVector<double,Eigen::Dynamic> activeDofValues;
    FullVector<double,Eigen::Dynamic> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(0,activeDofValues, dependentDofValues);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error extracting dof values from node.");
        throw e;
    }

//    rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    // resize output objects
    if (rMatrix.GetNumColumns()!=this->mNumActiveDofs || rMatrix.GetNumRows()!=this->mNumActiveDofs)
    {
        rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    }
    else
    {
        rMatrix.SetZeroEntries();
    }

    rVector.Resize(this->mNumActiveDofs);
    if (this->mConstraintMatrix.GetNumEntries() == 0)
    {
        //mLogger << "non-symmetric, zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesGeneral(rType, rMatrix, coefficientMatrixJK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
        /*        mLogger << "dependent " << "\n";
                dependentDofValues.Trans().Info(12,10);
                mLogger << "mConstraintRHS " << "\n";
                mConstraintRHS.Trans().Info(12,10);
                mLogger << "delta " << "\n";
                (dependentDofValues - this->mConstraintRHS).Trans().Info(12,10);
                mLogger << "coefficientMatrixJK" << "\n";
                (NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(coefficientMatrixJK)).Info(12,3);
        */
    }
    else
    {
        //mLogger << "non-symmetric, non-zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRVector2General<double> coefficientMatrixKJ(this->mNumDofs - this->mNumActiveDofs, this->mNumActiveDofs);
        SparseMatrixCSRVector2General<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesGeneral(rType, rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build global matrix
        SparseMatrixCSRVector2General<double> constraintMatrixVector2 (this->mConstraintMatrix);
        SparseMatrixCSRVector2General<double> transConstraintMatrixVector2 (constraintMatrixVector2.Transpose());
        rMatrix -= transConstraintMatrixVector2 * coefficientMatrixKJ + coefficientMatrixJK * constraintMatrixVector2;
        rMatrix += transConstraintMatrixVector2 * coefficientMatrixKK * constraintMatrixVector2;

        // build equivalent load vector
        rVector = (transConstraintMatrixVector2 * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - constraintMatrixVector2 * activeDofValues);
    }
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType rType, SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullVector<double,Eigen::Dynamic> activeDofValues;
    FullVector<double,Eigen::Dynamic> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(0,activeDofValues, dependentDofValues);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error extracting dof values from node.");
        throw e;
    }

//    rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    // resize output objects
    if (rMatrix.GetNumColumns()!=this->mNumActiveDofs || rMatrix.GetNumRows()!=this->mNumActiveDofs)
    {
        rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
    }
    else
    {
        rMatrix.SetZeroEntries();
    }

    rVector.Resize(this->mNumActiveDofs);
    if (this->mConstraintMatrix.GetNumEntries() == 0)
    {
        //mLogger << "non-symmetric, zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesSymmetric(rType, rMatrix, coefficientMatrixJK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
    }
    else
    {
        //mLogger << "non-symmetric, non-zero constraint matrix" << "\n";

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        //TODO this should actually be a symmetric matrix, but the multiplications are not implemented
        SparseMatrixCSRVector2General<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        Error::eError error = this->BuildGlobalCoefficientSubMatricesSymmetric(rType, rMatrix, coefficientMatrixJK, coefficientMatrixKK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build global matrix
        SparseMatrixCSRVector2General<double> constraintMatrixVector2 (this->mConstraintMatrix);
        SparseMatrixCSRVector2General<double> transConstraintMatrixVector2 (constraintMatrixVector2.Transpose());
        rMatrix -= (coefficientMatrixJK * constraintMatrixVector2).SymmetricPart();

        //this should be optimized to actually be performed as C^TKC returning a symmetric matrix - just don't have time right now to do that
        rMatrix += (transConstraintMatrixVector2 * coefficientMatrixKK * constraintMatrixVector2).SymmetricPart();

        // build equivalent load vector
        rVector = (transConstraintMatrixVector2 * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - constraintMatrixVector2 * activeDofValues);
    }
    return Error::SUCCESSFUL;
}


//! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (nonsymmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::STIFFNESS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}

//! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (symmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::STIFFNESS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (nonsymmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::STIFFNESS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::STIFFNESS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (damping) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (symmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix1(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::DAMPING, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix1] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix1] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (nonsymmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::MASS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (symmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::MASS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (nonsymmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::MASS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}


//! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureEnum::eMatrixType::MASS, rMatrix, rVector);
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix2] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return error;
}






// build global external load vector
void NuTo::StructureBase::BuildGlobalExternalLoadVector(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // check dof numbering
    if (this->mNodeNumberingRequired)
    {
        try
        {
            this->NodeBuildGlobalDofs();
        }
        catch (MechanicsException& e)
        {
            e.AddMessage("[NuTo::StructureBase::BuildGlobalExternalLoadVector] error building global dof numbering.");
            throw e;
        }
    }

    rVector.Resize(this->mNumActiveDofs);
    FullVector<double,Eigen::Dynamic> dependentDofLoadVector(this->mNumDofs - this->mNumActiveDofs);

    // loop over all loads
    boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
    while (loadIter != this->mLoadMap.end())
    {
        loadIter->second->AddLoadToGlobalSubVectors(rLoadCase, rVector, dependentDofLoadVector);
        loadIter++;
    }
    if (this->mConstraintMatrix.GetNumEntries() != 0)
    {
        rVector -=  this->mConstraintMatrix.TransMult(dependentDofLoadVector);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalExternalLoadVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

// build global external load vector
void NuTo::StructureBase::BuildGlobalExternalLoadVector(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rVector_j, NuTo::FullVector<double,Eigen::Dynamic>& rVector_k)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // check dof numbering
    if (this->mNodeNumberingRequired)
    {
        try
        {
            this->NodeBuildGlobalDofs();
        }
        catch (MechanicsException& e)
        {
            e.AddMessage("[NuTo::StructureBase::BuildGlobalExternalLoadVector] error building global dof numbering.");
            throw e;
        }
    }

    rVector_j.Resize(this->mNumActiveDofs);
    rVector_k.Resize(this->mNumDofs - this->mNumActiveDofs);

    // loop over all loads
    boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
    while (loadIter != this->mLoadMap.end())
    {
        loadIter->second->AddLoadToGlobalSubVectors(rLoadCase, rVector_j, rVector_k);
        loadIter++;
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalExternalLoadVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


// build global gradient of the internal potential (e.g. the internal forces)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector(NuTo::FullVector<double,Eigen::Dynamic>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    // check dof numbering
    if (this->mNodeNumberingRequired)
    {
        try
        {
            this->NodeBuildGlobalDofs();
        }
        catch (MechanicsException& e)
        {
            e.AddMessage("[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] error building global dof numbering.");
            throw e;
        }
    }
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] First update of tmp static data required.");
    }


    rVector.Resize(this->mNumActiveDofs);
    FullVector<double,Eigen::Dynamic> dependentDofGradientVector(this->mNumDofs - this->mNumActiveDofs);

    try
    {
        // build sub vectors, update of history variables afterwards=false
        Error::eError error = this->BuildGlobalGradientInternalPotentialSubVectors(rVector, dependentDofGradientVector, false);
        if (error!=Error::SUCCESSFUL)
            return error;
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] error building sub vectors.");
        throw e;
    }
    if (this->mConstraintMatrix.GetNumEntries() != 0)
    {
        rVector -=  this->mConstraintMatrix.TransMult(dependentDofGradientVector);
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return Error::SUCCESSFUL;
}

//! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
//! values smaller than that one will not be added to the global matrix
void NuTo::StructureBase::SetToleranceStiffnessEntries(double rToleranceStiffnessEntries)
{
    mToleranceStiffnessEntries = rToleranceStiffnessEntries;
}

//! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
//! values smaller than that one will not be added to the global matrix
double NuTo::StructureBase::GetToleranceStiffnessEntries()const
{
    return mToleranceStiffnessEntries;
}

//! @brief returns the number of degrees of freedom
//! @return ... number of degrees of freedom
int NuTo::StructureBase::GetNumDofs()const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::GetNumDofs] Build global Dofs first.");
    }
    return mNumDofs;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::GetNumDofs] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief returns the number of active degrees of freedom
//! @return ... number of active degrees of freedom
int NuTo::StructureBase::GetNumActiveDofs()const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::GetNumActiveDofs] Build global Dofs first.");
    }
    return mNumActiveDofs;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::GetNumActiveDofs] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief returns the number of dependent degrees of freedom
//! @return ... number of active degrees of freedom
int NuTo::StructureBase::GetNumDependentDofs()const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::GetNumDependentDofs] Build global Dofs first.");
    }
    return mNumDofs-mNumActiveDofs;
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::GetNumDependentDofs] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief returns the a reference to the constraint matrix
const NuTo::SparseMatrixCSRGeneral<double>& NuTo::StructureBase::GetConstraintMatrix()const
{
    return mConstraintMatrix;
}

#ifdef _OPENMP
//@brief determines the maximum independent sets and stores it at the structure
void NuTo::StructureBase::CalculateMaximumIndependentSets()
{
#define UNDONE 1
#define SELECTED 2
#define DELETED 3
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        std::vector<ElementBase*> elementVector;
        GetElementsTotal(elementVector);

        /*    //for test purpose
        	mMIS.resize(elementVector.size());
        	for (unsigned int elementCount = 0; elementCount< mMIS.size(); elementCount++)
        	{
        		mMIS[elementCount].resize(1);
        		mMIS[elementCount][0] = elementVector[elementCount];
        	}
        	return;
        	*/

        //Build the connectivity graph
        //First get for all nodes all the elements
        std::map<const NodeBase*, std::vector<unsigned int> > elementsPerNode;
        for (unsigned int elementCount = 0; elementCount< elementVector.size(); elementCount++)
        {
            for (int nodeCount = 0; nodeCount< elementVector[elementCount]->GetNumInfluenceNodes(); nodeCount++)
            {
                elementsPerNode[elementVector[elementCount]->GetInfluenceNode(nodeCount)].push_back(elementCount);
            }
        }

        //Get the neighboring elements (always referring to the location in the vector elementVector)
        std::vector< std::vector<int> > NeighborElements(elementVector.size());
        for (auto itNode = elementsPerNode.begin(); itNode!=elementsPerNode.end(); itNode++)
        {
            for (unsigned int elementCount1 = 0; elementCount1 < itNode->second.size(); elementCount1++)
            {
                for (unsigned int elementCount2 = elementCount1+1; elementCount2 < itNode->second.size(); elementCount2++)
                {
                    NeighborElements[itNode->second[elementCount1]].push_back(itNode->second[elementCount2]);
                    NeighborElements[itNode->second[elementCount2]].push_back(itNode->second[elementCount1]);
                }
            }
        }

        //build the maximum independent sets
        std::vector<int> elementState(elementVector.size());
        for (unsigned int count=0; count<elementState.size(); count++)
            elementState[count] = UNDONE;

        unsigned int numDeleted=0;
        unsigned int curMIS = 0;
        mMIS.resize(10);
        while (numDeleted<elementVector.size())
        {
            if (mMIS.size()<=curMIS)
                mMIS.resize(curMIS+1);
            for (unsigned int countElement=0; countElement<elementVector.size(); countElement++)
            {
                if (elementState[countElement]!=UNDONE)
                    continue;

                //add element to the set
                //std::cout << "add element " << ElementGetId(elementVector[countElement]) << " to set " << curMIS << std::endl;
                (mMIS[curMIS]).push_back(elementVector[countElement]);
                elementState[countElement] = DELETED;
                numDeleted++;

                //mark all the neighboring elements as selected, which prevents them to being added to this set
                for (unsigned int theNeighborCount=0; theNeighborCount<NeighborElements[countElement].size(); theNeighborCount++)
                {
                    if (elementState[NeighborElements[countElement][theNeighborCount]]==UNDONE)
                        elementState[NeighborElements[countElement][theNeighborCount]] = SELECTED;
                }
            }
            //reset the selected elements to be undone
            for (unsigned int countElement=0; countElement<elementVector.size(); countElement++)
            {
                if (elementState[countElement]==SELECTED)
                    elementState[countElement]=UNDONE;

            }
            curMIS++;
        }
        mMIS.resize(curMIS);
        /*    std::cout << "maximum number of independent sets " << mMIS.size() << std::endl;
            for (unsigned int count=0; count<mMIS.size(); count++)
            {
            	std::cout << "MIS " << count << " with " << mMIS[count].size() << " elements " << std::endl;
            	for (unsigned int count2=0 ; count2<mMIS[count].size(); count2++)
            		std::cout << ElementGetId(mMIS[count][count2]) << " ";
            	std::cout << std::endl;
            }
        */
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::CalculateMaximumIndependentSets] error calculating maximum independent sets.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::CalculateMaximumIndependentSets] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}
#else
//@brief determines the maximum independent sets and stores it at the structure, do nothing for applications without openmp
void NuTo::StructureBase::CalculateMaximumIndependentSets()
{
}
#endif


//! @brief is only true for structure used as multiscale (structure in a structure)
//! @parameters rTypeOfSpecimen 0 box, 1 dogbone
//! @parameters rBoundingBox box for the spheres (3*2 matrix)
//! @parameters rSeed seed for the random number generator
//! @parameters rRadiusBoundaryParticles radius particles simulated on the boundary
//! @parameters rDistanceBoundaryParticles distance of the boundary particles
//! @parameters rTypeOfSpecimen 0 box, 1 dogbone
//! @return ... matrix with spheres (coordinates x y z and radius)
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::CreateSpheresOnSpecimenBoundary(int rTypeOfSpecimen,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, int rSeed,
		double rRadiusBoundaryParticles, double rDistanceBoundaryParticles)
{
    if (rBoundingBox.GetNumRows()!=3 && rBoundingBox.GetNumColumns()!=2)
    	throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] bounding box has to have the dimension [3,2]");

	// calculate specimen length
    std::array<double,3> lBox;
    for (int count=0; count<3; count++)
    {
    	lBox[count] = rBoundingBox(count,1) - rBoundingBox(count,0);
    	if (lBox[count]<=0)
    		throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] box dimensions should be not negative.");
    }

	FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> particles;
	particles.Resize(1000,4);

	int numParticles(0);

	// random number generator
	dsfmt_t randomNumberGenerator;
	// init random number generator with milliseconds from ..
	dsfmt_init_gen_rand(&randomNumberGenerator, rSeed);

	for (int count=0; count<3; count++)
		if (rBoundingBox(0,1)-rBoundingBox(0,0)<10.*rRadiusBoundaryParticles)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] minimum diameter and box dimension do not match.");

	std::vector<boost::array<double,3> > nodes; //node coordinates
	std::vector<boost::array<double,3> > helpCenters; //coordinates of circumcenter of circular edges
    std::vector<boost::array<int   ,3> > edges; //edges that connect the nodes (third entry is the center of the circumcenter in the xy plane (helpCenters), otherwise -1
    std::vector<boost::array<int   ,4> > faces; //faces with 2 edges (only regular quads with opposite sides being parallel
                                                //implemented, one edge can be a circle
                                                //third and forth value give the common point of both edges (either 0 or zero)

	switch (rTypeOfSpecimen)
	{
	case 0:
	{
		//box
	    switch (mDimension)
		{
		case 3:
			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),rBoundingBox(2,1) }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,-1}});
			edges.push_back(boost::array<int, 3> {{1,3,-1}});
			edges.push_back(boost::array<int, 3> {{3,2,-1}});
			edges.push_back(boost::array<int, 3> {{2,0,-1}});
			edges.push_back(boost::array<int, 3> {{4,5,-1}});
			edges.push_back(boost::array<int, 3> {{5,7,-1}});
			edges.push_back(boost::array<int, 3> {{7,6,-1}});
			edges.push_back(boost::array<int, 3> {{6,4,-1}});
			edges.push_back(boost::array<int, 3> {{0,4,-1}});
			edges.push_back(boost::array<int, 3> {{1,5,-1}});
			edges.push_back(boost::array<int, 3> {{3,7,-1}});
			edges.push_back(boost::array<int, 3> {{2,6,-1}});

			//add the faces
			faces.push_back(boost::array<int, 4> {{0,1,1,0}});
			faces.push_back(boost::array<int, 4> {{0,9,1,0}});
			faces.push_back(boost::array<int, 4> {{1,10,1,0}});
			faces.push_back(boost::array<int, 4> {{2,11,1,0}});
			faces.push_back(boost::array<int, 4> {{3,8,1,0}});
			faces.push_back(boost::array<int, 4> {{4,5,1,0}});
			break;
		case 2:
		{
			double rCuttingPlane(0.5*(rBoundingBox(2,0)+rBoundingBox(2,1)));

			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),rCuttingPlane }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,-1}});
			edges.push_back(boost::array<int, 3> {{1,3,-1}});
			edges.push_back(boost::array<int, 3> {{3,2,-1}});
			edges.push_back(boost::array<int, 3> {{2,0,-1}});

		}
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] only implemented for 2D and 3D.");
		}
	}
	break;
	case 1:
	{
		double D = rBoundingBox(0,1) - rBoundingBox(1,0);
		if (fabs(rBoundingBox.at(1,1) - rBoundingBox.at(1,0)-1.5*D)>1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");
		//dogbone specimen
		switch (mDimension)
		{
		case 3:
			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1)-0.25*D,rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.25*D,rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.25*D,rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1)-0.25*D,rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),rBoundingBox(2,0) }});


			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1)-0.25*D,rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.25*D,rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.25*D,rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1)-0.25*D,rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),rBoundingBox(2,1) }});

			//add the help centers of the triangles
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)-0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rBoundingBox(2,0) }});
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,1)+0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rBoundingBox(2,0) }});
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)-0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rBoundingBox(2,1) }});
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,1)+0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rBoundingBox(2,1) }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,-1}});
			edges.push_back(boost::array<int, 3> {{2,1,0}});
			edges.push_back(boost::array<int, 3> {{2,3,-1}});
			edges.push_back(boost::array<int, 3> {{3,4,-1}});
			edges.push_back(boost::array<int, 3> {{4,5,-1}});
			edges.push_back(boost::array<int, 3> {{6,5,1}});
			edges.push_back(boost::array<int, 3> {{6,7,-1}});
			edges.push_back(boost::array<int, 3> {{7,0,-1}});

			edges.push_back(boost::array<int, 3> {{8,9,-1}});
			edges.push_back(boost::array<int, 3> {{10,9,2}});
			edges.push_back(boost::array<int, 3> {{10,11,-1}});
			edges.push_back(boost::array<int, 3> {{11,12,-1}});
			edges.push_back(boost::array<int, 3> {{12,13,-1}});
			edges.push_back(boost::array<int, 3> {{14,13,3}});
			edges.push_back(boost::array<int, 3> {{14,15,-1}});
			edges.push_back(boost::array<int, 3> {{15,8,-1}});

			edges.push_back(boost::array<int, 3> {{0,8,-1}});
			edges.push_back(boost::array<int, 3> {{1,9,-1}});
			edges.push_back(boost::array<int, 3> {{2,10,-1}});
			edges.push_back(boost::array<int, 3> {{3,11,-1}});
			edges.push_back(boost::array<int, 3> {{4,12,-1}});
			edges.push_back(boost::array<int, 3> {{5,13,-1}});
			edges.push_back(boost::array<int, 3> {{6,14,-1}});
			edges.push_back(boost::array<int, 3> {{7,15,-1}});

			//add the faces (only the vertical ones, the horizontal ones are tricky and done separately
			faces.push_back(boost::array<int, 4> {{0,16,0,0}});
			faces.push_back(boost::array<int, 4> {{1,17,1,0}});
			faces.push_back(boost::array<int, 4> {{2,18,0,0}});
			faces.push_back(boost::array<int, 4> {{3,19,0,0}});
			faces.push_back(boost::array<int, 4> {{4,20,0,0}});
			faces.push_back(boost::array<int, 4> {{5,21,1,0}});
			faces.push_back(boost::array<int, 4> {{6,22,0,0}});
			faces.push_back(boost::array<int, 4> {{7,23,0,0}});

			break;
		case 2:
		{
			double rCuttingPlane(0.5*(rBoundingBox(2,0)+rBoundingBox(2,1)));

			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1),      rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,1)-0.25*D,rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.25*D,rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0),      rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0),      rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.25*D,rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1)-0.25*D,rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,1),      rCuttingPlane }});

			//add the help centers of the triangles
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)-0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rCuttingPlane }});
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,1)+0.525*D,0.5*(rBoundingBox(1,0)+rBoundingBox(1,1)),rCuttingPlane }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,-1}});
			edges.push_back(boost::array<int, 3> {{2,1,0}});
			edges.push_back(boost::array<int, 3> {{2,3,-1}});
			edges.push_back(boost::array<int, 3> {{3,4,-1}});
			edges.push_back(boost::array<int, 3> {{4,5,-1}});
			edges.push_back(boost::array<int, 3> {{6,5,1}});
			edges.push_back(boost::array<int, 3> {{6,7,-1}});
			edges.push_back(boost::array<int, 3> {{7,0,-1}});
		}
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] only implemented for 2D and 3D.");
		}
	}
	break;
	case 2:
	{
		//brazilian splitting (cylinder)
		double D = rBoundingBox(0,1) - rBoundingBox(1,0);
		if (fabs(rBoundingBox.at(1,1) - rBoundingBox.at(1,0)-D)>1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] for the cylinder specimen, the y dimension should be the same as the x-dimension.");
		//dogbone specimen
		switch (mDimension)
		{
		case 3:
			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.5*D,rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,1),rBoundingBox(2,0) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.5*D,rBoundingBox(2,0) }});

			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.5*D,rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,1),rBoundingBox(2,1) }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.5*D,rBoundingBox(2,1) }});

			//add the help centers of the triangles
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0)+0.5*D,rBoundingBox(2,0) }});
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0)+0.5*D,rBoundingBox(2,1) }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,0}});
			edges.push_back(boost::array<int, 3> {{1,2,0}});
			edges.push_back(boost::array<int, 3> {{2,3,0}});
			edges.push_back(boost::array<int, 3> {{3,0,0}});

			edges.push_back(boost::array<int, 3> {{4,5,1}});
			edges.push_back(boost::array<int, 3> {{5,6,1}});
			edges.push_back(boost::array<int, 3> {{6,7,1}});
			edges.push_back(boost::array<int, 3> {{7,4,1}});

			edges.push_back(boost::array<int, 3> {{0,4,-1}});
			edges.push_back(boost::array<int, 3> {{1,5,-1}});
			edges.push_back(boost::array<int, 3> {{2,6,-1}});
			edges.push_back(boost::array<int, 3> {{3,7,-1}});

			//add the faces (only the vertical ones, the horizontal ones are tricky and done separately
			faces.push_back(boost::array<int, 4> {{0,8,0,0}});
			faces.push_back(boost::array<int, 4> {{1,9,0,0}});
			faces.push_back(boost::array<int, 4> {{2,10,0,0}});
			faces.push_back(boost::array<int, 4> {{3,11,0,0}});

			break;
		case 2:
		{
			double rCuttingPlane(0.5*(rBoundingBox(2,0)+rBoundingBox(2,1)));

			//add the corners
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0),rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,1),rBoundingBox(1,0)+0.5*D,rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,1),rCuttingPlane }});
			nodes.push_back(boost::array<double, 3> {{rBoundingBox(0,0),rBoundingBox(1,0)+0.5*D,rCuttingPlane }});

			//add the help centers of the triangles
			helpCenters.push_back(boost::array<double, 3> {{rBoundingBox(0,0)+0.5*D,rBoundingBox(1,0)+0.5*D,rCuttingPlane }});

			//add the edge
			edges.push_back(boost::array<int, 3> {{0,1,0}});
			edges.push_back(boost::array<int, 3> {{1,2,0}});
			edges.push_back(boost::array<int, 3> {{2,3,0}});
			edges.push_back(boost::array<int, 3> {{3,0,0}});
		}
		break;
		default:
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] only implemented for 2D and 3D.");
		}
	}
	break;
	default:
	throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] unknown specimen type.");
	}

	// first add the corners
	for (unsigned int count=0; count<nodes.size(); count++)
	{
		// check size of table
		if(numParticles == particles.GetNumRows())
		{
			particles.ConservativeResizeRows(particles.GetNumRows()+1000);
		}
		particles(numParticles,0) = nodes[count][0];
	    particles(numParticles,1) = nodes[count][1];
    	particles(numParticles,2) = nodes[count][2];
		particles(numParticles,3) = rRadiusBoundaryParticles;
		numParticles++;
	}

	// add particles along the edges
	for (unsigned int theEdge=0; theEdge<edges.size(); theEdge++)
	{
		std::cout << "edge " << theEdge+1 << "/" << edges.size() << "\n";
		boost::array<double,3>& node1(nodes[edges[theEdge][0]]);
		boost::array<double,3>& node2(nodes[edges[theEdge][1]]);

		double lEdge;
		double deltaAngle(0),angle1(0),radius(0);
		if (edges[theEdge][2]!=-1)
		{
			//cirle
			boost::array<double,3>& node3(helpCenters[edges[theEdge][2]]);
			boost::array<double,3> delta1;
			boost::array<double,3> delta2;
			for (int count=0; count<3; count++)
			{
			    delta1[count] = node1[count]-node3[count];
			    delta2[count] = node2[count]-node3[count];
			}
			double r1 = sqrt(delta1[0]*delta1[0]+delta1[1]*delta1[1]);
			double r2 = sqrt(delta2[0]*delta2[0]+delta2[1]*delta2[1]);
			if (fabs(delta1[2])>1e-10 || fabs(delta2[2])>1e-10)
				throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] circular edges only in the xy plane.");
			if (fabs(r1-r2)>1e-10)
				throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] circumcenter of the circular edge is not correct.");
			radius = r1;
			angle1 = atan2(delta1[1],delta1[0]);
			double angle2 = atan2(delta2[1],delta2[0]);
			deltaAngle = angle2-angle1;
			while (deltaAngle>2.*M_PI)
				deltaAngle-=2.*M_PI;
			while (deltaAngle<0)
				deltaAngle+=2.*M_PI;
			lEdge = r1*deltaAngle;
		}
		else
		{
			//straight line
			boost::array<double,3> delta;
			for (int count=0; count<3; count++)
			    delta[count] = node2[count]-node1[count];

			lEdge = sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
		}

		int numParticlesPerEdge = lEdge/rDistanceBoundaryParticles;

		int thisParticlesPerEdge(2); //the ends have already been inserted
		int numTries(0);
		while(thisParticlesPerEdge<numParticlesPerEdge)
		{
			// check size of table
			if(numParticles == particles.GetNumRows())
			{
				particles.ConservativeResizeRows(particles.GetNumRows()+1000);
			}

			double s = dsfmt_genrand_close_open(&randomNumberGenerator);
			if (edges[theEdge][2]==-1)
			{
				//straight line
				boost::array<double,3> delta;
				for (int count=0; count<3; count++)
				    delta[count] = node2[count]-node1[count];
				//straight edge
				for (int count=0; count<3; count++)
				{
					particles(numParticles,count) = node1[count]+s*delta[count];
				}
			}
			else
			{
				//circle
				boost::array<double,3>& node3(helpCenters[edges[theEdge][2]]);
				double angle=angle1+s*deltaAngle;

				particles(numParticles,0) = node3[0]+radius*cos(angle);
				particles(numParticles,1) = node3[1]+radius*sin(angle);
				particles(numParticles,2) = node1[2];
			}

			particles(numParticles,3) = rRadiusBoundaryParticles;
			//check for overlap
			bool noSeparation(true);
			for (int countParticle=0; countParticle<numParticles; countParticle++)
			{
				double deltaX = particles(countParticle,0)-particles(numParticles,0);
				double deltaY = particles(countParticle,1)-particles(numParticles,1);
				double deltaZ = particles(countParticle,2)-particles(numParticles,2);
				double sumR = particles(countParticle,3)+particles(numParticles,3);
				if (deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ<sumR*sumR)
				{
					noSeparation = false;
					//std::cout << "failure at " << particles(numParticles,0 ) << " " << particles(numParticles,1 ) << "\n";
					break;
				}
			}
			if (noSeparation==true)
			{
				thisParticlesPerEdge++;
				numParticles++;
				numTries=0;
				//std::cout << "Particles " << "\n" << particles.mEigenMatrix.block(0,0,numParticles,4) << "\n";
			}
			else
			{
				numTries++;
				if (numTries>10000)
					throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] could not add particles on a boundary edge.");
			}
		}
	}

	// add particles along the faces
	for (unsigned int theFace=0; theFace<faces.size(); theFace++)
	{
		std::cout << "face " << theFace+1 << "/" << faces.size() << "\n";

		int theEdge[2];
		theEdge[0] = faces[theFace][0];
		theEdge[1] = faces[theFace][1];
		int endPoint[2];
		endPoint[0] = faces[theFace][2]; //first or second node of edge is the common one
		endPoint[1] = faces[theFace][3];

		if (edges[theEdge[0]][2]!=-1 && edges[theEdge[1]][2]!=-1)
		{
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] only one edge might be curved.");
		}

		if (edges[theEdge[0]][endPoint[0]]!=edges[theEdge[1]][endPoint[1]])
		{
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] the edges of the faces do not have a common end point.");
		}


		//define the nodes per face here (0,1 and 0,2 define the edge and 3 is at most one help node (there is only one)
		boost::array<double,3>& node1(nodes[edges[theEdge[0]][endPoint[0]]]);

		double lEdge[2];
		double radius(0), angle1(0), deltaAngle(0);
		boost::array<boost::array<double,3> ,2> delta;
		for (int countEdge=0; countEdge<2; countEdge++)
		{
			boost::array<double,3>& node2(nodes[edges[theEdge[countEdge]][endPoint[countEdge]==1 ? 0 : 1]]);
			if (edges[theEdge[countEdge]][2]!=-1)
			{
				boost::array<double,3>& node3(helpCenters[edges[theEdge[countEdge]][2]]);

				//circle
				boost::array<double,3> delta1;
				boost::array<double,3> delta2;
				for (int count=0; count<3; count++)
				{
				    delta1[count] = node1[count]-node3[count];
				    delta2[count] = node2[count]-node3[count];
				}
				double r1 = sqrt(delta1[0]*delta1[0]+delta1[1]*delta1[1]);
				double r2 = sqrt(delta2[0]*delta2[0]+delta2[1]*delta2[1]);
				if (fabs(delta1[2])>1e-10 || fabs(delta2[2])>1e-10)
					throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] circular edges only in the xy plane.");
				if (fabs(r1-r2)>1e-10)
					throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] circumcenter of the circular edge is not correct.");
				radius = r1;
				angle1 = atan2(delta1[1],delta1[0]);
				double angle2 = atan2(delta2[1],delta2[0]);
				deltaAngle = angle2-angle1;
				if (endPoint[countEdge]==1)
				{
					angle1 = angle2;
					deltaAngle*=-1;
				}
				while (deltaAngle>2.*M_PI)
					deltaAngle-=2.*M_PI;
				while (deltaAngle<0)
					deltaAngle+=2.*M_PI;
				lEdge[countEdge] = r1*deltaAngle;
			}
			else
			{
				//straight line
				for (int count=0; count<3; count++)
				    delta[countEdge][count] = node2[count]-node1[count];

				lEdge[countEdge] = sqrt(delta[countEdge][0]*delta[countEdge][0]+delta[countEdge][1]*delta[countEdge][1]+delta[countEdge][2]*delta[countEdge][2]);
			}
		}

		double lArea = lEdge[0]*lEdge[1];

		int numParticlesPerFace = lArea/(rDistanceBoundaryParticles*rDistanceBoundaryParticles);

		//the edges have already been inserted
		int thisParticlesPerFace(2*((int)(lEdge[0]/rDistanceBoundaryParticles)+(int)(lEdge[1]/rDistanceBoundaryParticles)-2));

		int numTries(0);
		while(thisParticlesPerFace<numParticlesPerFace)
		{
			// check size of table
			if(numParticles == particles.GetNumRows())
			{
				particles.ConservativeResizeRows(particles.GetNumRows()+1000);
			}

			particles(numParticles,0) = node1[0];
			particles(numParticles,1) = node1[1];
			particles(numParticles,2) = node1[2];
			for (int countEdge=0; countEdge<2; countEdge++)
			{
				double s = dsfmt_genrand_close_open(&randomNumberGenerator);
				if (edges[theEdge[countEdge]][2]!=-1)
				{
					//circle
					boost::array<double,3>& node3(helpCenters[edges[theEdge[countEdge]][2]]);
					double angle=angle1+s*deltaAngle;

					//std::cout << "particle location " << node3[0] << " " << radius <<  " " << angle <<  " " << node1[0] << "\n";
					particles(numParticles,0) += node3[0]+radius*cos(angle)-node1[0];
					particles(numParticles,1) += node3[1]+radius*sin(angle)-node1[1];
					//std::cout << "particle location " << particles(numParticles,0) << " " << particles(numParticles,1) <<  " " << particles(numParticles,2) << "\n";
				}
				else
				{
					//straight line
					for (int count=0; count<3; count++)
					{
						particles(numParticles,count) += s*delta[countEdge][count];
					}
				}
			}

			particles(numParticles,3) = rRadiusBoundaryParticles;
			//check for overlap
			bool noSeparation(true);
			for (int countParticle=0; countParticle<numParticles; countParticle++)
			{
				double deltaX = particles(countParticle,0)-particles(numParticles,0);
				double deltaY = particles(countParticle,1)-particles(numParticles,1);
				double deltaZ = particles(countParticle,2)-particles(numParticles,2);
				double sumR = particles(countParticle,3)+particles(numParticles,3);
				if (deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ<sumR*sumR)
				{
					noSeparation = false;
					//std::cout << "failure at " << particles(numParticles,0 ) << " " << particles(numParticles,1 ) << "\n";
					break;
				}
			}
			if (noSeparation==true)
			{
				thisParticlesPerFace++;
				numParticles++;
				numTries=0;
				//std::cout << "Particles " << "\n" << particles.mEigenMatrix.block(0,0,numParticles,4) << "\n";
			}
			else
			{
				numTries++;
				if (numTries>10000)
					throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] could not add particles on a boundary face.");
			}
		}
	}

	// add particles along special faces not included in the general set up
	switch(rTypeOfSpecimen)
	{
	case 0:
		//box - nothing to be done
		break;
	case 1:
		//dogbone - add the top and bottom surface
		for (unsigned int theFace=0; theFace<2; theFace++)
		{
			double D = rBoundingBox(0,1) - rBoundingBox(0,0);
			double lArea = 1.5*D*D;
			if (fabs(rBoundingBox.at(1,1) - rBoundingBox.at(1,0)-1.5*D)>1e-10)
				throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");
			//subtract the circles
			double radius = 0.725*D; ;
			double deltaAngle = 2.*0.2/0.525;
			lArea -= 2.*(deltaAngle/(2.*M_PI)*M_PI*radius*radius);

			int numParticlesPerFace = lArea/(rDistanceBoundaryParticles*rDistanceBoundaryParticles);

			//the edges have already been inserted
			int thisParticlesPerFace(2*((int)(D/rDistanceBoundaryParticles)+
					                    (int)(0.2*D/rDistanceBoundaryParticles)+
					                    (int)(0.2*D/rDistanceBoundaryParticles)+
					                    (int)(deltaAngle*radius/rDistanceBoundaryParticles))-8);

			int numTries(0);
			while(thisParticlesPerFace<numParticlesPerFace)
			{
				// check size of table
				if(numParticles == particles.GetNumRows())
				{
					particles.ConservativeResizeRows(particles.GetNumRows()+1000);
				}

				//create random coordinate
				particles(numParticles,0) = rRadiusBoundaryParticles+(D-2.*rRadiusBoundaryParticles)*dsfmt_genrand_close_open(&randomNumberGenerator);
				particles(numParticles,1) = rRadiusBoundaryParticles+(1.5*D-2.*rRadiusBoundaryParticles)*dsfmt_genrand_close_open(&randomNumberGenerator);
				particles(numParticles,2) = rBoundingBox(2,theFace);
				particles(numParticles,3) = rRadiusBoundaryParticles;

				bool noSeparation(true);
				//check for overlapping with the boundary
		    	//right circle
				double deltaX = particles(numParticles,0)-(rBoundingBox(0,1)+0.525*D);
				double deltaY = particles(numParticles,1)-(rBoundingBox(1,0)+0.75*D);
				double sumR = particles(numParticles,3)+0.725*D;
				if (deltaX*deltaX+deltaY*deltaY<sumR*sumR)
					noSeparation=false;
				//left circle
				deltaX = particles(numParticles,0)-(rBoundingBox(0,0)-0.525*D);
				deltaY = particles(numParticles,1)-(rBoundingBox(1,0)+0.75*D);
				if (deltaX*deltaX+deltaY*deltaY<sumR*sumR)
					noSeparation=false;;

				//check for overlap
				if (noSeparation==true)
				{
					for (int countParticle=0; countParticle<numParticles; countParticle++)
					{
						double deltaX = particles(countParticle,0)-particles(numParticles,0);
						double deltaY = particles(countParticle,1)-particles(numParticles,1);
						double deltaZ = particles(countParticle,2)-particles(numParticles,2);
						double sumR = particles(countParticle,3)+particles(numParticles,3);
						if (deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ<sumR*sumR)
						{
							noSeparation = false;
							//std::cout << "failure at " << particles(numParticles,0 ) << " " << particles(numParticles,1 ) << "\n";
							break;
						}
					}
				}

				if (noSeparation==true)
				{
					thisParticlesPerFace++;
					numParticles++;
					numTries=0;
					//std::cout << "Particles " << "\n" << particles.mEigenMatrix.block(0,0,numParticles,4) << "\n";
				}
				else
				{
					numTries++;
					if (numTries>10000)
						throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] could not add particles on a boundary face.");
				}
			}
		}
		break;
	case 2:
		//cylinder - add the top and bottom surface
		for (unsigned int theFace=0; theFace<2; theFace++)
		{
			double D = rBoundingBox(0,1) - rBoundingBox(0,0);
			double lArea = 0.25*M_PI*D*D;
			if (fabs(rBoundingBox.at(1,1) - rBoundingBox.at(1,0)-D)>1e-10)
				throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");

			int numParticlesPerFace = lArea/(rDistanceBoundaryParticles*rDistanceBoundaryParticles);

			//the edges have already been inserted
			int thisParticlesPerFace(4*((int)(0.5*M_PI*D/rDistanceBoundaryParticles))-4);

			int numTries(0);
			while(thisParticlesPerFace<numParticlesPerFace)
			{
				// check size of table
				if(numParticles == particles.GetNumRows())
				{
					particles.ConservativeResizeRows(particles.GetNumRows()+1000);
				}

				//create random coordinate
				particles(numParticles,0) = rRadiusBoundaryParticles+(D-2.*rRadiusBoundaryParticles)*dsfmt_genrand_close_open(&randomNumberGenerator);
				particles(numParticles,1) = rRadiusBoundaryParticles+(D-2.*rRadiusBoundaryParticles)*dsfmt_genrand_close_open(&randomNumberGenerator);
				particles(numParticles,2) = rBoundingBox(2,theFace);
				particles(numParticles,3) = rRadiusBoundaryParticles;

				bool noSeparation(true);
				//check for overlapping with the boundary
		    	//right circle
				double deltaX = particles(numParticles,0)-(rBoundingBox(0,0)+0.5*D);
				double deltaY = particles(numParticles,1)-(rBoundingBox(1,0)+0.5*D);
				double sumR = 0.5*D-particles(numParticles,3);
				if (sumR<0)
					throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] that should not have happend.");
				if (deltaX*deltaX+deltaY*deltaY>sumR*sumR)
					noSeparation=false;

				//check for overlap
				if (noSeparation==true)
				{
					for (int countParticle=0; countParticle<numParticles; countParticle++)
					{
						double deltaX = particles(countParticle,0)-particles(numParticles,0);
						double deltaY = particles(countParticle,1)-particles(numParticles,1);
						double deltaZ = particles(countParticle,2)-particles(numParticles,2);
						double sumR = particles(countParticle,3)+particles(numParticles,3);
						if (deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ<sumR*sumR)
						{
							noSeparation = false;
							//std::cout << "failure at " << particles(numParticles,0 ) << " " << particles(numParticles,1 ) << "\n";
							break;
						}
					}
				}

				if (noSeparation==true)
				{
					thisParticlesPerFace++;
					numParticles++;
					numTries=0;
					//std::cout << "Particles " << "\n" << particles.mEigenMatrix.block(0,0,numParticles,4) << "\n";
				}
				else
				{
					numTries++;
					if (numTries>10000)
						throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] could not add particles on a boundary face.");
				}
			}
		}
		break;
	default:
		throw MechanicsException("[NuTo::StructureBase::CreateSpheresOnSpecimenBoundary] unknown specimen type.");
	}



	particles.ConservativeResizeRows(numParticles);
	return particles;
}



//! @brief is only true for structure used as multiscale (structure in a structure)
//! @parameters rTypeOfSpecimen 0 box, 1 dogbone, 2 cylinder
//! @parameters rBoundingBox box for the spheres (3*2 matrix)
//! @parameters rRelParticleVolume percentage of particle volume inside the box
//! @parameters rGradingCurve matrix with each line min_diameter, max_diameter, volume percentage of that sieve size
//! @parameters relativeDistance scaling factor to increase the diameter when inserting the sphere to ensure a minimum distance
//! @parameters absoluteDistance scaling value to increase the diameter when inserting the sphere to ensure a minimum distance
//! @parameters rSeed seed for the random number generator
//! @parameters rSpheresBoundary particles simulated on the boundary e.g. created with CreateSpheresOnBoxBoundary (they do not contribute to the grading curve)
//! @return ... matrix with spheres (coordinates x y z and radius)
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::CreateSpheresInSpecimen(int rTypeOfSpecimen,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, double rRelParticleVolume,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rGradingCurve, double relativeDistance, double absoluteDistance,
		int rSeed, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheresBoundary)
{
    if (rBoundingBox.GetNumRows()!=3 && rBoundingBox.GetNumColumns()!=2)
    	throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] bounding box has to have the dimension [3,2]");

    if (rGradingCurve.GetNumRows()<1)
    	throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] at least one class in the grading curve should be defined.");

	// calculate specimen length
    std::array<double,3> lBox;
    for (int count=0; count<3; count++)
    {
    	lBox[count] = rBoundingBox(count,1) - rBoundingBox(count,0);
    	if (lBox[count]<=0)
    		throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] box dimensions should be not negative.");
    }

    // volume of particles per class
    std::vector<double> Vsoll, Vist;
    Vsoll.resize(rGradingCurve.GetNumRows());
    Vist.resize(rGradingCurve.GetNumRows());
    std::vector<int> numParticlesPerClass;
    numParticlesPerClass.resize(rGradingCurve.GetNumRows());

    // calculating volume of the specimen
    double Vspecimen;
    switch (rTypeOfSpecimen)
    {
    case 0:
    	//box
    	Vspecimen = lBox[0] * lBox[1] * lBox[2];
    break;
    case 1:
    {
    	Vspecimen = lBox[0] * lBox[1] * lBox[2];
		double D = rBoundingBox(0,1) - rBoundingBox(0,0);
		if (fabs(rBoundingBox.at(1,1) - rBoundingBox.at(1,0)-1.5*D)>1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");
		//subtract the circles
		double radius = 0.725*D; ;
		double deltaAngle = 2.*0.2/0.525;
		Vspecimen -= 2.*(deltaAngle/(2.*M_PI)*M_PI*radius*radius)* lBox[2];
    }
    break;
    case 2:
    {
		double D = rBoundingBox(0,1) - rBoundingBox(0,0);
    	Vspecimen =  M_PI*0.25*D*D*lBox[2];
		if (fabs(lBox[0]-lBox[1])>1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the cylindern, the x and y dimension should be identical (Diameter).");
		if (D<1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the cylindern, the x,y dimension should be positive (Diameter).");
    }
    break;
    default:
    	throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] specimen type not implemented.");
    }

    if(Vspecimen < 1.0e-14)
    {
        throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] negative volume of the box.");
    }

	// calculating mass of the aggregates */
	double volumeSumParticles = Vspecimen * rRelParticleVolume;
	int numParticles(rSpheresBoundary.GetNumRows());
	FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> particles(0,4);
	if (numParticles>0)
		particles = rSpheresBoundary;

	// random number generator
	dsfmt_t randomNumberGenerator;
	// init random number generator with milliseconds from ..
	dsfmt_init_gen_rand(&randomNumberGenerator, rSeed);
	double lastPrintedFraction = 0.;
	for(int gc=0;gc<rGradingCurve.GetNumRows();gc++)
	{
		double dMin(rGradingCurve(gc,0));
		double dMax(rGradingCurve(gc,1));
        if (dMin>dMax)
		    throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] the minimum radius is larger than the maximum radius.");
		double volumeFrac(rGradingCurve(gc,2));
		if (volumeFrac<0 || volumeFrac>1)
	        throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] the mass fraction should be in the range [0,1].");

		numParticlesPerClass[gc]=0;

		// calculate reference volume fraction of the mineral-size-class
		Vsoll[gc] = volumeSumParticles * volumeFrac;
		if (gc>0)
		{
			Vsoll[gc] = Vsoll[gc] + (Vsoll[gc-1] - Vist[gc-1]);
		}
		Vist[gc] = 0.0;

		// generate particles until the reference volume fraction is reached
		bool finished(false);
		while (!finished)
		{
			// check size of table
			if(numParticles == particles.GetNumRows())
			{
				particles.ConservativeResizeRows(particles.GetNumRows()+1000);
			}

			// calculate radius and volume of the particle
			double randomNumber(dsfmt_genrand_close_open(&randomNumberGenerator));// = geometry_rng.rand();

			double radius = 0.5 * dMin*dMax/pow((1.0-randomNumber)*(dMax*dMax*dMax) +
					randomNumber*(dMin*dMin*dMin),1.0/3.0);

			// volume
			double volumeParticle (4.0/3.0*M_PI*radius*radius*radius);

			Vist[gc] += volumeParticle;

			//create new particle
			if (Vist[gc] < Vsoll[gc])
			{
			    //std::cout << "sphere " << numParticles+1 << " " << radius << " volume " << volumeParticle << std::endl;
				particles(numParticles,3) = radius;
			    numParticles++;
			    numParticlesPerClass[gc]++;
			}
			else
			{
				finished=true;
				Vist[gc] -= volumeParticle;
			}
		}
		std::cout << "Volume for class " << gc+1 << " : " <<  Vist[gc]/Vspecimen << "(" << Vsoll[gc]/Vspecimen << ")" << std::endl;
		//sort only the newly introduced radii
		std::sort(((double*)&particles.data()[3*particles.GetNumRows()+numParticles-numParticlesPerClass[gc]]),((double*)&particles.data()[3*particles.GetNumRows()+numParticles]), std::greater<double>( ));

		//create boxes for the previously inserted particles
		//width of each box = largest diameter
		std::array<int,3> nSubBox;
		std::array<double,3> lSubBox;
		for (int count=0; count<3; count++)
		{
			nSubBox[count] = lBox[count]/(dMax);
			lSubBox[count] = lBox[count]/nSubBox[count];
		}

		int numberOfSubBoxes = nSubBox[0]*nSubBox[1]*nSubBox[2];
		std::vector<std::vector<int> > subBox(numberOfSubBoxes);

		// now start inserting existing particles into the box
		for (int countParticle=0; countParticle<numParticles-numParticlesPerClass[gc]; countParticle++)
		{
			InsertParticleIntoBox(particles,countParticle,subBox,nSubBox,lSubBox,rBoundingBox);
		}

		// now start inserting new particles into the box
		for (int countParticle=numParticles-numParticlesPerClass[gc]; countParticle<numParticles; countParticle++)
		{
			bool inserted(false);
			int numTries(0);
			while(!inserted)
			{
				//create random coordinate
				std::array<int,3> cSubBox;
				for (int count=0; count<3; count++)
				{
					particles(countParticle,count) = particles(countParticle,3)+(lBox[count]-2*particles(countParticle,3))*dsfmt_genrand_close_open(&randomNumberGenerator);
					cSubBox[count] = (particles(countParticle,count)-rBoundingBox(count,0))/lSubBox[count];
				}

				//check for overlapping with the boundary
			    switch (rTypeOfSpecimen)
			    {
			    case 0:
			    	//box, do nothing
			    break;
			    case 1:
			    {
			    	//dogbone specimen right circle
					double D = rBoundingBox(0,1) - rBoundingBox(0,0);
					double deltaX = particles(countParticle,0)-(rBoundingBox(0,1)+0.525*D);
					double deltaY = particles(countParticle,1)-(rBoundingBox(1,0)+0.75*D);
					double sumR = particles(countParticle,3)+0.725*D;
					if (deltaX*deltaX+deltaY*deltaY<sumR*sumR)
						continue;
					//left circle
					deltaX = particles(countParticle,0)-(rBoundingBox(0,0)-0.525*D);
					deltaY = particles(countParticle,1)-(rBoundingBox(1,0)+0.75*D);
					if (deltaX*deltaX+deltaY*deltaY<sumR*sumR)
						continue;
			    }
				break;
			    case 2:
			    {
					//check for overlapping with the boundary
					//circle
			    	double D = rBoundingBox(0,1) - rBoundingBox(0,0);
					double deltaX = particles(countParticle,0)-(rBoundingBox(0,0)+0.5*D);
					double deltaY = particles(countParticle,1)-(rBoundingBox(1,0)+0.5*D);
					double sumR = 0.5*D-particles(countParticle,3)*(1.+relativeDistance)-absoluteDistance;
					if (sumR<0)
						throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] that should not have happend.");
					if (deltaX*deltaX+deltaY*deltaY>sumR*sumR)
						continue;
			    }
			    break;
			    default:
			    	throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] specimen type not implemented.");
			    }

				//calculate the corresponding box
				int theBox = cSubBox[0]*nSubBox[1]*nSubBox[2]+cSubBox[1]*nSubBox[2]+cSubBox[2];

				//check for overlap with all the ellipses in that box
				bool noSeparation(true);
				for (unsigned int countOtherEllipse=0; countOtherEllipse<subBox[theBox].size(); countOtherEllipse++)
				{
					int otherEllipse=subBox[theBox][countOtherEllipse];
					double deltaX = particles(countParticle,0)-particles(otherEllipse,0);
					double deltaY = particles(countParticle,1)-particles(otherEllipse,1);
					double deltaZ = particles(countParticle,2)-particles(otherEllipse,2);
					double sumR = particles(countParticle,3)*(1.+relativeDistance)+absoluteDistance+particles(otherEllipse,3);
					if (deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ<sumR*sumR)
					{
						noSeparation = false;
						break;
					}
				}

				if (noSeparation)
				{
					//insert
					inserted = true;
					InsertParticleIntoBox(particles,countParticle,subBox,nSubBox,lSubBox,rBoundingBox);
					if ((double)(countParticle)/numParticles-lastPrintedFraction>0.1)
					{
						std::cout << (double)(countParticle)/numParticles*100. << "% of particles already inserted." << "\n";
						lastPrintedFraction = (double)(countParticle)/numParticles;
					}
				}
				else
				{
					numTries++;
					if (numTries>1e7)
						throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] unable to insert sphere after a 1e7 tries.");
				}
			}
		}
	}

	return particles.GetBlock(rSpheresBoundary.GetNumRows(),0,numParticles-rSpheresBoundary.GetNumRows(),4);
}

//! @brief cut spheres at a given z-coordinate to create circles (in 2D)
//! @parameters rSpheres matrix with the spheres (x,y,z,r)
//! @parameters rZCoord z coordinate (where to cut)
//! @parameters rMinRadius minimal radius of the circle
//! @return ... matrix with the circles (x,y,r)
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::CutSpheresZ(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheres, double rZCoord, double rMinRadius)
{
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> circles(1000,3);
	int numCircles(0);
	for (int countSphere=0; countSphere<rSpheres.GetNumRows(); countSphere++)
	{
        double delta=rSpheres(countSphere,2)-rZCoord;
		if (fabs(delta)<rSpheres(countSphere,3))
        {
			double radius=sqrt(rSpheres.at(countSphere,3)*rSpheres.at(countSphere,3)-delta*delta);
        	if (radius>rMinRadius)
        	{
        		//add circle
        		if (numCircles==circles.GetNumRows())
        		{
        			circles.ConservativeResizeRows(numCircles+1000);
        		}
        		circles(numCircles,0) = rSpheres(countSphere,0);
        		circles(numCircles,1) = rSpheres(countSphere,1);
        		circles(numCircles,2) = radius;
        		numCircles++;
        	}
        }
	}
	circles.ConservativeResizeRows(numCircles);

    return circles;
}


//! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
void NuTo::StructureBase::InsertParticleIntoBox(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rParticles, int rTheParticle, std::vector<std::vector<int > >& rSubBox, std::array<int,3>& rNSubBox,std::array<double,3>& rLSubBox, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox)
{

	//calculate current coordinate box of the center
	std::array<int,3> cSubBoxMin;
	std::array<int,3> cSubBoxMax;
	for (int count=0; count<3; count++)
	{
		cSubBoxMax[count] = (rParticles(rTheParticle,count)+rParticles(rTheParticle,3)-rBoundingBox(count,0))/rLSubBox[count]+1;
		if (cSubBoxMax[count]>=rNSubBox[count])
			cSubBoxMax[count]=rNSubBox[count]-1;

		cSubBoxMin[count] = (rParticles(rTheParticle,count)-rParticles(rTheParticle,3)-rBoundingBox(count,0))/rLSubBox[count]-1;
		if (cSubBoxMin[count]<0)
			cSubBoxMin[count]=0;
	}

	//insert the center box + all the surroundings
	for (int countx=cSubBoxMin[0]; countx<=cSubBoxMax[0]; countx++)
	{
		for (int county=cSubBoxMin[1]; county<=cSubBoxMax[1]; county++)
		{
			for (int countz=cSubBoxMin[2]; countz<=cSubBoxMax[2]; countz++)
			{
				int theBox = countx*rNSubBox[1]*rNSubBox[2]+county*rNSubBox[2]+countz;
				rSubBox[theBox].push_back(rTheParticle);
			}
		}
	}
}

//@brief determines if in the omp parallelization the maximum independent sets are used (parallel assembly of the stiffness, generally faster)
// or sequential insertion of the element stiffness using a barrier (faster for different load balancing of the elements)
// is only relevant for openmp, otherwise the routine is just empty
void NuTo::StructureBase::UseMaximumIndependentSets(bool rUseMIS)
{
#ifdef _OPENMP
	mUseMIS=rUseMIS;
#endif //_OPENMP
}

//@brief set the number of processors for openmp parallelization
void NuTo::StructureBase::SetNumProcessors(int rNumProcessors)
{
#ifdef _OPENMP
	mNumProcessors = rNumProcessors;
#endif //_OPENMP
}
//@brief get the number of processors for openmp parallelization
int NuTo::StructureBase::GetNumProcessors() const
{
#ifdef _OPENMP
    return mNumProcessors;
#endif //_OPENMP
    return 1;
}

void NuTo::StructureBase::SetOMPNested(bool rNested)
{
#ifdef _OPENMP
	omp_set_nested(rNested);
#endif //_OPENMP
}

//! @brief sets the Hessian to be constant or variable
//! @parameters rTimeDerivative (0 = stiffness, 1 damping, 2 mass)
//! @parameters rValue (true = const false=variable)
void NuTo::StructureBase::SetHessianConstant(int rTimeDerivative, bool rValue)
{
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<3);
	mHessianConstant[rTimeDerivative] = rValue;
}

//! @brief sets the Hessian to be constant or variable
//! @parameters rTimeDerivative (0 = stiffness, 1 damping, 2 mass)
//! @return (true = const false=variable)
bool NuTo::StructureBase::GetHessianConstant(int rTimeDerivative)const
{
	assert(rTimeDerivative>=0);
	assert(rTimeDerivative<3);
	return mHessianConstant[rTimeDerivative];
}

bool NuTo::StructureBase::InterpolationTypeIsConstitutiveInput(NuTo::Node::eAttributes rDofType)
{
	InterpolationType* interpolationType;
    for (auto interpolation = mInterpolationTypeMap.begin(); interpolation != mInterpolationTypeMap.end(); interpolation++){
    	interpolationType = interpolation->second;
    	if(interpolationType->IsConstitutiveInput(rDofType)){
    		return true;
    	}
    }

    return false;
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
