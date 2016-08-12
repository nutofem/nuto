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


# ifdef _OPENMP
#include <omp.h>
# endif

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/base/Timer.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
#include "nuto/math/SparseDirectSolverPardiso.h"

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
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include "nuto/visualize/VisualizeComponent.h"
#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#endif // ENABLE_VISUALIZE


NuTo::StructureBase::StructureBase(int rDimension)  : NuTo::NuToObject::NuToObject(),
        mConstraintMatrix(mDofStatus, false),
        mConstraintMappingRHS(mDofStatus, false),
        mConstraintRHS(mDofStatus)
{
    if (rDimension!=1 && rDimension!=2 && rDimension!=3)
    {
        throw MechanicsException("[StructureBase::StructureBase] The dimension of a structure is either 1, 2 or 3.");
    }
    mDimension = rDimension;
    mPrevTime = 0.;
    mTime = 0.;
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

#ifdef _OPENMP
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
    ar & BOOST_SERIALIZATION_NVP(mPrevTime);
    ar & BOOST_SERIALIZATION_NVP(mTime);
    ar & BOOST_SERIALIZATION_NVP(mDimension);
    ar & BOOST_SERIALIZATION_NVP(mConstitutiveLawMap);
    ar & BOOST_SERIALIZATION_NVP(mConstraintMap);
    ar & BOOST_SERIALIZATION_NVP(mNumLoadCases);
    ar & BOOST_SERIALIZATION_NVP(mLoadMap);
    ar & BOOST_SERIALIZATION_NVP(mGroupMap);
    ar & BOOST_SERIALIZATION_NVP(mIntegrationTypeMap);
    ar & BOOST_SERIALIZATION_NVP(mSectionMap);
    ar & BOOST_SERIALIZATION_NVP(mInterpolationTypeMap);
    ar & BOOST_SERIALIZATION_NVP(mMappingIntEnum2String);
#ifdef ENABLE_VISUALIZE
//    ar & BOOST_SERIALIZATION_NVP(mGroupVisualizeComponentsMap);
    std::cout << "WARNING: Visualization components are not serialized!\n";
#endif
    ar & BOOST_SERIALIZATION_NVP(mDofStatus);
    ar & BOOST_SERIALIZATION_NVP(mNodeNumberingRequired);
    ar & BOOST_SERIALIZATION_NVP(mConstraintMatrix);
    ar & BOOST_SERIALIZATION_NVP(mConstraintMappingRHS);
    ar & BOOST_SERIALIZATION_NVP(mConstraintRHS);
    ar & BOOST_SERIALIZATION_NVP(mHaveTmpStaticData);
    ar & BOOST_SERIALIZATION_NVP(mUpdateTmpStaticDataRequired);
    ar & BOOST_SERIALIZATION_NVP(mToleranceStiffnessEntries);
#ifdef _OPENMP
    // mMIS contains Element-Ptr, so its serialized, by serializing the Pointer-Addresses and updating them afterwards
    // see Structure::loadImplement and Structure::saveImplement
    ar & BOOST_SERIALIZATION_NVP(mNumProcessors);
#endif // _OPENMP
    ar & boost::serialization::make_nvp("mLogger", mLogger);
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

    mLogger << "num dofs : " << GetNumTotalDofs() << "\n";
    mLogger << "num active dofs : " << GetNumTotalActiveDofs() << "\n";
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

//! @brief set number of cycles to be extrapolated in the cycle jump routine
//! @brief ... rNumber[0] is the number of extrapolated cycles itself Njump
//! @brief ... rNumber[1] is the weighting coefficient of the implicit term
//! @brief ... rNumber[2] is the weighting coefficient of the explicit term
//! @brief ... rNumber[3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
void NuTo::StructureBase::SetNumExtrapolatedCycles(NuTo::FullVector<double,Eigen::Dynamic> rNumber)
{
	if (rNumber.size()<3)
	        throw NuTo::MechanicsException("[NuTo::StructureBase::SetNumExtrapolatedCycles] at least number of extrapolation cycles and weighting coefficient for explicit and implicit terms are required.");

	// check the number of extrapolation cycles
	if (rNumber[0] < 0)
			throw NuTo::MechanicsException("[NuTo::StructureBase::SetNumExtrapolatedCycles] number of extrapolation cycles is negative.");

	// check the weighting coefficients
	for (int i = 1; i < rNumber.size(); ++i) {
		if (rNumber[i] < 0 || rNumber[i]>1)
				throw NuTo::MechanicsException("[NuTo::StructureBase::SetNumExtrapolatedCycles] the " + std::to_string(i) + " weighting coefficient is out of the [0,1] range.");
	}
	mNumExtrapolatedCycles = rNumber;
}

//! @brief get the number of cycles to be extrapolated. Returns
//! @brief ... [0] is the number of extrapolated cycles itself Njump
//! @brief ... [1] is the weighting coefficient of the implicit term
//! @brief ... [2] is the weighting coefficient of the explicit term
//! @brief ... [3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
NuTo::FullVector<double,Eigen::Dynamic> NuTo::StructureBase::GetNumExtrapolatedCycles() const
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
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Element group does not exist.");

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

    if (mVerboseLevel > 5)
    {
        for (auto const &iPair : mGroupVisualizeComponentsMap)
        {
            std::cout << "ele group: \t" << iPair.first << std::endl;
            for (auto const &iComponentPtr : iPair.second)
            {
                std::cout << "components: \t " << iComponentPtr->GetComponentName() << std::endl;
            }
        }
    }
#endif // ENABLE_VISUALIZE
}


// add visualization components for an element group
void NuTo::StructureBase::AddVisualizationComponent(int rElementGroup, const std::string& rVisualizeComponent)
{
    std::cout << rVisualizeComponent << std::endl;
    if (rVisualizeComponent == "Accelerations")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ACCELERATION);
    else if (rVisualizeComponent == "AngularAccelerations")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ANGULAR_ACCELERATION);
    else if (rVisualizeComponent == "AngularVelocities")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ANGULAR_VELOCITY);
    else if (rVisualizeComponent == "BondStress")
        AddVisualizationComponent(rElementGroup, VisualizeBase::BOND_STRESS);
    else if (rVisualizeComponent == "Constitutive")
        AddVisualizationComponent(rElementGroup, VisualizeBase::CONSTITUTIVE);
    else if (rVisualizeComponent == "Crack")
        AddVisualizationComponent(rElementGroup, VisualizeBase::CRACK);
    else if (rVisualizeComponent == "Damage")
        AddVisualizationComponent(rElementGroup, VisualizeBase::DAMAGE);
    else if (rVisualizeComponent == "Displacements")
        AddVisualizationComponent(rElementGroup, VisualizeBase::DISPLACEMENTS);
    else if (rVisualizeComponent == "Element")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ELEMENT);
    else if (rVisualizeComponent == "EngineeringPlasticStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_PLASTIC_STRAIN);
    else if (rVisualizeComponent == "EngineeringStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_STRAIN);
    else if (rVisualizeComponent == "EngineeringStress")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ENGINEERING_STRESS);
    else if (rVisualizeComponent == "HeatFlux")
        AddVisualizationComponent(rElementGroup, VisualizeBase::HEAT_FLUX);
    else if (rVisualizeComponent == "LatticeStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::LATTICE_STRAIN);
    else if (rVisualizeComponent == "LatticeStress")
        AddVisualizationComponent(rElementGroup, VisualizeBase::LATTICE_STRESS);
    else if (rVisualizeComponent == "LocalEqStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::LOCAL_EQ_STRAIN);
    else if (rVisualizeComponent == "NonlocalEqStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::NONLOCAL_EQ_STRAIN);
    else if (rVisualizeComponent == "ParticleRadius")
        AddVisualizationComponent(rElementGroup, VisualizeBase::PARTICLE_RADIUS);
    else if (rVisualizeComponent == "PrincipalEngineeringStress")
        AddVisualizationComponent(rElementGroup, VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);
    else if (rVisualizeComponent == "RelativeHumidity")
        AddVisualizationComponent(rElementGroup, VisualizeBase::RELATIVE_HUMIDITY);
    else if (rVisualizeComponent == "Rotations")
        AddVisualizationComponent(rElementGroup, VisualizeBase::ROTATION);
    else if (rVisualizeComponent == "Section")
        AddVisualizationComponent(rElementGroup, VisualizeBase::SECTION);
    else if (rVisualizeComponent == "ShrinkageStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::SHRINKAGE_STRAIN);
    else if (rVisualizeComponent == "Slip")
        AddVisualizationComponent(rElementGroup, VisualizeBase::SLIP);
    else if (rVisualizeComponent == "Temperature")
        AddVisualizationComponent(rElementGroup, VisualizeBase::TEMPERATURE);
    else if (rVisualizeComponent == "TotalInelasticEqStrain")
        AddVisualizationComponent(rElementGroup, VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN);
    else if (rVisualizeComponent == "Velocities")
        AddVisualizationComponent(rElementGroup, VisualizeBase::VELOCITY);
    else if (rVisualizeComponent == "WaterVolumeFraction")
        AddVisualizationComponent(rElementGroup, VisualizeBase::WATER_VOLUME_FRACTION);
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Visualization component not implemented or misspelled.");



}
void NuTo::StructureBase::AddVisualizationComponentNonlocalWeights(int rElementGroup, int rElementId, int rIp)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Element group does not exist.");

    const ElementBase *elementBase = ElementGetElementPtr(rElementId);
    int numIp = elementBase->GetNumIntegrationPoints();

    if (rIp < 0 or rIp >= numIp)
        throw MechanicsException(__PRETTY_FUNCTION__, "Integration point number is out of range.");

    try
    {
        // create a new visualization list for an element group or add components to an already existing list
        if (mGroupVisualizeComponentsMap.find(rElementGroup) == mGroupVisualizeComponentsMap.end())
        {
            std::list<std::shared_ptr<VisualizeComponent>> visualizationPtrList;
            visualizationPtrList.push_back(std::make_shared<VisualizeComponentNonlocalWeight>(VisualizeComponentNonlocalWeight(rElementId, rIp)));

            mGroupVisualizeComponentsMap.insert(std::pair<int,std::list<std::shared_ptr<VisualizeComponent>>>(rElementGroup, visualizationPtrList));
            // mGroupVisualizeComponentsMap.emplace(rElementGroup, visualizationPtrList);       //<- use this for gcc version 4.9 or higher!

            mGroupVisualizationType.insert(std::pair<int, VisualizeBase::eVisualizationType>(rElementGroup, VisualizeBase::VORONOI_CELL));
            // mGroupVisualizationType.emplace(rElementGroup, VisualizeBase::VORONOI_CELL);     //<- use this for gcc version 4.9 or higher!


        } else
        {
            mGroupVisualizeComponentsMap.at(rElementGroup).push_back(std::make_shared<VisualizeComponentNonlocalWeight>(VisualizeComponentNonlocalWeight(rElementId, rIp)));
        }

    }
    catch (NuTo::MechanicsException &e)
     {
        e.AddMessage(__PRETTY_FUNCTION__, "error setting element and local ip number.");
        throw e;
     }
    catch(...)
     {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "error setting element and local ip number.");
     }
#endif // ENABLE_VISUALIZE
}


void NuTo::StructureBase::SetVisualizationType(const int rElementGroup, const VisualizeBase::eVisualizationType rVisualizationType)
{
    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Element group does not exist.");

    // check if the element group exists
    if (mGroupVisualizationType.find(rElementGroup) == mGroupVisualizationType.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Please add a visualization component first before setting the visualization type.");

    mGroupVisualizationType.at(rElementGroup) = rVisualizationType;
}


void NuTo::StructureBase::ClearVisualizationComponents()
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    mGroupVisualizeComponentsMap.clear();

#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ExportVtkDataFileNodes(const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

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
#endif // ENABLE_VISUALIZE
}



void NuTo::StructureBase::ExportVtkDataFileElements(const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

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

#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rResultFileName, bool rXML)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeElementData(visualize,mGroupVisualizeComponentsMap.at(rGroupIdent));
    this->ElementGroupAddToVisualize(rGroupIdent,visualize,mGroupVisualizeComponentsMap.at(rGroupIdent));

    if (rXML)
        visualize.ExportVtuDataFile(rResultFileName);
    else
        visualize.ExportVtkDataFile(rResultFileName);

#endif // ENABLE_VISUALIZE
}

std::map<int, std::list<std::shared_ptr<NuTo::VisualizeComponent>>>& NuTo::StructureBase::GetGroupVisualizeComponentsMap(void)
{
    return mGroupVisualizeComponentsMap;
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

        case NuTo::VisualizeBase::HEAT_FLUX:
        case NuTo::VisualizeBase::SLIP:
        case NuTo::VisualizeBase::CRACK:
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
        case NuTo::VisualizeBase::LATTICE_STRESS:
        case NuTo::VisualizeBase::LATTICE_STRAIN:
        case NuTo::VisualizeBase::LATTICE_PLASTIC_STRAIN:
            rVisualize.DefineCellDataVector(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::BOND_STRESS:
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        case NuTo::VisualizeBase::SHRINKAGE_STRAIN:
        case NuTo::VisualizeBase::THERMAL_STRAIN:
            rVisualize.DefineCellDataTensor(it.get()->GetComponentName());
            break;

        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
        case NuTo::VisualizeBase::TEMPERATURE:
        case NuTo::VisualizeBase::CRACK_PHASE_FIELD:
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
        	throw MechanicsException(__PRETTY_FUNCTION__, "undefined visualize components.");
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





//! @brief ... evaluates the structure
void NuTo::StructureBase::Evaluate(const ConstitutiveInputMap& rInput, std::map<StructureEnum::eOutput, StructureOutputBase *> &rStructureOutput)
{
    throw MechanicsException("[NuTo::StructureBase::Evaluate] Not implemented.");
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian(StructureEnum::eOutput rOutput)
{
    Timer timer(std::string(__FUNCTION__) + ": " + StructureEnum::OutputToString(rOutput), GetShowTime(), GetLogger());
    if (mNodeNumberingRequired) NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    std::set<StructureEnum::eOutput> supportedTypes({StructureEnum::HESSIAN0, StructureEnum::HESSIAN1, StructureEnum::HESSIAN2, StructureEnum::HESSIAN2_LUMPED});
    if (supportedTypes.find(rOutput) == supportedTypes.end())
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] " + StructureEnum::OutputToString(rOutput) + " is not a matrix type or is not supported right now.");

    StructureOutputBlockMatrix hessian(mDofStatus, true);

    std::map<StructureEnum::eOutput, StructureOutputBase *> evaluateMap;
    evaluateMap[rOutput] = &hessian;

    ConstitutiveInputMap input;
    input[Constitutive::Input::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
            CalculateStaticData::EULER_BACKWARD);

    Evaluate(input, evaluateMap);

    return hessian;
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::BuildGlobalInternalGradient()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    if (mNodeNumberingRequired) NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    StructureOutputBlockVector internalGradient(mDofStatus, true);

    std::map<StructureEnum::eOutput, StructureOutputBase *> evaluateMap;
    evaluateMap[StructureEnum::INTERNAL_GRADIENT] = &internalGradient;

    ConstitutiveInputMap input;
    input[Constitutive::Input::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
            CalculateStaticData::EULER_BACKWARD);

    Evaluate(input, evaluateMap);

    return internalGradient;
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian0_CDF(double rDelta)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    bool showTime = GetShowTime();
    SetShowTime(false);


    /*
     * TT:
     * Warning: this code is ugly, repetitive and error-prone.
     * This is mainly due to the J/K  and JJ, JK, KJ, KK stuff.
     * Please be very careful and test properly.
     */

    StructureOutputBlockMatrix hessian0_CDF(GetDofStatus(), true);
    auto internalGradient0 = BuildGlobalInternalGradient();
    auto dofValues = NodeExtractDofValues(0);


    for (auto dofCol : GetDofStatus().GetActiveDofTypes())
    {
        // modify active dof values --> entries JJ (R.J) and KJ (R.K)
        auto& columnDofValues = dofValues.J[dofCol];

        for (int iCol = 0; iCol < columnDofValues.rows(); ++iCol)
        {
            columnDofValues[iCol] += rDelta;
            NodeMergeDofValues(0, dofValues);
            auto internalGradient1 = BuildGlobalInternalGradient();

            StructureOutputBlockVector column = (internalGradient1 - internalGradient0) * (1/rDelta);
            for (auto dofRow : GetDofStatus().GetActiveDofTypes())
            {
                // set JJ entries
                auto& colJ = column.J[dofRow];
                auto& matrixJJ = hessian0_CDF.JJ(dofRow, dofCol);
                for (int i = 0; i < colJ.rows(); ++i)
                    matrixJJ.AddValue(i, iCol, colJ.at(i,0));

                // set KJ entries
                auto& colK = column.K[dofRow];
                auto& matrixKJ = hessian0_CDF.KJ(dofRow, dofCol);
                for (int i = 0; i < colK.rows(); ++i)
                    matrixKJ.AddValue(i, iCol, colK.at(i,0));
            }
            columnDofValues[iCol] -= rDelta;
        }


        // modify dependent dof values --> entries JK (R.J) and KK (R.K)
        auto& rowDofValues = dofValues.K[dofCol];

        for (int iCol = 0; iCol < rowDofValues.rows(); ++iCol)
        {
            rowDofValues[iCol] += rDelta;
            NodeMergeDofValues(0, dofValues);
            auto internalGradient1 = BuildGlobalInternalGradient();

            StructureOutputBlockVector column = (internalGradient1 - internalGradient0) * (1/rDelta);
            for (auto dofRow : GetDofStatus().GetActiveDofTypes())
            {
                // set JK entries
                auto& colJ = column.J[dofRow];
                auto& matrixJK = hessian0_CDF.JK(dofRow, dofCol);
                for (int i = 0; i < colJ.rows(); ++i)
                    matrixJK.AddValue(i, iCol, colJ.at(i,0));

                // set KJ entries
                auto& colK = column.K[dofRow];
                auto& matrixKK = hessian0_CDF.KK(dofRow, dofCol);
                for (int i = 0; i < colK.rows(); ++i)
                    matrixKK.AddValue(i, iCol, colK.at(i,0));
            }
            rowDofValues[iCol] -= rDelta;
        }
    }
    NodeMergeDofValues(0, dofValues);




    SetShowTime(showTime);

    return hessian0_CDF;
}

bool NuTo::StructureBase::CheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    bool isHessianCorrect = true;

    NodeBuildGlobalDofs(__FUNCTION__);

    bool hasInteractingConstraints = mDofStatus.HasInteractingConstraints();
    mDofStatus.SetHasInteractingConstraints(true); // this ensures the full assembly of KJ and KK, which could be skipped if CMat.Entries = 0


    auto hessian0 = BuildGlobalHessian0();
    auto hessian0_CDF = BuildGlobalHessian0_CDF(rDelta);

    isHessianCorrect = isHessianCorrect && CheckHessian0_Submatrix(hessian0.JJ, hessian0_CDF.JJ, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect && CheckHessian0_Submatrix(hessian0.JK, hessian0_CDF.JK, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect && CheckHessian0_Submatrix(hessian0.KJ, hessian0_CDF.KJ, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect && CheckHessian0_Submatrix(hessian0.KK, hessian0_CDF.KK, rRelativeTolerance, rPrintWrongMatrices);

    mDofStatus.SetHasInteractingConstraints(hasInteractingConstraints);

    return isHessianCorrect;
}

bool NuTo::StructureBase::CheckHessian0_Submatrix(const BlockSparseMatrix& rHessian0, BlockSparseMatrix& rHessian0_CDF, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    int row, col;
    assert(rHessian0.GetNumRows() == rHessian0_CDF.GetNumRows());
    assert(rHessian0.GetNumColumns() == rHessian0_CDF.GetNumColumns());

    if (rHessian0.GetNumRows() == 0 or rHessian0.GetNumColumns() == 0)
        return true;

    bool isSubmatrixCorrect = true;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");
    for (auto dofRow : GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : GetDofStatus().GetActiveDofTypes())
        {
            FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian0_CDF_Full(rHessian0_CDF(dofRow, dofCol));

            double scaling = 1./rHessian0_CDF(dofRow, dofCol).AbsMax();

            auto& diff = rHessian0_CDF(dofRow, dofCol);
            diff.AddScal(rHessian0(dofRow, dofCol), -1.);

            diff *= scaling;
            double error = diff.AbsMax(row, col);
            if (error > rRelativeTolerance)
            {
                GetLogger() << "[" << __FUNCTION__ << "] max error in (" << Node::DofToString(dofRow)<< "," << Node::DofToString(dofCol) << ") "
                        << error << " at entry (" << row << "," << col << ")\n";
                GetLogger() << "hessian0(" << row << "," << col <<") = " << FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(rHessian0(dofRow, dofCol))(row, col) << "\n";
                GetLogger() << "hessian0_CDF(" << row << "," << col <<") = " << hessian0_CDF_Full(row, col) << "\n";
                isSubmatrixCorrect = false;
                if (rPrintWrongMatrices)
                {
                    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> diffPrint(diff);
                    diffPrint.SetSmallEntriesZero(1.e-10);
                    GetLogger() << "####### relative difference\n" << diffPrint.format(fmt) << "\n";
                    GetLogger() << "####### hessian0\n" << FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(rHessian0(dofRow, dofCol)).format(fmt) << "\n";
                    GetLogger() << "####### hessian0_CDF\n" << hessian0_CDF_Full.format(fmt) << "\n";
                }
            }
        }
    }
    return isSubmatrixCorrect;
}

void NuTo::StructureBase::SolveGlobalSystemStaticElastic(int rLoadCase)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    if (GetNumTimeDerivatives() > 0)
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Only use this method for a system with 0 time derivatives.");

    if (mNodeNumberingRequired) NodeBuildGlobalDofs(__PRETTY_FUNCTION__);


    StructureOutputBlockVector deltaDof_dt0(GetDofStatus(), true);
    deltaDof_dt0.J.SetZero();
    deltaDof_dt0.K = ConstraintGetRHSAfterGaussElimination();

    auto hessian0 = BuildGlobalHessian0();
    auto residual = hessian0 * deltaDof_dt0 - BuildGlobalExternalLoadVector(rLoadCase) + BuildGlobalInternalGradient();

    hessian0.ApplyCMatrix(GetConstraintMatrix());
    residual.ApplyCMatrix(GetConstraintMatrix());

    // reuse deltaDof_dt0
    deltaDof_dt0.J = SolveBlockSystem(hessian0.JJ, residual.J);

    deltaDof_dt0.K = NodeCalculateDependentDofValues(deltaDof_dt0.J);
    NodeMergeDofValues(0, deltaDof_dt0);
}


NuTo::BlockFullVector<double> NuTo::StructureBase::SolveBlockSystem(const BlockSparseMatrix& rMatrix, const BlockFullVector<double>& rVector) const
{
    NuTo::FullVector<double, Eigen::Dynamic> resultForSolver;
    NuTo::SparseMatrixCSRGeneral<double>     matrixForSolver = rMatrix.ExportToCSRGeneral();
    matrixForSolver.SetOneBasedIndexing();

    //allocate solver
#if defined(HAVE_PARDISO) && defined(_OPENMP)
    NuTo::SparseDirectSolverPardiso mySolver(GetNumProcessors(), GetVerboseLevel()); // note: not the MKL version
#else
    NuTo::SparseDirectSolverMUMPS mySolver;
#endif

#ifdef SHOW_TIME
    mySolver.SetShowTime(GetShowTime());
#endif

    mySolver.Solve(matrixForSolver, rVector.Export(), resultForSolver);

    return BlockFullVector<double>(-resultForSolver, GetDofStatus());
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::BuildGlobalExternalLoadVector(int rLoadCase)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    if (mNodeNumberingRequired) NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    StructureOutputBlockVector externalLoad(GetDofStatus(), true);

    if (DofTypeIsActive(Node::DISPLACEMENTS))
    {
        auto& vectorJ = externalLoad.J[Node::DISPLACEMENTS];
        auto& vectorK = externalLoad.K[Node::DISPLACEMENTS];


        // loop over all loads
        boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
        while (loadIter != this->mLoadMap.end())
        {
            loadIter->second->AddLoadToGlobalSubVectors(rLoadCase, vectorJ, vectorK);
            loadIter++;
        }
    }
    if (DofTypeIsActive(Node::TEMPERATURE))
    {
        auto& vectorJ = externalLoad.J[Node::TEMPERATURE];
        auto& vectorK = externalLoad.K[Node::TEMPERATURE];

        // loop over all loads
        boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
        while (loadIter != this->mLoadMap.end())
        {
            loadIter->second->AddLoadToGlobalSubVectors(rLoadCase, vectorJ, vectorK);
            loadIter++;
        }
    }
    return externalLoad;
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

std::set<NuTo::Node::eDof> NuTo::StructureBase::DofTypesGet() const
{
    return GetDofStatus().GetDofTypes();
}

void NuTo::StructureBase::DofTypeDeactivateAll()
{
    std::set<NuTo::Node::eDof> activeDofs = DofTypesGetActive();
    for (auto dof : activeDofs)
        DofTypeSetIsActive(dof, false);

    assert(DofTypesGetActive().empty());
}

void NuTo::StructureBase::DofTypeActivateAll()
{
    DofTypeSetIsActive(DofTypesGet());
}

void NuTo::StructureBase::DofTypeSetIsActive(Node::eDof rDofType, bool rIsActive)
{
    BOOST_FOREACH(auto interpolationTypePair, mInterpolationTypeMap)
    {
        auto& interpolationType = interpolationTypePair.second;
        if (interpolationType->IsDof(rDofType))
            interpolationType->SetIsActive(rIsActive, rDofType);
    }
    UpdateDofStatus();
}

void NuTo::StructureBase::DofTypeSetIsActive(const std::set<Node::eDof>& rActiveDofTypes)
{
    DofTypeDeactivateAll();
    for (auto activeDofType : rActiveDofTypes)
        DofTypeSetIsActive(activeDofType, true);
}

bool NuTo::StructureBase::DofTypeIsActive(Node::eDof rDofType) const
{
	const auto& activeDofTypes = DofTypesGetActive();
	return activeDofTypes.find(rDofType) != activeDofTypes.end();
}

void NuTo::StructureBase::DofTypeSetIsConstitutiveInput(Node::eDof rDofType, bool rIsConstitutiveInput)
{
    BOOST_FOREACH(auto interpolationTypePair, mInterpolationTypeMap)
    {
        auto& interpolationType = interpolationTypePair.second;
        if (interpolationType->IsDof(rDofType))
            interpolationType->SetIsConstitutiveInput(rIsConstitutiveInput, rDofType);
    }
}

void NuTo::StructureBase::DofTypeSetIsSymmetric(Node::eDof rDofType, bool rIsSymmetric)
{
    mDofStatus.SetIsSymmetric(rDofType, rIsSymmetric);
}

bool NuTo::StructureBase::DofTypeIsSymmetric(Node::eDof rDofType) const
{
    return mDofStatus.IsSymmetric(rDofType);
}

const NuTo::DofStatus& NuTo::StructureBase::GetDofStatus() const
{
    return mDofStatus;
}

void NuTo::StructureBase::UpdateDofStatus()
{
    std::set<Node::eDof> dofTypes;
    std::set<Node::eDof> activeDofTypes;
    BOOST_FOREACH(auto interpolationTypePair, mInterpolationTypeMap)
    {
        const std::set<Node::eDof>& dofs = interpolationTypePair.second->GetDofs();
        dofTypes.insert(dofs.begin(), dofs.end());

        const std::set<Node::eDof>& activeDofs = interpolationTypePair.second->GetActiveDofs();
        activeDofTypes.insert(activeDofs.begin(), activeDofs.end());

    }
    dofTypes.erase(Node::COORDINATES);
    activeDofTypes.erase(Node::COORDINATES);

    mDofStatus.SetDofTypes(dofTypes);
    mDofStatus.SetActiveDofTypes(activeDofTypes);

//    std::cout << "CMat entries : " << mConstraintMatrixDof.GetNumActiveEntires() << std::endl;
    //std::cout << "CMat entries : " << mConstraintMatrixDof.GetNumActiveEntires() << std::endl;
    mDofStatus.SetHasInteractingConstraints(mConstraintMatrix.GetNumActiveEntires() != 0);
//    mDofStatus.SetHasInteractingConstraints(true);
}

int NuTo::StructureBase::GetNumTotalDofs() const
{
    return GetNumTotalActiveDofs() + GetNumTotalDependentDofs();
}

int NuTo::StructureBase::GetNumTotalActiveDofs() const
{
    int numTotalActiveDofs = 0;
    for (auto pair : mDofStatus.GetNumActiveDofsMap())
        numTotalActiveDofs += pair.second;
    return numTotalActiveDofs;
}

int NuTo::StructureBase::GetNumTotalDependentDofs() const
{
    int numTotalActiveDofs = 0;
    for (auto pair : mDofStatus.GetNumDependentDofsMap())
        numTotalActiveDofs += pair.second;
    return numTotalActiveDofs;
}

std::set<NuTo::Node::eDof> NuTo::StructureBase::DofTypesGetActive() const
{
    return GetDofStatus().GetActiveDofTypes();
}


int NuTo::StructureBase::GetNumDofs(Node::eDof rDofType) const
{
    return GetNumActiveDofs(rDofType) + GetNumDependentDofs(rDofType);
}

int NuTo::StructureBase::GetNumActiveDofs(Node::eDof rDofType) const
{
    auto it = mDofStatus.GetNumActiveDofsMap().find(rDofType);
    if (it == mDofStatus.GetNumActiveDofsMap().end())
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] There are no " + Node::DofToString(rDofType) + " dofs.");

    return it->second;
}

int NuTo::StructureBase::GetNumDependentDofs(Node::eDof rDofType) const
{
    auto it = mDofStatus.GetNumDependentDofsMap().find(rDofType);
    if (it == mDofStatus.GetNumDependentDofsMap().end())
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] There are no " + Node::DofToString(rDofType) + " dofs.");

    return it->second;
}

int NuTo::StructureBase::GetNumDofs(std::string rDofType) const
{
    return GetNumDofs(Node::DofToEnum(rDofType));
}

int NuTo::StructureBase::GetNumActiveDofs(std::string rDofType) const
{
    return GetNumActiveDofs(Node::DofToEnum(rDofType));
}

int NuTo::StructureBase::GetNumDependentDofs(std::string rDofType) const
{
    return GetNumDependentDofs(Node::DofToEnum(rDofType));
}


const NuTo::BlockSparseMatrix& NuTo::StructureBase::GetConstraintMatrix() const
{
    return mConstraintMatrix;
}

void NuTo::StructureBase::DofTypeSetIsActive(std::string rDofType, bool rIsActive)
{
    DofTypeSetIsActive(Node::DofToEnum(rDofType), rIsActive);
}

void NuTo::StructureBase::DofTypeSetIsConstitutiveInput(std::string rDofType, bool rIsConstitutiveInput)
{
    DofTypeSetIsConstitutiveInput(Node::DofToEnum(rDofType), rIsConstitutiveInput);
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
        mMIS.clear();
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
//            std::cout << "maximum number of independent sets " << mMIS.size() << std::endl;
//            for (unsigned int count=0; count<mMIS.size(); count++)
//            {
//            	std::cout << "MIS " << count << " with " << mMIS[count].size() << " elements " << std::endl;
//            	for (unsigned int count2=0 ; count2<mMIS[count].size(); count2++)
//            		std::cout << ElementGetId(mMIS[count][count2]) << " ";
//            	std::cout << std::endl;
//            }

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

bool NuTo::StructureBase::InterpolationTypeIsConstitutiveInput(NuTo::Node::eDof rDofType)
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
