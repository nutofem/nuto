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

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

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
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include "nuto/visualize/VisualizeComponentConstitutive.h"
#include "nuto/visualize/VisualizeComponentCrack.h"
#include "nuto/visualize/VisualizeComponentDamage.h"
#include "nuto/visualize/VisualizeComponentDisplacement.h"
#include "nuto/visualize/VisualizeComponentElement.h"
#include "nuto/visualize/VisualizeComponentEngineeringPlasticStrain.h"
#include "nuto/visualize/VisualizeComponentEngineeringStrain.h"
#include "nuto/visualize/VisualizeComponentEngineeringStress.h"
#include "nuto/visualize/VisualizeComponentLatticeStrain.h"
#include "nuto/visualize/VisualizeComponentLatticeStress.h"
#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
#include "nuto/visualize/VisualizeComponentParticleRadius.h"
#include "nuto/visualize/VisualizeComponentPrincipalEngineeringStress.h"
#include "nuto/visualize/VisualizeComponentSection.h"
#endif // ENABLE_VISUALIZE


NuTo::StructureBase::StructureBase(int rDimension)  : NuTo::NuToObject::NuToObject()
{
    if (rDimension!=1 && rDimension!=2 && rDimension!=3)
    {
        throw MechanicsException("[StructureBase::StructureBase] The dimension of a structure is either 1, 2 or 3.");
    }
    mDimension = rDimension;
    mNumActiveDofs = 0;
    mNumDofs   = 0;
    mNodeNumberingRequired = true;

    mMappingIntEnum2String.resize(NuTo::IntegrationType::NumIntegrationTypes);
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss1Ip]=
        NuTo::IntegrationType1D2NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss2Ip]=
        NuTo::IntegrationType1D2NGauss2Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType1D2NGauss3Ip]=
        NuTo::IntegrationType1D2NGauss3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss1Ip]=
        NuTo::IntegrationType2D3NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D3NGauss3Ip]=
        NuTo::IntegrationType2D3NGauss3Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NGauss1Ip]=
        NuTo::IntegrationType2D4NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType2D4NGauss4Ip]=
        NuTo::IntegrationType2D4NGauss4Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NGauss1Ip]=
        NuTo::IntegrationType3D8NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip]=
        NuTo::IntegrationType3D8NGauss2x2x2Ip::GetStrIdentifierStatic();

    mHaveTmpStaticData = false;
    mUpdateTmpStaticDataRequired = true;
    mToleranceStiffnessEntries = 0.;
    //parameters for Newton Raphson iteration
    mToleranceResidualForce = 1e-6;
    mAutomaticLoadstepControl=true;
    mMaxDeltaLoadFactor=1.;
    mMaxNumNewtonIterations=20;
    mDecreaseFactor=0.5;
    mMinNumNewtonIterations=7;
    mIncreaseFactor=1.5;
    mMinDeltaLoadFactor=1e-6;
    mMinLineSearchFactor=1e-3;

#ifdef _OPENMP
    mUseMIS = true;
    // then the environment variable is used
    mNumProcessors = 1;
#endif // _OPENMP

}

int NuTo::StructureBase::GetDimension()
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
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
    & BOOST_SERIALIZATION_NVP(mDimension)
    & BOOST_SERIALIZATION_NVP(mConstitutiveLawMap)
    & BOOST_SERIALIZATION_NVP(mConstraintMap)
    & BOOST_SERIALIZATION_NVP(mLoadMap)
    & BOOST_SERIALIZATION_NVP(mGroupMap)
    & BOOST_SERIALIZATION_NVP(mIntegrationTypeMap)
    & BOOST_SERIALIZATION_NVP(mSectionMap)
    & BOOST_SERIALIZATION_NVP(mMappingIntEnum2String)
#ifdef ENABLE_VISUALIZE
    & BOOST_SERIALIZATION_NVP(mVisualizeComponents)
#endif //  ENABLE_VISUALIZE
    & BOOST_SERIALIZATION_NVP(mNumDofs)
    & BOOST_SERIALIZATION_NVP(mNumActiveDofs)
    & BOOST_SERIALIZATION_NVP(mNodeNumberingRequired)
    & BOOST_SERIALIZATION_NVP(mConstraintMatrix)
    & BOOST_SERIALIZATION_NVP(mConstraintRHS)
    & BOOST_SERIALIZATION_NVP(mHaveTmpStaticData)
    & BOOST_SERIALIZATION_NVP(mUpdateTmpStaticDataRequired)
    & BOOST_SERIALIZATION_NVP(mToleranceStiffnessEntries)
    & BOOST_SERIALIZATION_NVP(mToleranceResidualForce)
    & BOOST_SERIALIZATION_NVP(mAutomaticLoadstepControl)
    & BOOST_SERIALIZATION_NVP(mMaxDeltaLoadFactor)
    & BOOST_SERIALIZATION_NVP(mMaxNumNewtonIterations)
    & BOOST_SERIALIZATION_NVP(mDecreaseFactor)
    & BOOST_SERIALIZATION_NVP(mMinNumNewtonIterations)
    & BOOST_SERIALIZATION_NVP(mIncreaseFactor)
    & BOOST_SERIALIZATION_NVP(mMinDeltaLoadFactor)
    & BOOST_SERIALIZATION_NVP(mMinLineSearchFactor)
#ifdef _OPENMP
    & BOOST_SERIALIZATION_NVP(mUseMIS)
    & BOOST_SERIALIZATION_NVP(mMIS)
    & BOOST_SERIALIZATION_NVP(mNumProcessors)
#endif // _OPENMP
    & BOOST_SERIALIZATION_NVP(mLogger);
#ifdef DEBUG_SERIALIZATION
    mLogger << "finish serialization of structure base" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureBase::Info()const
{
    mLogger << "dimension : " << mDimension << "\n";

    mLogger << "num dofs : " << mNumDofs << "\n";
    mLogger << "num active dofs : " << mNumActiveDofs << "\n";
    // print info for sections
    SectionInfo(mVerboseLevel);

    // print info for groups
    GroupInfo(mVerboseLevel);
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

// Export to Vtk data file ////////////////////////////////////////////////////
#ifdef ENABLE_VISUALIZE
//! @brief ... Add visualization displacements to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentDisplacements()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentDisplacement());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add element ID to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentElement()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentElement());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentElement] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add engineering strains to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringStrain()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringStrain());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add engineering plastic strain stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringPlasticStrain()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringPlasticStrain());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add engineering stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringStress()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringStress());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add section to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentSection()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentSection());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}
//! @brief ... Add constitutive id to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentConstitutive()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentConstitutive());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add nonlocal weights to the internal list, which is finally exported via the ExportVtkDataFile command
//! @param rElementId ... Element id
//! @param rIp ... local ip number
void NuTo::StructureBase::AddVisualizationComponentNonlocalWeights(int rElementId, int rIp)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        const ElementBase *elementBase = ElementGetElementPtr(rElementId);
        int numIp = elementBase->GetNumIntegrationPoints();
        if (rIp<0 || rIp>=numIp)
            throw MechanicsException("[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] Integration point number is out of range.");
        mVisualizeComponents.push_back(new NuTo::VisualizeComponentNonlocalWeight(elementBase,rElementId,rIp));
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] error setting element and local ip number.");
        throw e;
    }
    catch(...)
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] error setting element and local ip number.");
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add the damage variable to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentDamage()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentDamage());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentDamage] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add visualization of principal stresses to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentPrincipalEngineeringStress()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentPrincipalEngineeringStress());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentPrincipalEngineeringStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add crack id vector to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentCracks()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentCrack());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentCracks] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add visualization particle radius to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentParticleRadius()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentParticleRadius());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentParticleRadius] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add lattice stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentLatticeStress()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentLatticeStress());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentLatticeStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//! @brief ... Add lattice stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentLatticeStrain()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.push_back(new NuTo::VisualizeComponentLatticeStrain());
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentLatticeStrain] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}
void NuTo::StructureBase::ClearVisualizationComponents()
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    mVisualizeComponents.clear();
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ClearVisualizationComponents] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


void NuTo::StructureBase::ExportVtkDataFile(const std::string& rFileName) const
{
    mLogger<<"[NuTo::StructureBase::ExportVtkDataFile] this routine is deprecated, use ExportVtkDataFileElements instead." << "\n";
}

void NuTo::StructureBase::ExportVtkDataFileNodes(const std::string& rFileName) const
{
	ExportVtkDataFileNodes(rFileName, false);
}

void NuTo::StructureBase::ExportVtkDataFileElements(const std::string& rFileName) const
{
	ExportVtkDataFileElements(rFileName, false);
}

void NuTo::StructureBase::ExportVtkDataFileNodes(const std::string& rFileName, bool rXML) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeNodeData(visualize,mVisualizeComponents);
    this->NodeTotalAddToVisualize(visualize,mVisualizeComponents);
    if (rXML)
        visualize.ExportVtuDataFile(rFileName);
    else
        visualize.ExportVtkDataFile(rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

void NuTo::StructureBase::ExportVtkDataFileElements(const std::string& rFileName, bool rXML) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeElementData(visualize,mVisualizeComponents);
    this->ElementTotalAddToVisualize(visualize,mVisualizeComponents);
    if (rXML)
        visualize.ExportVtuDataFile(rFileName);
    else
        visualize.ExportVtkDataFile(rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rFileName) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeElementData(visualize,mVisualizeComponents);
    this->ElementGroupAddToVisualize(rGroupIdent,visualize,mVisualizeComponents);
    visualize.ExportVtkDataFile(rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementGroupExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//Visualize for all integration points the fine scale structure(either damage or homogeneous part)
void NuTo::StructureBase::ElementGroupVisualizeIpMultiscale(int rGroupIdent, const std::string& rFileName, bool rVisualizeDamage)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeElementData(visualize,mVisualizeComponents);
    this->ElementGroupAddToVisualizeIpMultiscale(rGroupIdent,visualize,mVisualizeComponents,rVisualizeDamage);
    visualize.ExportVtkDataFile(rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::ElementGroupVisualizeIpMultiscale] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}

//Visualize for all integration points the fine scale structure
void NuTo::StructureBase::ElementGroupVisualizeIpMultiscaleDamage(int rGroupIdent, const std::string& rFileName)const
{
    ElementGroupVisualizeIpMultiscale(rGroupIdent, rFileName, true);
}

//Visualize for all integration points the fine scale structure
void NuTo::StructureBase::ElementGroupVisualizeIpMultiscaleHomogeneous(int rGroupIdent, const std::string& rFileName)const
{
    ElementGroupVisualizeIpMultiscale(rGroupIdent, rFileName, false);
}

void NuTo::StructureBase::DefineVisualizeElementData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const
{
    boost::ptr_list<NuTo::VisualizeComponentBase>::const_iterator itWhat = mVisualizeComponents.begin();
    while (itWhat != mVisualizeComponents.end())
    {
        switch (itWhat->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DAMAGE:
            rVisualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
            rVisualize.DefinePointDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ELEMENT:
            rVisualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
            rVisualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            rVisualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
            rVisualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
            rVisualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::SECTION:
            rVisualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::CONSTITUTIVE:
            rVisualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::CRACK:
            rVisualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
            rVisualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::LATTICE_STRESS:
            rVisualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::LATTICE_STRAIN:
            rVisualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::LATTICE_PLASTIC_STRAIN:
            rVisualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::PARTICLE_RADIUS:
            //do nothing;
            break;
        default:
        	throw MechanicsException("[NuTo::StructureBase::DefineVisualizeElementData] undefined visualize components.");
        }
        itWhat++;
    }
}

void NuTo::StructureBase::DefineVisualizeNodeData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const
{
    boost::ptr_list<NuTo::VisualizeComponentBase>::const_iterator itWhat = mVisualizeComponents.begin();
    while (itWhat != mVisualizeComponents.end())
    {
        switch (itWhat->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
            rVisualize.DefinePointDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::PARTICLE_RADIUS:
        	rVisualize.DefinePointDataScalar(itWhat->GetComponentName());
            break;
        default:
        	;
         }
        itWhat++;
    }
}


#endif // ENABLE_VISUALIZE

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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, SparseMatrixCSRGeneral<double>& rMatrix, FullMatrix<double>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullMatrix<double> activeDofValues;
    FullMatrix<double> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(activeDofValues, dependentDofValues);
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

    rVector.Resize(this->mNumActiveDofs, 1);
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
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, SparseMatrixCSRSymmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullMatrix<double> activeDofValues;
    FullMatrix<double> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(activeDofValues, dependentDofValues);
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
    rVector.Resize(this->mNumActiveDofs, 1);
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
        FullMatrix<double> deltaRHS = this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues;
        FullMatrix<double> Kdd_Mult_DeltaRHS = coefficientMatrixKK * deltaRHS;
        rVector = this->mConstraintMatrix.TransMult(Kdd_Mult_DeltaRHS) - coefficientMatrixJK * deltaRHS;
    }
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, SparseMatrixCSRVector2General<double>& rMatrix, FullMatrix<double>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullMatrix<double> activeDofValues;
    FullMatrix<double> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(activeDofValues, dependentDofValues);
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

    rVector.Resize(this->mNumActiveDofs, 1);
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
                (NuTo::FullMatrix<double>(coefficientMatrixJK)).Info(12,3);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
    //check for dof numbering and build of tmp static data
    BuildGlobalCoefficientMatrixCheck();

    // get dof values stored at the nodes
    FullMatrix<double> activeDofValues;
    FullMatrix<double> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(activeDofValues, dependentDofValues);
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

    rVector.Resize(this->mNumActiveDofs, 1);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::STIFFNESS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::STIFFNESS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::STIFFNESS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::STIFFNESS, rMatrix, rVector);
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


//! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
//! @param rMatrix ... global coefficient matrix (nonsymmetric)
//! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::MASS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::MASS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::MASS, rMatrix, rVector);
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix2(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    Error::eError error = BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::MASS, rMatrix, rVector);
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
void NuTo::StructureBase::BuildGlobalExternalLoadVector(NuTo::FullMatrix<double>& rVector)
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

    rVector.Resize(this->mNumActiveDofs, 1);
    FullMatrix<double> dependentDofLoadVector(this->mNumDofs - this->mNumActiveDofs,1);

    // loop over all loads
    boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
    while (loadIter != this->mLoadMap.end())
    {
        loadIter->second->AddLoadToGlobalSubVectors(rVector, dependentDofLoadVector);
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
void NuTo::StructureBase::BuildGlobalExternalLoadVector(NuTo::FullMatrix<double>& rVector_j, NuTo::FullMatrix<double>& rVector_k)
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

    rVector_j.Resize(this->mNumActiveDofs, 1);
    rVector_k.Resize(this->mNumDofs - this->mNumActiveDofs,1);

    // loop over all loads
    boost::ptr_map<int,LoadBase>::const_iterator loadIter = this->mLoadMap.begin();
    while (loadIter != this->mLoadMap.end())
    {
        loadIter->second->AddLoadToGlobalSubVectors(rVector_j, rVector_k);
        loadIter++;
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::BuildGlobalExternalLoadVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
}


// build global gradient of the internal potential (e.g. the internal forces)
NuTo::Error::eError NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector(NuTo::FullMatrix<double>& rVector)
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


    rVector.Resize(this->mNumActiveDofs, 1);
    FullMatrix<double> dependentDofGradientVector(this->mNumDofs - this->mNumActiveDofs,1);

    try
    {
        // build sub vectors
        Error::eError error = this->BuildGlobalGradientInternalPotentialSubVectors(rVector, dependentDofGradientVector);
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

//! @brief returns the a reference to the constraint matrix
const NuTo::SparseMatrixCSRGeneral<double>& NuTo::StructureBase::GetConstraintMatrix()const
{
    return mConstraintMatrix;
}

//! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
//! @param load factor
void NuTo::StructureBase::SetLoadFactor(double rLoadFactor)
{
    throw MechanicsException("[NuTo::StructureBase::SetLoadFactor] not implemented - overload this function in your derived Structure class.");
}

//! @brief do a postprocessing step after each updated load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataAfterUpdate(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const
{
}

//! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const
{
}

//! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, const NuTo::FullMatrix<double>& rResidualVector)const
{
}

//! @brief do a postprocessing step in each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataInLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, double rPrevResidual)const
{
}

//! @brief initialize some stuff before a new load step (e.g. output directories for visualization, if required)
void NuTo::StructureBase::InitBeforeNewLoadStep(int rLoadStep)
{
}

//! @brief performs a Newton Raphson iteration (displacement and/or load control) - structure is not saved before update
NuTo::Error::eError NuTo::StructureBase::NewtonRaphson()
{
    std::stringstream saveStringStream;
    bool saveStructureBeforeUpdate(false);
    bool isSaved(false);
    bool initialStateInEquilibrium(true);
    return NewtonRaphson(saveStructureBeforeUpdate,saveStringStream,isSaved,initialStateInEquilibrium);
}

//! @brief performs a Newton Raphson iteration (displacement and/or load control)
//! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
//!             be careful, store it only once
NuTo::Error::eError NuTo::StructureBase::NewtonRaphson(bool rSaveStructureBeforeUpdate,
        std::stringstream& rSaveStringStream,
        bool& rIsSaved)
{
    bool initialStateInEquilibrium(true);
    return NewtonRaphson(rSaveStructureBeforeUpdate,rSaveStringStream,rIsSaved,initialStateInEquilibrium);
}

//! @brief performs a Newton Raphson iteration (displacement and/or load control)
//! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
//!             be careful, store it only once
NuTo::Error::eError NuTo::StructureBase::NewtonRaphson(bool rSaveStructureBeforeUpdate,
        std::stringstream& rSaveStringStream,
        bool& rIsSaved, bool rInitialStateInEquilibrium)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    try
    {
        // start analysis
        double deltaLoadFactor(mMaxDeltaLoadFactor);
        double curLoadFactor(0);
        int loadStep(1);

        this->SetLoadFactor(0);
        this->NodeBuildGlobalDofs();
        this->ElementTotalUpdateTmpStaticData();

        NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
        this->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);

        InitBeforeNewLoadStep(loadStep);
        if (mNumActiveDofs==0)
        {
            bool convergenceTotal = false;
            while (convergenceTotal==false)
            {
                try
                {
                    curLoadFactor+=deltaLoadFactor;
                    InitBeforeNewLoadStep(loadStep);
                    this->SetLoadFactor(curLoadFactor);
                    if (curLoadFactor>1)
                    {
                        deltaLoadFactor = 1.-(curLoadFactor-deltaLoadFactor);
                        curLoadFactor = 1.;
                    }
                    this->NodeBuildGlobalDofs();
                    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
                    this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                    Error::eError error = this->ElementTotalUpdateTmpStaticData();
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (error==Error::NO_CONVERGENCE)
                        {
                            mLogger << "No convergence for current load factor  " << curLoadFactor << "\n";
                            curLoadFactor-=deltaLoadFactor;
                            //decrease load step
                            deltaLoadFactor*=mDecreaseFactor;

                            //check for minimum delta (this mostly indicates an error in the software
                            if (deltaLoadFactor<mMinDeltaLoadFactor)
                            {
                                mLogger << "No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                                return Error::NO_CONVERGENCE;
                            }
                            continue;
                        }
                        else
                            return error;
                    }
                    AdaptModel();
                    this->PostProcessDataAfterConvergence(loadStep, 0, curLoadFactor, deltaLoadFactor, 0);
                    error = this->ElementTotalUpdateStaticData();
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (error==Error::NO_CONVERGENCE)
                        {
                            mLogger << "No convergence for current load factor  " << curLoadFactor << "\n";
                            curLoadFactor-=deltaLoadFactor;
                            //decrease load step
                            deltaLoadFactor*=mDecreaseFactor;

                            //check for minimum delta (this mostly indicates an error in the software
                            if (deltaLoadFactor<mMinDeltaLoadFactor)
                            {
                                mLogger << "No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                                return Error::NO_CONVERGENCE;
                            }
                        }
                        else
                            return error;
                        continue;
                    }
                    this->PostProcessDataAfterUpdate(loadStep, 0, curLoadFactor, deltaLoadFactor, 0);
                    loadStep++;
                }
                catch(MechanicsException& e)
                {
                    e.AddMessage("[NuTo::StructureBase::NewtonRaphson] Error in Newton-Raphson iteration.");
                    throw e;
                }
                if (curLoadFactor>1-1e-8)
                    convergenceTotal = true;
            }
            return Error::SUCCESSFUL;
        }

        //init some auxiliary variables
        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
        NuTo::FullMatrix<double> dispForceVector;
        NuTo::FullMatrix<double> intForceVector;
        NuTo::FullMatrix<double> extForceVector;
        NuTo::FullMatrix<double> rhsVector;


        NuTo::FullMatrix<double> intForceVectorInit;
        //calculate the initial out of balance force
        if (rInitialStateInEquilibrium==false)
        {
            try
            {
                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                        mLogger << "[NuTo::StructureBase::NewtonRaphson] exception due to no convergence in initial unequilibrated state " << "\n";
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] exception due to error in initial unequilibrated state");
                }
                error = this->BuildGlobalGradientInternalPotentialVector(intForceVectorInit);
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                        mLogger << "[NuTo::StructureBase::NewtonRaphson] exception due to no convergence in initial unequilibrated state " << "\n";
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] exception due to error in initial unequilibrated state");
                }

                mLogger << "out of balance force (norm) for zero load is " << intForceVectorInit.Norm() << "\n";
            }
            catch(MechanicsException& e)
            {
                e.AddMessage("[NuTo::StructureBase::NewtonRaphson] [NuTo::StructureBase::NewtonRaphson] exception for constitutive model in unloaded state.");
                throw e;
            }
            catch(...)
            {
                throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] exception for constitutive model in unloaded state.");
            }
        }

        //allocate solver
        NuTo::SparseDirectSolverMUMPS mySolver;
        //NuTo::SparseDirectSolverMKLPardiso mySolver;
#ifdef SHOW_TIME
        if (mShowTime==true)
            mySolver.SetShowTime(true);
        else
            mySolver.SetShowTime(false);
#endif
        bool convergenceInitialLoadStep(false);
        while (convergenceInitialLoadStep==false)
        {
            try
            {
                //calculate stiffness
                curLoadFactor=deltaLoadFactor;
                this->SetLoadFactor(curLoadFactor);
                this->NodeBuildGlobalDofs();
                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                    {
                        //decrease load step
                        deltaLoadFactor*=mDecreaseFactor;

                        //restore initial state
                        this->SetLoadFactor(0);
                        this->NodeBuildGlobalDofs();
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);

                        //check for minimum delta (this mostly indicates an error in the software
                        if (deltaLoadFactor<mMinDeltaLoadFactor)
                        {
                            mLogger << "[NuTo::StructureBase::NewtonRaphson] No convergence for ElementTotalUpdateTmpStaticData calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                            return Error::NO_CONVERGENCE;
                        }
                        continue;
                    }
                    else
                        return error;
                }

                error = this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                    {
                        //decrease load step
                        deltaLoadFactor*=mDecreaseFactor;

                        //restore initial state
                        this->SetLoadFactor(0);
                        this->NodeBuildGlobalDofs();
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);

                        //check for minimum delta (this mostly indicates an error in the software
                        if (deltaLoadFactor<mMinDeltaLoadFactor)
                        {
                            mLogger << "[NuTo::StructureBase::NewtonRaphson] No convergence for BuildGlobalCoefficientMatrix0 calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                            return Error::NO_CONVERGENCE;
                        }
                        continue;
                    }
                    else
                        return error;
                }
                //   NuTo::FullMatrix<double>(stiffnessMatrixCSRVector2).Info(12,3);
                //    mLogger << "disp force vector "<< "\n";
                //    dispForceVector.Trans().Info(12,10);
                //Check the stiffness matrix
                //CheckStiffness();
                //mLogger << "stiffness is calculated in Newton Raphson " << "\n";
                //mLogger << "total energy of system " << ElementTotalGetTotalEnergy() << "\n";

                //update displacements of all nodes according to the new conre mat
                {
                    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
                    this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                    error = this->ElementTotalUpdateTmpStaticData();
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (error==Error::NO_CONVERGENCE)
                        {
                            //decrease load step
                            deltaLoadFactor*=mDecreaseFactor;

                            //restore initial state
                            this->SetLoadFactor(0);
                            this->NodeBuildGlobalDofs();
                            this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);

                            //check for minimum delta (this mostly indicates an error in the software
                            if (deltaLoadFactor<mMinDeltaLoadFactor)
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson] return with no convergence" << "\n";
                                return Error::NO_CONVERGENCE;
                            }
                            continue;
                        }
                        else
                            return error;
                    }
                }

                // build global external load vector and RHS vector
                this->BuildGlobalExternalLoadVector(extForceVector);
                //mLogger<<" calculate gradient 1163" << "\n";
                error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                    {
                        //decrease load step
                        deltaLoadFactor*=mDecreaseFactor;

                        //restore initial state
                        this->SetLoadFactor(0);
                        this->NodeBuildGlobalDofs();
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);

                        //check for minimum delta (this mostly indicates an error in the software
                        if (deltaLoadFactor<mMinDeltaLoadFactor)
                        {
                            mLogger << "No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                            return Error::NO_CONVERGENCE;
                        }
                        continue;
                    }
                    else
                        return error;
                }
                convergenceInitialLoadStep = true;
            }
            catch(MechanicsException& e)
            {
                e.AddMessage("[NuTo::StructureBase::NewtonRaphson] Error in Newton-Raphson iteration.");
                throw e;
            }
        }

        rhsVector = extForceVector - intForceVector;
        //attention this is only different for the first iteration step
        //since the internal force due to the applied constraints is not considered for the first iteration
        //in order to balance it (no localization in the boundary region)
        //for the linesearch this internal force has to be considered in order to obtain for a linesearch
        //factor of zero the normRHS
        double normRHS = rhsVector.Norm();
//    rhsVector.Trans().Info(12,10);
        rhsVector = extForceVector + dispForceVector;
        if (rInitialStateInEquilibrium==false)
        {
            rhsVector -= intForceVectorInit;
        }
//    rhsVector.Trans().Info(12,10);

        //calculate absolute tolerance for matrix entries to be not considered as zero
        double maxValue, minValue, ToleranceZeroStiffness;
        if (stiffnessMatrixCSRVector2.GetNumColumns()==0)
        {
            maxValue = 1.;
            minValue = 1.;
        }
        else
        {
            stiffnessMatrixCSRVector2.Max(maxValue);
            stiffnessMatrixCSRVector2.Min(minValue);
        }
        //mLogger << "min and max " << minValue << " , " << maxValue << "\n";

        ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
        this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
        //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
        //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
        //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

        //mySolver.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale") + std::string("0") + std::string(".vtk"));

        //store the structure only once in order to be able to restore the situation before entering the routine
        rIsSaved = false;

        //repeat until max displacement is reached
        bool convergenceStatusLoadSteps(false);
        while (!convergenceStatusLoadSteps)
        {
            double normResidual(1);
            double maxResidual(1);
            int numNewtonIterations(0);
            double alpha(1.);
            NuTo::FullMatrix<double> displacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsDependentDOFs;
            int convergenceStatus(0);
            //0 - not converged, continue Newton iteration
            //1 - converged
            //2 - stop iteration, decrease load step
            while(convergenceStatus==0)
            {
                numNewtonIterations++;

                if (numNewtonIterations>mMaxNumNewtonIterations && alpha<0.25)
                {
                    if (mVerboseLevel>5)
                    {
                        mLogger << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << mMaxNumNewtonIterations << ")" << "\n";
                    }
                    convergenceStatus = 2; //decrease load step
                    break;
                }

                // solve
                NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
                NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
                this->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);
                NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
                stiffnessMatrixCSR.SetOneBasedIndexing();
                try
                {
                    mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
                }
                catch(...)
                {
                    mLogger << "Error solving system of equations using mumps." << "\n";
                    if (mNumActiveDofs<1000)
                    {
                        NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                        if (mNumActiveDofs<30)
                        {
                            mLogger << "stiffness full" << "\n";
                            mLogger.Out(stiffnessMatrixFull,12,3);
                        }
                        NuTo::FullMatrix<double> eigenValues;
                        stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                        mLogger << "eigenvalues" << "\n";
                        mLogger.Out(eigenValues.Trans(),12,3);
                        NuTo::FullMatrix<double> eigenVectors;
                        stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                        mLogger << "eigenvector 1" << "\n";
                        mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
                    }
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] Error solving system of equations using mumps.");
                }

                //mLogger << " rhsVector" << "\n";
                //rhsVector.Trans().Info(10,3);
                //std::cout << " delta_disp" << "\n";
                //deltaDisplacementsActiveDOFs.Trans().Info(10,3);

                //perform a linesearch
                alpha = 1.;
                do
                {
                    //add new displacement state
                    displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;

                    //mLogger << " displacementsActiveDOFs" << "\n";
                    //displacementsActiveDOFs.Trans().Info(10,3);
                    this->NodeMergeActiveDofValues(displacementsActiveDOFs);
                    try
                    {
                        Error::eError error = this->ElementTotalUpdateTmpStaticData();
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                convergenceStatus=2;
                                mLogger << "Constitutive model is not converging, try with smaller load step" << "\n";
                            }
                            else
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson] error while calling the gradient (resforce) routine." << "\n";
                                return error;
                            }
                        }
                        else
                        {
/*if (mNumActiveDofs<1000)
{
    this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
	NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
	if (mNumActiveDofs<30)
	{
		mLogger << "stiffness full" << "\n";
		mLogger.Out(stiffnessMatrixFull,12,3);
	}
	NuTo::FullMatrix<double> eigenValues;
	stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
	mLogger << "eigenvalues" << "\n";
	mLogger.Out(eigenValues.Trans(),12,3);
	NuTo::FullMatrix<double> eigenVectors;
	stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
	mLogger << "eigenvector 1" << "\n";
	mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
}
*/
                        	// calculate residual
                            Error::eError error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                            if (error!=Error::SUCCESSFUL)
                            {
                                if (error==Error::NO_CONVERGENCE)
                                {
                                    convergenceStatus=2;
                                    mLogger << "Constitutive model is not converging, try with smaller load step" << "\n";
                                }
                                else
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson] error while calling the gradient (resforce) routine." << "\n";
                                    return error;
                                }
                            }
                            else
                            {
                                //mLogger << "intForceVector "  << "\n";
                                //intForceVector.Trans().Info(10,3);
                                rhsVector = extForceVector - intForceVector;
                                normResidual = rhsVector.Norm();
                                maxResidual = rhsVector.Abs().Max();

                                //mLogger << "total energy of system " << ElementTotalGetTotalEnergy() << "\n";
                                PostProcessDataInLineSearch(loadStep, numNewtonIterations, alpha, curLoadFactor, normResidual, normRHS);
                            }
                        }
                    }
                    catch(MechanicsException& e)
                    {
                        e.AddMessage("[NuTo::StructureBase::NewtonRaphson] NuTo::MechanicsException while calling the gradient (resforce) routine.");
                        throw e;
                    }
                    catch(...)
                    {
                        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] Error calling the gradient (resforce) routine.");
                    }

                    alpha*=0.5;
                }
                while(convergenceStatus!=2 && alpha>mMinLineSearchFactor && normResidual>normRHS*(1-0.5*alpha) && normResidual>mToleranceResidualForce && maxResidual>mToleranceResidualForce);

                if (convergenceStatus==2)
                    break;

                this->PostProcessDataAfterLineSearch(loadStep, numNewtonIterations, 2.*alpha, curLoadFactor, normResidual, rhsVector);
                //std::string str;
                //getline (std::cin,str);

                if (normResidual>normRHS*(1-0.5*alpha) && normResidual>mToleranceResidualForce && maxResidual>mToleranceResidualForce)
                {
                    convergenceStatus=2;
                    {
                        if (mNumActiveDofs<1000)
                        {
                            mLogger << "System is not converging." << "\n";
                            NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                            if (mNumActiveDofs<30)
                            {
                                mLogger << "stiffness full" << "\n";
                                mLogger.Out(stiffnessMatrixFull,12,3);
                            }
                            NuTo::FullMatrix<double> eigenValues;
                            stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                            mLogger << "eigenvalues" << "\n";
                            mLogger.Out(eigenValues.Trans(),12,3);
                            NuTo::FullMatrix<double> eigenVectors;
                            stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                            mLogger << "eigenvector 1" << "\n";
                            mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
                        }
                    }
                    break;
                }

                //mLogger << "\n" << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<"\n";

                //check convergence
                if (normResidual<mToleranceResidualForce || maxResidual<mToleranceResidualForce)
                {
                    //this test is only relevant for problems with a model adaptation, otherwise, just assume a converged solution and continue
                    if(AdaptModel())
                    {
                        mLogger << "adaptation is performed " << "\n";
                        numNewtonIterations=0;
                        try
                        {
                            Error::eError error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                            if (error!=Error::SUCCESSFUL)
                            {
                                if (error==Error::NO_CONVERGENCE)
                                {
                                    //decrease load step
                                    deltaLoadFactor*=mDecreaseFactor;

                                    //restore initial state
                                    this->SetLoadFactor(0);
                                    this->NodeBuildGlobalDofs();
                                    this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                                    convergenceStatus=2;

                                    mLogger << "********************************************************************************" << "\n";
                                    mLogger << "**************** reduce load step after adaptation is performed ****************" << "\n";
                                    mLogger << "********************************************************************************" << "\n";

                                    //check for minimum delta (this mostly indicates an error in the software
                                    if (deltaLoadFactor<mMinDeltaLoadFactor)
                                    {
                                        mLogger << "[NuTo::StructureBase::NewtonRaphson] No convergence for gradient calculation after adaptation." << "\n";
                                        return Error::NO_CONVERGENCE;
                                    }
                                }
                                else
                                    return error;
                            }
                        }
                        catch(MechanicsException& e)
                        {
                            e.AddMessage("[NuTo::StructureBase::NewtonRaphson] error in gradient calculation after adaptation.");
                            throw e;
                        }

                        //mLogger << "intForceVector "  << "\n";
                        //intForceVector.Trans().Info(10,3);
                        rhsVector = extForceVector - intForceVector;
                        mLogger<<" normRHS after adaptation" << rhsVector.Norm() << "\n";
                        //exit(0);
                    }
                    else
                    {
                        this->PostProcessDataAfterConvergence(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                        convergenceStatus=1;
                        //CheckStiffness();
                        //NodeInfo(12);
                        break;
                    }
                }

                //convergence status == 0 (continue Newton iteration)
                normRHS = rhsVector.Norm();
                //build new stiffness matrix
                //there should not be a problem here, since the gradient has already been successfully caclulated for that loading scenario
                Error::eError error = this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
                if (error!=Error::SUCCESSFUL)
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] stiffness for first step of new load increment could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");

                //mLogger << dispForceVector.Norm() << "\n";
//check stiffness
//CheckStiffness();
                //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
                //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
                //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";
            }

            if (convergenceStatus==1)
            {
                // the update is only required to allow for a stepwise solution procedure in the fine scale model
                // a final update is only required for an update on the macroscale, otherwise,the original state has
                // to be reconstructed.

                if (curLoadFactor>1-1e-8)
                {
                    if (rSaveStructureBeforeUpdate==false)
                    {
                        Error::eError error = this->ElementTotalUpdateStaticData();
                        if (error!=Error::SUCCESSFUL)
                            throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] update could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");
                        this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                    }
                    convergenceStatusLoadSteps=true;
                }
                else
                {
                    Error::eError error = this->ElementTotalUpdateStaticData();
                    if (error!=Error::SUCCESSFUL)
                        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] update could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");

                    this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                    //store the last converged step in order to be able to go back to that state
                    displacementsActiveDOFsLastConverged  = displacementsActiveDOFs;

                    //eventually increase load step
                    if (mAutomaticLoadstepControl)
                    {
                        if (numNewtonIterations<mMinNumNewtonIterations)
                        {
                            deltaLoadFactor*=mIncreaseFactor;
                        }
                        if (deltaLoadFactor>mMaxDeltaLoadFactor)
                            deltaLoadFactor = mMaxDeltaLoadFactor;
                    }

                    //increase displacement
                    curLoadFactor+=deltaLoadFactor;
                    if (curLoadFactor>1)
                    {
                        deltaLoadFactor -= curLoadFactor -1.;
                        curLoadFactor=1;
                    }
                }
                loadStep++;
                //initialize some stuff before a new load step (e.g. output directories for visualization, if required)
                InitBeforeNewLoadStep(loadStep);
            }
            else
            {
                assert(convergenceStatus==2);
                if (mAutomaticLoadstepControl==false)
                    return Error::NO_CONVERGENCE;

                mLogger << "no convergence with current step size (" << deltaLoadFactor << "), current not converging load factor " << curLoadFactor << "\n";
                mLogger << "check stiffness " << "\n";
                CheckStiffness();
                mLogger << "and continue with smaller load step " << "\n";

                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                curLoadFactor-=deltaLoadFactor;

                //set the previous displacement state
                this->SetLoadFactor(curLoadFactor);

                // build global dof numbering
                this->NodeBuildGlobalDofs();

                //set previous converged displacements
                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                if (error!=Error::SUCCESSFUL)
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] last converged state could not be recalculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");
                //this first part of the routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
                //for an update on the coarse scale
                //as a consequence, in an iterative solution with updates in between, the initial state has to be restored after leaving the routine
                if (rSaveStructureBeforeUpdate==true && rIsSaved==false)
                {
                    assert(curLoadFactor==0);
                    //store the structure only once in order to be able to restore the situation before entering the routine
                    this->SaveStructure(rSaveStringStream);
                    rIsSaved = true;
                }

                //decrease load step
                deltaLoadFactor*=mDecreaseFactor;
                curLoadFactor+=deltaLoadFactor;

                //check for minimum delta (this mostly indicates an error in the software)
                if (deltaLoadFactor<mMinDeltaLoadFactor)
                {
                    mLogger << "return with a MechanicsNoConvergenceException " << "\n";
                    return Error::NO_CONVERGENCE;
                }
            }

            if (!convergenceStatusLoadSteps)
            {
                try
                {
                    bool convergenceConstitutive(false);
                    while (convergenceConstitutive==false)
                    {
                        this->SetLoadFactor(curLoadFactor);
                        this->NodeBuildGlobalDofs();

                        //update stiffness in order to calculate new dispForceVector, but still with previous displacement state
                        Error::eError error = this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (mAutomaticLoadstepControl==false)
                                return Error::NO_CONVERGENCE;
                            mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                            curLoadFactor-=deltaLoadFactor;

                            //set the previous displacement state
                            this->SetLoadFactor(curLoadFactor);

                            // build global dof numbering
                            this->NodeBuildGlobalDofs();

                            //set previous converged displacements
                            this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                            Error::eError error = this->ElementTotalUpdateTmpStaticData();
                            if (error!=Error::SUCCESSFUL)
                                throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] previous converged load step could not be recalculated, check your implementation.");

                            //decrease load step
                            deltaLoadFactor*=mDecreaseFactor;
                            curLoadFactor+=deltaLoadFactor;

                            //check for minimum delta (this mostly indicates an error in the software
                            if (deltaLoadFactor<mMinDeltaLoadFactor)
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                return Error::NO_CONVERGENCE;
                            }
                            continue;
                        }
                        //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
                        //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
                        //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

                        //update displacements of all nodes according to the new conre mat
                        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
                        this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                        error = this->ElementTotalUpdateTmpStaticData();
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                if (mAutomaticLoadstepControl==false)
                                    return Error::NO_CONVERGENCE;
                                mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                                curLoadFactor-=deltaLoadFactor;

                                //set the previous displacement state
                                this->SetLoadFactor(curLoadFactor);

                                // build global dof numbering
                                this->NodeBuildGlobalDofs();

                                //set previous converged displacements
                                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                                if (error!=Error::SUCCESSFUL)
                                    return error;

                                //decrease load step
                                deltaLoadFactor*=mDecreaseFactor;
                                curLoadFactor+=deltaLoadFactor;

                                //check for minimum delta (this mostly indicates an error in the software
                                if (deltaLoadFactor<mMinDeltaLoadFactor)
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                    mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                    return Error::NO_CONVERGENCE;
                                }
                                continue;
                            }
                            else
                                return error;
                        }

                        // calculate initial residual for next load step
                        //mLogger<<" calculate gradient 1532" << "\n";
                        error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                if (mAutomaticLoadstepControl==false)
                                    return Error::NO_CONVERGENCE;
                                mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                                curLoadFactor-=deltaLoadFactor;

                                //set the previous displacement state
                                this->SetLoadFactor(curLoadFactor);

                                // build global dof numbering
                                this->NodeBuildGlobalDofs();

                                //set previous converged displacements
                                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                                if (error!=Error::SUCCESSFUL)
                                    return error;

                                //decrease load step
                                deltaLoadFactor*=mDecreaseFactor;
                                curLoadFactor+=deltaLoadFactor;

                                //check for minimum delta (this mostly indicates an error in the software
                                if (deltaLoadFactor<mMinDeltaLoadFactor)
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                    mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                    return Error::NO_CONVERGENCE;
                                }
                                continue;
                            }
                            else
                                return error;
                        }

                        convergenceConstitutive = true;
                    }

                }
                catch(NuTo::MechanicsException& e)
                {
                    e.AddMessage("[NuTo::StructureBase::NewtonRaphson]: error calculating new displacement increment.");
                    throw e;
                }
                catch(...)
                {
                    throw NuTo::MechanicsException("[NuTo::StructureBase::NewtonRaphson]: Exception after the initial stiffness calculation after a load decrease.");
                }

                //update rhs vector for next Newton iteration
                rhsVector = extForceVector - intForceVector;
                normRHS = rhsVector.Norm();
                //attention this is only different for the first iteration step (load application)
                //since the internal force due to the applied constraints is not considered for the first iteration
                //in order to balance it (no localization in the boundary region)
                //for the linesearch this internal force has to be considered in order to obtain for a linesearch
                //factor of zero the normRHS
                rhsVector = dispForceVector + extForceVector;
                if (rInitialStateInEquilibrium==false && curLoadFactor==deltaLoadFactor)
                {
                    rhsVector -= intForceVectorInit;
                }
            }
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::NewtonRaphson] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::NewtonRaphson] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mLogger<<"[NuTo::StructureBase::NewtonRaphson] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;
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
        std::map<NodeBase*, std::vector<unsigned int> > elementsPerNode;
        for (unsigned int elementCount = 0; elementCount< elementVector.size(); elementCount++)
        {
            for (int nodeCount = 0; nodeCount< elementVector[elementCount]->GetNumNodes(); nodeCount++)
            {
                elementsPerNode[elementVector[elementCount]->GetNode(nodeCount)].push_back(elementCount);
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
NuTo::FullMatrix<double> NuTo::StructureBase::CreateSpheresOnSpecimenBoundary(int rTypeOfSpecimen,
		FullMatrix<double>& rBoundingBox, int rSeed,
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

	FullMatrix<double> particles;
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
		if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-1.5*D)>1e-10)
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
		if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-D)>1e-10)
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
			if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-1.5*D)>1e-10)
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
			if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-D)>1e-10)
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
//! @parameters rTypeOfSpecimen 0 box, 1 dogbone
//! @parameters rBoundingBox box for the spheres (3*2 matrix)
//! @parameters rRelParticleMass percentage of particle mass inside the box
//! @parameters rGradingCurve matrix with each line min_diameter, max_diameter, mass percentage of that sieve size and density of particles
//! @parameters relativeDistance scaling factor to increase the diameter when inserting the sphere to ensure a minimum distance
//! @parameters rDensity density of the mixture (concrete)
//! @parameters rSeed seed for the random number generator
//! @parameters rSpheresBoundary particles simulated on the boundary e.g. created with CreateSpheresOnBoxBoundary (they do not contribute to the grading curve)
//! @return ... matrix with spheres (coordinates x y z and radius)
NuTo::FullMatrix<double> NuTo::StructureBase::CreateSpheresInSpecimen(int rTypeOfSpecimen, FullMatrix<double>& rBoundingBox, double rRelParticleMass, FullMatrix<double>& rGradingCurve,
		double relativeDistance, double rDensity, int rSeed, NuTo::FullMatrix<double>& rSpheresBoundary)
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
		if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-1.5*D)>1e-10)
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
		if (fabs(rBoundingBox(1,1) - rBoundingBox(1,0)-D)>1e-10)
			throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");
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
	double massSumParticles = Vspecimen * rDensity * rRelParticleMass;

	FullMatrix<double> particles(rSpheresBoundary);

	int numParticles(rSpheresBoundary.GetNumRows());

	// random number generator
	dsfmt_t randomNumberGenerator;
	// init random number generator with milliseconds from ..
	dsfmt_init_gen_rand(&randomNumberGenerator, rSeed);
	double lastPrintedFraction = 0.;
	for(int gc=0;gc<rGradingCurve.GetNumRows();gc++)
	{
		double dMin(rGradingCurve(gc,0));
		double dMax(rGradingCurve(gc,1));
		double massFrac(rGradingCurve(gc,2));
		double particleDensity(rGradingCurve(gc,3));
		numParticlesPerClass[gc]=0;

		// calculate reference volume fraction of the mineral-size-class
		Vsoll[gc] = massSumParticles * massFrac / particleDensity;
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
		std::sort(((double*)&particles.mEigenMatrix.data()[3*particles.GetNumRows()+numParticles-numParticlesPerClass[gc]]),((double*)&particles.mEigenMatrix.data()[3*particles.GetNumRows()+numParticles]), std::greater<double>( ));

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
					double sumR = 0.5*D-particles(countParticle,3);
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
					double sumR = particles(countParticle,3)*(1.+relativeDistance)+particles(otherEllipse,3);
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
					if ((double)(countParticle)/numParticles-lastPrintedFraction>0.01)
					{
						std::cout << (double)(countParticle)/numParticles*100. << "% of particles already inserted." << "\n";
						lastPrintedFraction = (double)(countParticle)/numParticles;
					}
				}
				else
				{
					numTries++;
					if (numTries>100000)
						throw MechanicsException("[NuTo::StructureBase::CreateSpheresInBox] unable to insert sphere after a 100000 tries.");
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
NuTo::FullMatrix<double> NuTo::StructureBase::CutSpheresZ(NuTo::FullMatrix<double>& rSpheres, double rZCoord, double rMinRadius)
{
	NuTo::FullMatrix<double> circles(1000,3);
	int numCircles(0);
	for (int countSphere=0; countSphere<rSpheres.GetNumRows(); countSphere++)
	{
        double delta=rSpheres(countSphere,2)-rZCoord;
		if (fabs(delta)<rSpheres(countSphere,3))
        {
			double radius=sqrt(rSpheres(countSphere,3)*rSpheres(countSphere,3)-delta*delta);
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
void NuTo::StructureBase::InsertParticleIntoBox(NuTo::FullMatrix<double>& rParticles, int rTheParticle, std::vector<std::vector<int > >& rSubBox, std::array<int,3>& rNSubBox,std::array<double,3>& rLSubBox, FullMatrix<double>& rBoundingBox)
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

void NuTo::StructureBase::SetOMPNested(bool rNested)
{
#ifdef _OPENMP
	omp_set_nested(rNested);
#endif //_OPENMP
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
