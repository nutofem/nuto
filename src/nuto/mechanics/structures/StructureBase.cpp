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
#include "nuto/visualize/VisualizeComponentNonlocalWeight.h"
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
        mLogger<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
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
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    VisualizeUnstructuredGrid visualize;
    this->DefineVisualizeData(visualize,mVisualizeComponents);
    this->ElementTotalAddToVisualize(visualize,mVisualizeComponents);
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
    this->DefineVisualizeData(visualize,mVisualizeComponents);
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
    this->DefineVisualizeData(visualize,mVisualizeComponents);
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

void NuTo::StructureBase::DefineVisualizeData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const
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
        default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] invalid data description.");
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
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRGeneral<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK);
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);
        if (error!=Error::SUCCESSFUL)
            return error;

        // build global matrix
        SparseMatrixCSRGeneral<double> transConstraintMatrix = this->mConstraintMatrix.Transpose();
        rMatrix -= transConstraintMatrix * coefficientMatrixKJ + coefficientMatrixJK * this->mConstraintMatrix;
        rMatrix += transConstraintMatrix * coefficientMatrixKK * this->mConstraintMatrix;

        // build equivalent load vector
        rVector = (transConstraintMatrix * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues);
    }
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
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRSymmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0Symmetric(rMatrix, coefficientMatrixJK);
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0Symmetric(rMatrix, coefficientMatrixJK, coefficientMatrixKK);
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
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2General<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK);
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);
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
    return Error::SUCCESSFUL;
}

// build global coefficient matrix0
NuTo::Error::eError NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK);
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
        Error::eError error = this->BuildGlobalCoefficientSubMatrices0Symmetric(rMatrix, coefficientMatrixJK, coefficientMatrixKK);
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
    return Error::SUCCESSFUL;
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

                //mLogger << "no convergence with current step size (" << deltaLoadFactor << "), current not converging load factor " << curLoadFactor << "\n";
                //mLogger << "check stiffness " << "\n";
                //CheckStiffness();
                //mLogger << "and continue with smaller load step " << "\n";

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
