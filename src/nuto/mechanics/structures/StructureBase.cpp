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

#ifdef SHOW_TIME
    #include <ctime>
#endif

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>


#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
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
    std::cout << "start serialization of structure base" << std::endl;
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
       & BOOST_SERIALIZATION_NVP(mVisualizeComponents)
       & BOOST_SERIALIZATION_NVP(mNumDofs)
       & BOOST_SERIALIZATION_NVP(mNumActiveDofs)
       & BOOST_SERIALIZATION_NVP(mNodeNumberingRequired)
       & BOOST_SERIALIZATION_NVP(mConstraintMatrix)
       & BOOST_SERIALIZATION_NVP(mConstraintRHS)
       & BOOST_SERIALIZATION_NVP(mHaveTmpStaticData)
       & BOOST_SERIALIZATION_NVP(mUpdateTmpStaticDataRequired)
       & BOOST_SERIALIZATION_NVP(mToleranceStiffnessEntries);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structure base" << std::endl;
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureBase::Info()const
{
    std::cout << "dimension : " << mDimension << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentElement] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentNonlocalWeights] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentDamage] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentPrincipalEngineeringStress] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::AddVisualizationComponentCracks] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::ClearVisualizationComponents] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


void NuTo::StructureBase::ExportVtkDataFile(const std::string& rFileName) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsTotal(ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rFileName) const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    // find group by name
    const Group<ElementBase>* ElementGroup = dynamic_cast<const Group<ElementBase>*>( this->GroupGetGroupPtr(rGroupIdent));
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsByGroup(ElementGroup,ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ExportVtkDataFile] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

void NuTo::StructureBase::ExportVtkDataFile(const std::vector<const ElementBase*>& rElements, const std::string& rFileName) const
{
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] First update of tmp static data required.");
    }
    VisualizeUnstructuredGrid Visualize;
    boost::ptr_list<NuTo::VisualizeComponentBase>::const_iterator itWhat = mVisualizeComponents.begin();
    while (itWhat != mVisualizeComponents.end())
    {
        switch (itWhat->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DAMAGE:
            Visualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
            Visualize.DefinePointDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ELEMENT:
            Visualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
            Visualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            Visualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
            Visualize.DefineCellDataTensor(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
            Visualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::SECTION:
            Visualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::CONSTITUTIVE:
            Visualize.DefineCellDataScalar(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::CRACK:
            Visualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
            Visualize.DefineCellDataVector(itWhat->GetComponentName());
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] invalid data description.");
        }
        itWhat++;
    }

    if (mHaveTmpStaticData && mUpdateTmpStaticDataRequired)
    {
    	throw NuTo::MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] Update of tmpStaticData required first.");
    }

    for (unsigned int ElementCount = 0; ElementCount < rElements.size(); ElementCount++)
    {
        rElements[ElementCount]->Visualize(Visualize, mVisualizeComponents);
    }
    Visualize.ExportVtkDataFile(rFileName);
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
void NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRGeneral<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
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

    std::cout << "active " << mNumActiveDofs << " dependent " << mNumDofs-mNumActiveDofs << std::endl;

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
        //std::cout << "non-symmetric, zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK);

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
    }
    else
    {
        //std::cout << "non-symmetric, non-zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRGeneral<double> coefficientMatrixKJ(this->mNumDofs - this->mNumActiveDofs, this->mNumActiveDofs);
        SparseMatrixCSRGeneral<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);

        // build global matrix
        SparseMatrixCSRGeneral<double> transConstraintMatrix = this->mConstraintMatrix.Transpose();
        rMatrix -= transConstraintMatrix * coefficientMatrixKJ + coefficientMatrixJK * this->mConstraintMatrix;
        rMatrix += transConstraintMatrix * coefficientMatrixKK * this->mConstraintMatrix;

        // build equivalent load vector
        rVector = (transConstraintMatrix * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

// build global coefficient matrix0
void NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRSymmetric<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
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
        //std::cout << "symmetric, zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0Symmetric(rMatrix, coefficientMatrixJK);

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
    }
    else
    {
        //std::cout << "symmetric, non-zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRGeneral<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRSymmetric<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0Symmetric(rMatrix, coefficientMatrixJK, coefficientMatrixKK);

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
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

// build global coefficient matrix0
void NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2General<double>& rMatrix, FullMatrix<double>& rVector)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
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
        //std::cout << "non-symmetric, zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK);

        // build equivalent load vector
        rVector = coefficientMatrixJK * (dependentDofValues - this->mConstraintRHS);
/*        std::cout << "dependent " << std::endl;
        dependentDofValues.Trans().Info(12,10);
        std::cout << "mConstraintRHS " << std::endl;
        mConstraintRHS.Trans().Info(12,10);
        std::cout << "delta " << std::endl;
        (dependentDofValues - this->mConstraintRHS).Trans().Info(12,10);
        std::cout << "coefficientMatrixJK" << std::endl;
        (NuTo::FullMatrix<double>(coefficientMatrixJK)).Info(12,3);
*/
    }
    else
    {
    	//std::cout << "non-symmetric, non-zero constraint matrix" << std::endl;

        // define additional submatrix
        SparseMatrixCSRVector2General<double> coefficientMatrixJK(this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);
        SparseMatrixCSRVector2General<double> coefficientMatrixKJ(this->mNumDofs - this->mNumActiveDofs, this->mNumActiveDofs);
        SparseMatrixCSRVector2General<double> coefficientMatrixKK(this->mNumDofs - this->mNumActiveDofs, this->mNumDofs - this->mNumActiveDofs);

        // build submatrices
        this->BuildGlobalCoefficientSubMatrices0General(rMatrix, coefficientMatrixJK, coefficientMatrixKJ, coefficientMatrixKK);

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
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
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
        std::cout<<"[NuTo::StructureBase::BuildGlobalExternalLoadVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

// build global gradient of the internal potential (e.g. the internal forces)
void NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector(NuTo::FullMatrix<double>& rVector)
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
        this->BuildGlobalGradientInternalPotentialSubVectors(rVector, dependentDofGradientVector);
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
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
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
        std::cout<<"[NuTo::StructureBase::GetNumDofs] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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
        std::cout<<"[NuTo::StructureBase::GetNumActiveDofs] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
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

//! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor)const
{
    throw MechanicsException("[NuTo::StructureBase::PostProcessDataAfterConvergence] not implemented - overload this function in your derived Structure class.");
}

//! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
void NuTo::StructureBase::PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor)const
{
    throw MechanicsException("[NuTo::StructureBase::PostProcessDataAfterLineSearch] not implemented - overload this function in your derived Structure class.");
}

//! @brief performs a Newton Raphson iteration (displacement and/or load control) - structure is not saved before update
//! @parameters rToleranceResidualForce  convergence criterion for the norm of the residual force vector
//! @parameters rAutomaticLoadstepControl yes, if the step length should be adapted
//! @parameters rMaxNumNewtonIterations maximum number of iterations per Newton iteration
//! @parameters rDecreaseFactor factor to decrease the load factor in case of no convergence with the prescribed number of Newton iterations
//! @parameters rMinNumNewtonIterations if convergence is achieved in less than rMinNumNewtonIterations, the step length is increased
//! @parameters rIncreaseFactor by this factor
//! @parameters rMinLoadFactor if the load factor is smaller the procedure is assumed to diverge (throwing an exception)
void NuTo::StructureBase::NewtonRaphson(double rToleranceResidualForce,
        bool rAutomaticLoadstepControl,
        double rMaxDeltaLoadFactor,
        int rMaxNumNewtonIterations,
        double rDecreaseFactor,
        int rMinNumNewtonIterations,
        double rIncreaseFactor,
        double rMinDeltaLoadFactor)
{
	std::stringstream saveStringStream;
	bool saveStructureBeforeUpdate(false);
	bool isSaved(false);
	NewtonRaphson(rToleranceResidualForce,rAutomaticLoadstepControl,rMaxDeltaLoadFactor,rMaxNumNewtonIterations,rDecreaseFactor,
			rMinNumNewtonIterations,rIncreaseFactor,rMinDeltaLoadFactor,saveStructureBeforeUpdate,saveStringStream,isSaved);
}

//! @brief performs a Newton Raphson iteration (displacement and/or load control)
//! @parameters rToleranceResidualForce  convergence criterion for the norm of the residual force vector
//! @parameters rAutomaticLoadstepControl yes, if the step length should be adapted
//! @parameters rMaxNumNewtonIterations maximum number of iterations per Newton iteration
//! @parameters rDecreaseFactor factor to decrease the load factor in case of no convergence with the prescribed number of Newton iterations
//! @parameters rMinNumNewtonIterations if convergence is achieved in less than rMinNumNewtonIterations, the step length is increased
//! @parameters rIncreaseFactor by this factor
//! @parameters rMinLoadFactor if the load factor is smaller the procedure is assumed to diverge (throwing an exception)
//! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
//!             be careful, store it only once, although the routine is called before every update
void NuTo::StructureBase::NewtonRaphson(double rToleranceResidualForce,
        bool rAutomaticLoadstepControl,
        double rMaxDeltaLoadFactor,
        int rMaxNumNewtonIterations,
        double rDecreaseFactor,
        int rMinNumNewtonIterations,
        double rIncreaseFactor,
        double rMinDeltaLoadFactor,
        bool rSaveStructureBeforeUpdate,
        std::stringstream& rSaveStringStream,
        bool& rIsSaved)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
try
{
    //check the parameters
    if (rToleranceResidualForce<1e-16)
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] tolerance should be larger than the accuracy of the system.");
    }

    if (rMaxNumNewtonIterations<1)
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] number of Newton iterations should be larger than zero.");
    }

    if (rDecreaseFactor>=1 || rDecreaseFactor<=0)
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] decrease factor should be in the range (0,1).");
    }

    if (rMinNumNewtonIterations<1 || rMinNumNewtonIterations>=rMaxNumNewtonIterations)
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] number of Newton iterations to decide whether the load step can be increased has to be positive and smaller than max Newton iterations.");
    }

    if (rIncreaseFactor<=1 )
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] increase factor should be larger than 1.");
    }

    if (rMinDeltaLoadFactor<=0 )
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] MinDeltaLoadFactor should be positive.");
    }
    if (rMaxDeltaLoadFactor<rMinDeltaLoadFactor )
    {
        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] MaxLoadFactor should equal or greater than minDeltaLoadFactor.");
    }

    //! @brief initializes some variables etc. before the Newton-Raphson routine is executed
    this->InitBeforeNewtonRaphson();

    // start analysis
    double deltaLoadFactor(rMaxDeltaLoadFactor);
    double curLoadFactor(deltaLoadFactor);

    //init some auxiliary variables
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> intForceVector;
    NuTo::FullMatrix<double> extForceVector;
    NuTo::FullMatrix<double> rhsVector;

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
#ifdef SHOW_TIME
   mySolver.SetShowTime(false);
#endif
    //calculate stiffness
    this->SetLoadFactor(curLoadFactor);
    this->NodeBuildGlobalDofs();
    if (mNumActiveDofs==0)
    {
        this->SetLoadFactor(1);
        this->NodeBuildGlobalDofs();
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        this->ElementTotalUpdateTmpStaticData();
        return;
    }
    this->ElementTotalUpdateTmpStaticData();

    this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
//    std::cout << "initial stiffness" << std::endl;
//    NuTo::FullMatrix<double>(stiffnessMatrixCSRVector2).Info(12,3);
//    std::cout << "disp force vector "<< std::endl;
//    dispForceVector.Trans().Info(12,10);
//Check the stiffness matrix
//CheckStiffness();
//std::cout << "total energy of system " << ElementTotalGetTotalEnergy() << std::endl;

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        this->ElementTotalUpdateTmpStaticData();
    }

    // build global external load vector and RHS vector
    this->BuildGlobalExternalLoadVector(extForceVector);
    this->BuildGlobalGradientInternalPotentialVector(intForceVector);
    rhsVector = extForceVector - intForceVector;
    //attention this is only different for the first iteration step
    //since the internal force due to the applied constraints is not considered for the first iteration
    //in order to balance it (no localization in the boundary region)
    //for the linesearch this internal force has to be considered in order to obtain for a linesearch
    //factor of zero the normRHS
    double normRHS = rhsVector.Norm();
//    rhsVector.Trans().Info(12,10);
    rhsVector = extForceVector + dispForceVector;
//    rhsVector.Trans().Info(12,10);

/*    {
        double energy;
        energy = this->ElementTotalGetTotalEnergy();
        energy += this->ConstraintTotalGetTotalEnergy();
        std::cout << "energy " << energy << std::endl;
    }
*/
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
    //std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

    //mySolver.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale") + std::string("0") + std::string(".vtk"));

    //store the structure only once in order to be able to restore the situation before entering the routine
    rIsSaved = false;

    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    int loadStep(1);
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    this->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
    while (!convergenceStatusLoadSteps)
    {
        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
        double alpha(1.);
        int convergenceStatus(0);
        //0 - not converged, continue Newton iteration
        //1 - converged
        //2 - stop iteration, decrease load step
        while(convergenceStatus==0)
        {
            numNewtonIterations++;

            if (numNewtonIterations>rMaxNumNewtonIterations)
            {
                if (mVerboseLevel>5)
                {
                    std::cout << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << rMaxNumNewtonIterations << ")" << std::endl;
                }
                convergenceStatus = 2; //decrease load step
                break;
            }

            // solve
            NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsDependentDOFs;
            this->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);
            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            try
            {
                mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
            }
            catch(...)
            {
                if (mNumActiveDofs<1000)
                {
                    NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                    std::cout << "stiffness full" << std::endl;
                    stiffnessMatrixFull.Info(12,3);
                    NuTo::FullMatrix<double> eigenValues;
                    stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                    std::cout << "eigenvalues" << std::endl;
                    eigenValues.Trans().Info(12,3);
                    NuTo::FullMatrix<double> eigenVectors;
                    stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                    std::cout << "eigenvector 1" << std::endl;
                    eigenVectors.GetColumn(0).Trans().Info(12,3);
                }
                throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] Error solving system of equations using mumps.");
            }

            //std::cout << " rhsVector" << std::endl;
            //rhsVector.Trans().Info(10,3);
            //std::cout << " delta_disp" << std::endl;
            //deltaDisplacementsActiveDOFs.Trans().Info(10,3);

            //perform a linesearch
            alpha = 1.;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;

                //std::cout << " displacementsActiveDOFs" << std::endl;
                //displacementsActiveDOFs.Trans().Info(10,3);
                this->NodeMergeActiveDofValues(displacementsActiveDOFs);
                this->ElementTotalUpdateTmpStaticData();

/*                {
                    double energy;
                    energy = this->ElementTotalGetTotalEnergy();
                    energy += this->ConstraintTotalGetTotalEnergy();
                    std::cout << "energy " << energy << std::endl;
                }
*/
                // calculate residual
                this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                //std::cout << "intForceVector "  << std::endl;
                //intForceVector.Trans().Info(10,3);

                rhsVector = extForceVector - intForceVector;
                normResidual = rhsVector.Norm();
                maxResidual = rhsVector.Abs().Max();

                //std::cout << "total energy of system " << ElementTotalGetTotalEnergy() << std::endl;

//double energyElement(fineScaleStructure->ElementTotalGetTotalEnergy());
//double energyConstraint(fineScaleStructure->ConstraintTotalGetTotalEnergy());
//std::cout << "alpha " << alpha << " normResidual " << normResidual << " energy " << energyElement+energyConstraint <<"(" << energyElement << "," << energyConstraint << ")" << std::endl;
std::cout << "alpha " << alpha << " normResidual " << normResidual << " normInit " << normRHS << std::endl;

                alpha*=0.5;
            }
            while(alpha>1e-5 && normResidual>normRHS*(1-0.5*alpha) && normResidual>rToleranceResidualForce && maxResidual>rToleranceResidualForce);

            this->PostProcessDataAfterLineSearch(loadStep, numNewtonIterations, 2.*alpha, curLoadFactor);

            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>rToleranceResidualForce && maxResidual>rToleranceResidualForce)
            {
                convergenceStatus=2;
                {
                    if (mNumActiveDofs<1000)
                    {
                        NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSRVector2);
                        std::cout << "stiffness full" << std::endl;
                        stiffnessMatrixFull.Info(12,3);
                        NuTo::FullMatrix<double> eigenValues;
                        stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
                        std::cout << "eigenvalues" << std::endl;
                        eigenValues.Trans().Info(12,3);
                        NuTo::FullMatrix<double> eigenVectors;
                        stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
                        std::cout << "eigenvector 1" << std::endl;
                        eigenVectors.GetColumn(0).Trans().Info(12,3);
                    }
                }
                break;
            }


            //std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;

            //check convergence
            if (normResidual<rToleranceResidualForce || maxResidual<rToleranceResidualForce)
            {
                this->PostProcessDataAfterConvergence(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor);
                convergenceStatus=1;
                //CheckStiffness();
                //NodeInfo(12);
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            normRHS = rhsVector.Norm();
            //build new stiffness matrix
            this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            //std::cout << dispForceVector.Norm() << std::endl;
//check stiffness
//CheckStiffness();
            //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
        }

        if (deltaLoadFactor<rMinDeltaLoadFactor)
            throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence, delta strain factor smaller than given tolerance.");

        if (convergenceStatus==1)
        {
            this->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
            if (curLoadFactor>=1-rMinDeltaLoadFactor)
            {
                convergenceStatusLoadSteps=true;
            }
            else
            {
                //this routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
                //for an update on the coarse scale
                //as a consequence, in an iterative solution with updates in between, the initial state has to be restored after leaving the routine
                if (rSaveStructureBeforeUpdate==true && rIsSaved==false)
                {
                    //store the structure only once in order to be able to restore the situation before entering the routine
                	this->SaveStructure(rSaveStringStream);
                	rIsSaved = true;
                }

                // the update is only required to allow for a stepwise solution procedure in the fine scale model
                // a final update is only required for an update on the macroscale, otherwise,the original state has
                // to be reconstructed.
                this->ElementTotalUpdateStaticData();

                //eventually increase load step
                if (rAutomaticLoadstepControl)
                {
                    if (numNewtonIterations<rMinNumNewtonIterations)
                    {
                        deltaLoadFactor*=rIncreaseFactor;
                    }
                    if (deltaLoadFactor>rMaxDeltaLoadFactor)
                        deltaLoadFactor = rMaxDeltaLoadFactor;
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
        }
        else
        {
            assert(convergenceStatus==2);
            if (rAutomaticLoadstepControl==false)
                throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence with the prescibed number of Newton iterations.");

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
            this->ElementTotalUpdateTmpStaticData();

            //decrease load step
            deltaLoadFactor*=rDecreaseFactor;
            curLoadFactor+=deltaLoadFactor;

            //check for minimum delta (this mostly indicates an error in the software
            if (deltaLoadFactor<rMinDeltaLoadFactor)
            {
                throw NuTo::MechanicsException("[NuTo::StructureBase::NewtonRaphson]: No convergence, delta strain factor smaller than minimum.");
            }
        }

        if (!convergenceStatusLoadSteps)
        {
            //update new displacement of RHS

            // build global dof numbering
            this->SetLoadFactor(curLoadFactor);
            this->NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector, but still with previous displacement state
            this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
//CheckStiffness();
            //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            this->ElementTotalUpdateTmpStaticData();

            // calculate initial residual for next load step
            this->BuildGlobalGradientInternalPotentialVector(intForceVector);

            //update rhs vector for next Newton iteration
            rhsVector = extForceVector - intForceVector;
            normRHS = rhsVector.Norm();
            //attention this is only different for the first iteration step (load application)
            //since the internal force due to the applied constraints is not considered for the first iteration
            //in order to balance it (no localization in the boundary region)
            //for the linesearch this internal force has to be considered in order to obtain for a linesearch
            //factor of zero the normRHS
            rhsVector = dispForceVector + extForceVector;
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
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::NewtonRaphson] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
