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
#include <string>

#include "nuto/mechanics/structures/StructureBase.h"
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
#include "nuto/visualize/VisualizeComponentDamage.h"
#include "nuto/visualize/VisualizeComponentDisplacement.h"
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
       & BOOST_SERIALIZATION_NVP(mNumDofs)
       & BOOST_SERIALIZATION_NVP(mNumActiveDofs)
       & BOOST_SERIALIZATION_NVP(mNodeNumberingRequired)
       & BOOST_SERIALIZATION_NVP(mConstraintMatrix)
       & BOOST_SERIALIZATION_NVP(mConstraintRHS)
       & BOOST_SERIALIZATION_NVP(mHaveTmpStaticData)
       & BOOST_SERIALIZATION_NVP(mUpdateTmpStaticDataRequired)
       & BOOST_SERIALIZATION_NVP(mToleranceStiffnessEntries)
       & BOOST_SERIALIZATION_NVP(mVisualizeComponents);
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
    Group<ElementBase>::iterator ElementIter = rElementGroup->begin();
    while (ElementIter != rElementGroup->end())
    {
        rElements.push_back(*ElementIter);
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


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
