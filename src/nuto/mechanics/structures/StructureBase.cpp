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
#endif // ENABLE_SERIALIZATION

#include <algorithm>
#include <sstream>
#include <string>

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
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
    std::cout << "start serialization of structure base" << std::endl;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject);
    ar & BOOST_SERIALIZATION_NVP(mDimension);
    //ar & boost::serialization::make_nvp("constitutiveLawMap", mConstitutiveLawMap);
//       & boost::serialization::make_nvp("constraintMap", mConstraintMap)
//       & boost::serialization::make_nvp("loadMap", mLoadMap)
//       & boost::serialization::make_nvp("groupMap", mGroupMap)
//       & boost::serialization::make_nvp("integrationTypeMap", mIntegrationTypeMap)
//       & boost::serialization::make_nvp("sectionMap", mSectionMap);
//         mVisualizeComponents ifdef VISUALIZATION
    std::cout << "finish serialization of structure base" << std::endl;
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
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentDisplacement());
}

//! @brief ... Add engineering strains to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringStrain()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringStrain());
}

//! @brief ... Add engineering plastic strain stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringPlasticStrain()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringPlasticStrain());
}

//! @brief ... Add engineering stress to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentEngineeringStress()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentEngineeringStress());
}

//! @brief ... Add section to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentSection()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentSection());
}
//! @brief ... Add constitutive id to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentConstitutive()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentConstitutive());
}

//! @brief ... Add nonlocal weights to the internal list, which is finally exported via the ExportVtkDataFile command
//! @param rElementId ... Element id
//! @param rIp ... local ip number
void NuTo::StructureBase::AddVisualizationComponentNonlocalWeights(int rElementId, int rIp)
{
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
}

//! @brief ... Add the damage variable to the internal list, which is finally exported via the ExportVtkDataFile command
void NuTo::StructureBase::AddVisualizationComponentDamage()
{
	mVisualizeComponents.push_back(new NuTo::VisualizeComponentDamage());
}


void NuTo::StructureBase::ClearVisualizationComponents()
{
	mVisualizeComponents.clear();
}


void NuTo::StructureBase::ExportVtkDataFile(const std::string& rFileName) const
{
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsTotal(ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName);
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rFileName) const
{
    // find group by name
    const Group<ElementBase>* ElementGroup = dynamic_cast<const Group<ElementBase>*>( this->GroupGetGroupPtr(rGroupIdent));
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsByGroup(ElementGroup,ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName);
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

// build global coefficient matrix0
void NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRGeneral<double>& rMatrix, FullMatrix<double>& rVector)
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
            e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error building global dof numbering.");
            throw e;
        }
    }

    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] First update of tmp static data required.");
    }

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
    rMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
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
        SparseMatrixCSRGeneral<double> transConstraintMatrix = this->mConstraintMatrix.transpose();
        rMatrix -= transConstraintMatrix * coefficientMatrixKJ + coefficientMatrixJK * this->mConstraintMatrix;
        rMatrix += transConstraintMatrix * coefficientMatrixKK * this->mConstraintMatrix;

        // build equivalent load vector
        rVector = (transConstraintMatrix * coefficientMatrixKK - coefficientMatrixJK) * (this->mConstraintRHS - dependentDofValues - this->mConstraintMatrix * activeDofValues);
    }
}

// build global coefficient matrix0
void NuTo::StructureBase::BuildGlobalCoefficientMatrix0(SparseMatrixCSRSymmetric<double>& rMatrix, FullMatrix<double>& rVector)
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
            e.AddMessage("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] error building global dof numbering.");
            throw e;
        }
    }
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
    	throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientMatrix0] First update of tmp static data required.");
    }


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
    rMatrix.Resize(this->mNumActiveDofs);
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
}

// build global external load vector
void NuTo::StructureBase::BuildGlobalExternalLoadVector(NuTo::FullMatrix<double>& rVector)
{
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
}

// build global gradient of the internal potential (e.g. the internal forces)
void NuTo::StructureBase::BuildGlobalGradientInternalPotentialVector(NuTo::FullMatrix<double>& rVector)
{
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
}
