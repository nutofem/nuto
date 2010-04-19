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

#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
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

    mMappingIntEnum2String.resize(NuTo::IntegrationTypeBase::NumIntegrationTypes);
    mMappingIntEnum2String[NuTo::IntegrationTypeBase::IntegrationType1D2NGauss1Ip]=
        NuTo::IntegrationType1D2NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationTypeBase::IntegrationType1D2NGauss2Ip]=
        NuTo::IntegrationType1D2NGauss2Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationTypeBase::IntegrationType3D8NGauss1Ip]=
        NuTo::IntegrationType3D8NGauss1Ip::GetStrIdentifierStatic();
    mMappingIntEnum2String[NuTo::IntegrationTypeBase::IntegrationType3D8NGauss2x2x2Ip]=
        NuTo::IntegrationType3D8NGauss2x2x2Ip::GetStrIdentifierStatic();
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
void NuTo::StructureBase::ExportVtkDataFile(const std::string& rFileName, const std::string& rWhat) const
{
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsTotal(ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName, rWhat);
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(const std::string& rGroupIdent, const std::string& rFileName, const std::string& rWhat) const
{
    // find group by name
    const Group<ElementBase>* ElementGroup = dynamic_cast<const Group<ElementBase>*>( this->GroupGetGroupPtr(rGroupIdent));
    std::vector<const ElementBase*> ElementVec;
    this->GetElementsByGroup(ElementGroup,ElementVec);
    this->ExportVtkDataFile(ElementVec, rFileName, rWhat);
}

void NuTo::StructureBase::ExportVtkDataFile(const std::vector<const ElementBase*>& rElements, const std::string& rFileName, const std::string& rWhat) const
{
    std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat> enumWhatMap;
    // split string
    std::vector<std::string> tokens;
    std::istringstream iss(rWhat);
    std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter<std::vector<std::string> >(tokens));

    std::vector<std::string>::iterator tokenIter = tokens.begin();
    while (tokenIter != tokens.end())
    {
        // convert to uppercase
        std::transform(tokenIter->begin(), tokenIter->end(), tokenIter->begin(), (int(*)(int)) toupper);
        // find enum type
        if (*tokenIter == "DISPLACEMENTS")
        {
            enumWhatMap[*tokenIter] = NuTo::VisualizeBase::DISPLACEMENTS;
        }
        else if (*tokenIter == "ENGINEERING_STRESS")
        {
            enumWhatMap[*tokenIter] = NuTo::VisualizeBase::ENGINEERING_STRESS;
        }
        else if (*tokenIter == "ENGINEERING_STRAIN")
        {
            enumWhatMap[*tokenIter] = NuTo::VisualizeBase::ENGINEERING_STRAIN;
        }
        else
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] invalid data description.");
        }
        tokenIter++;
    }

    // call export routine
    this->ExportVtkDataFile(rElements, rFileName, enumWhatMap);
}

void NuTo::StructureBase::ExportVtkDataFile(const std::vector<const ElementBase*>& rElements, const std::string& rFileName, const std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat>& rWhat) const
{
    VisualizeUnstructuredGrid Visualize;
    std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat>::const_iterator itWhat = rWhat.begin();
    while (itWhat != rWhat.end())
    {
        switch (itWhat->second)
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
            Visualize.DefinePointDataVector(itWhat->first);
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            Visualize.DefineCellDataTensor(itWhat->first);
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::ExportVtkDataFile] invalid data description.");
        }
        itWhat++;
    }
    for (unsigned int ElementCount = 0; ElementCount < rElements.size(); ElementCount++)
    {
        rElements[ElementCount]->Visualize(Visualize, rWhat);
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
