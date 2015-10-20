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

#include <boost/assign/ptr_map_inserter.hpp>

# ifdef _OPENMP
#include <omp.h>
# endif

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputDummy.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputFullVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"

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
// serializes the class
template void NuTo::Structure::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Structure::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of structure" << std::endl;
#endif
    std::vector<ElementDataBase*> elementDataVector;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase)
    & boost::serialization::make_nvp ("elementMap", mElementMap)
    & boost::serialization::make_nvp ("nodeMap", mNodeMap)
    & BOOST_SERIALIZATION_NVP(elementDataVector);
    std::vector<ElementDataBase*>::iterator itElementData = elementDataVector.begin();
    for (boost::ptr_map<int,ElementBase>::iterator itElement=mElementMap.begin(); itElement!=mElementMap.end(); itElement++,itElementData++)
    {
        itElement->second->SetDataPtr(*itElementData);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structure" << std::endl;
#endif
}

template void NuTo::Structure::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Structure::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Structure::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of structure" << std::endl;
#endif
    std::vector<ElementDataBase*> elementDataVector(mElementMap.size());
    std::vector<ElementDataBase*>::iterator itElementData = elementDataVector.begin();
    for (boost::ptr_map<int,ElementBase>::const_iterator itElement=mElementMap.begin(); itElement!=mElementMap.end(); itElement++,itElementData++)
    {
        *itElementData = itElement->second->GetDataPtr();
    }
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(StructureBase)
    & boost::serialization::make_nvp ("elementMap", mElementMap)
    & boost::serialization::make_nvp ("nodeMap", mNodeMap)
    & BOOST_SERIALIZATION_NVP(elementDataVector);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structure" << std::endl;
#endif
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Save (const std::string &filename, std::string rType )const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
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
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            std::string tmpString(this->GetTypeId());
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
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
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
    std::cout<<"[NuTo::Structure::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Structure::Restore (const std::string &filename, std::string rType )
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
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
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
    std::cout<<"[NuTo::Structure::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

#endif // ENABLE_SERIALIZATION

// based on the global dofs build submatrices of the global coefficent matrix0
NuTo::Error::eError NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType rType, SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
//    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
//    assert(rMatrixJK.GetNumEntries() == 0);

// define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_0_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::DAMPING)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_1_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::MASS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_2_TIME_DERIVATIVE);
    }
    else
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] matrix type not implemented (either stiffness, damping or mass.");

    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_COLUMN_DOF);

    NuTo::Element::eOutput matrixHessianOrderType;
    switch (rType)
    {
    case NuTo::StructureEnum::eMatrixType::STIFFNESS:
        matrixHessianOrderType = Element::HESSIAN_0_TIME_DERIVATIVE;
        break;
    case NuTo::StructureEnum::eMatrixType::DAMPING:
        matrixHessianOrderType = Element::HESSIAN_1_TIME_DERIVATIVE;
        break;
    case NuTo::StructureEnum::eMatrixType::MASS:
        matrixHessianOrderType = Element::HESSIAN_2_TIME_DERIVATIVE;
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] only stiffness, damping and mass matrices handled.");
    }

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] maximum independent set not calculated.");
        if (rMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixJK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
//					bool symmetryFlag = false;
                    Error::eError error;
                    error = elementPtr->Evaluate(elementOutput);

                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

                        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  elementMatrix = (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS) ?
                        //		elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second->GetFullMatrixDouble() :
                        //		elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second->GetFullMatrixDouble();

                        std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                        std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                        //std::cout << "elementMatrix.GetNumRows() " << elementMatrix.GetNumRows() << std::endl;
                        //std::cout << "elementVectorGlobalDofsRow.size() " << elementVectorGlobalDofsRow.size() << std::endl;
                        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                        // write element contribution to global matrix
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                {
                                    if (fabs(elementMatrix(rowCount, colCount))>mToleranceStiffnessEntries)
                                    {
                                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                        if (globalColumnDof < this->mNumActiveDofs)
                                        {
                                            rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                        }
                                        else
                                        {
                                            rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //dont use MIS
        // loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                //std::cout << "thread in structure " << omp_get_thread_num() << "\n";
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                //bool symmetryFlag = false;

                Error::eError error;
                error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)
                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

                    //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  elementMatrix = (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS) ?
                    //		elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second->GetFullMatrixDouble() :
                    //		elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second->GetFullMatrixDouble();

                    std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                    std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                    assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                    // write element contribution to global matrix
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
                            for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                            {
                                if (fabs(elementMatrix(rowCount, colCount))>mToleranceStiffnessEntries)
                                {
                                    int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                    if (globalColumnDof < this->mNumActiveDofs)
                                    {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                        rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                    }
                                    else
                                    {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                        rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        Error::eError error;
        error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)
        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {

            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  elementMatrix = (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS) ?
            //		elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second->GetFullMatrixDouble() :
            //		elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second->GetFullMatrixDouble();

            std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
            std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
            assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());


            /*
             //check stiffness matrix
             FullVector<double, Eigen::Dynamic> check_disp_j1,check_disp_j2,check_disp_k1,check_disp_k2;
             this->NodeExtractDofValues(0,check_disp_j1,check_disp_k1);
             //std::cout << "active dof values " << check_disp_j1 << std::endl;
             FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> elementMatrix_cdf(elementMatrix);
             elementMatrix_cdf.setZero();

             boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutputCDF;
             boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>( elementOutputCDF )( Element::INTERNAL_GRADIENT );
             elementPtr->Evaluate(elementOutputCDF);
             NuTo::FullVector<double,Eigen::Dynamic> elementVector1 = elementOutputCDF.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble();
             //std::cout << "elementVector1 \n" << elementVector1 << std::endl;
             double delta(1e-9);
             for (unsigned int countCol=0; countCol<elementVectorGlobalDofsColumn.size(); countCol++)
             {
             check_disp_j2 = check_disp_j1;
             check_disp_k2 = check_disp_k1;
             if (elementVectorGlobalDofsColumn[countCol]<check_disp_j1.GetNumRows())
             {
             check_disp_j2(elementVectorGlobalDofsColumn[countCol])+=delta;
             }
             else
             {
             check_disp_k2(elementVectorGlobalDofsColumn[countCol]-check_disp_j1.GetNumRows())+=delta;
             }

             this->NodeMergeDofValues(check_disp_j2,check_disp_k2);
             elementPtr->Evaluate(elementOutputCDF);
             NuTo::FullVector<double,Eigen::Dynamic> elementVector2 = elementOutputCDF.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble();
             elementMatrix_cdf.col(countCol) = ((elementVector2-elementVector1))*(1./delta);
             }
             this->NodeMergeDofValues(check_disp_j1,check_disp_k1);
             if ((elementMatrix_cdf-elementMatrix).cwiseAbs().maxCoeff()>1e-0)
             {
             std::cout << "elem stiffness exact\n" <<  elementMatrix << std::endl<< std::endl;;
             std::cout << "dofs \n"; for (unsigned int count=0; count<elementVectorGlobalDofsColumn.size(); count++) std::cout << elementVectorGlobalDofsColumn[count] << "  ";
             std::cout << std::endl << std::endl;
             std::cout << "elem stiffness cdf\n" <<  elementMatrix_cdf << std::endl<< std::endl;;
             std::cout << "delta \n" <<  elementMatrix_cdf-elementMatrix << std::endl<< std::endl;;
             std::cout << std::endl << std::endl;;
             exit(-1);
             }
             else
             {
             std::cout << "error element stiffness is " << (elementMatrix_cdf-elementMatrix).cwiseAbs().maxCoeff() << std::endl;
             }
             */
            // write element contribution to global matrix
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        if (fabs(elementMatrix(rowCount, colCount)) > mToleranceStiffnessEntries)
                        {
                            int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                            if (globalColumnDof < this->mNumActiveDofs)
                            {
                                rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                //std::cout << "add at (" << globalRowDof << "," << globalColumnDof << ") " << elementMatrix(rowCount, colCount) << std::endl;
                            }
                            else
                            {
                                rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        //write contribution of Lagrange Multipliers
        ConstraintsBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK);
    }
    return errorGlobal;
}

// based on the global dofs build submatrices of the global coefficent matrix
NuTo::Error::eError NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureEnum::eMatrixType rType, SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK)
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    //assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    //assert(rMatrixJK.GetNumEntries() == 0);
    assert(rMatrixKJ.IsSymmetric() == false);
    assert(rMatrixKJ.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKJ.GetNumColumns() == this->mNumActiveDofs);
    //assert(rMatrixKJ.GetNumEntries() == 0);
    assert(rMatrixKK.IsSymmetric() == false);
    assert(rMatrixKK.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    //assert(rMatrixKK.GetNumEntries() == 0);

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_0_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::DAMPING)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_1_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::MASS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_2_TIME_DERIVATIVE);
    }
    else
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] matrix type not implemented (either stiffness or mass.");

    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_COLUMN_DOF);

    NuTo::Element::eOutput matrixHessianOrderType;
    switch (rType)
    {
    case NuTo::StructureEnum::eMatrixType::STIFFNESS:
        matrixHessianOrderType = Element::HESSIAN_0_TIME_DERIVATIVE;
        break;
    case NuTo::StructureEnum::eMatrixType::DAMPING:
        matrixHessianOrderType = Element::HESSIAN_1_TIME_DERIVATIVE;
        break;
    case NuTo::StructureEnum::eMatrixType::MASS:
        matrixHessianOrderType = Element::HESSIAN_2_TIME_DERIVATIVE;
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] only stiffness, damping and mass matrices handled.");
    }

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] maximum independent set not calculated.");
        if (rMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixJK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixKJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixKJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixKK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] MatrixKK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            //here != had to replaced by < in order to make it compile under openmp
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
                    //bool symmetryFlag = false;

                    Error::eError error;
                    error = elementPtr->Evaluate(elementOutput);

                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

                        std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                        std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                        // write element contribution to global matrix
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                {
                                    int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                    if (globalColumnDof < this->mNumActiveDofs)
                                    {
                                        rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                    }
                                    else
                                    {
                                        rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                    }

                                }
                            }
                            else
                            {
                                for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                {
                                    int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                    if (globalColumnDof < this->mNumActiveDofs)
                                    {
                                        rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                                    }
                                    else
                                    {
                                        rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //dont use MIS
        // loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                //bool symmetryFlag = false;
                Error::eError error;
                error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)
                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

                    std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                    std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                    assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                    // write element contribution to global matrix
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
                            for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                            {
                                int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                if (globalColumnDof < this->mNumActiveDofs)
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                }
                                else
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                }

                            }
                        }
                        else
                        {
                            for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                            {
                                int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                if (globalColumnDof < this->mNumActiveDofs)
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                                }
                                else
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        //bool symmetryFlag = false;
        Error::eError error;
        error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)
        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> elementMatrix(elementOutput.find(matrixHessianOrderType)->second->GetFullMatrixDouble());

            std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
            std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
            assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

            // write element contribution to global matrix
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }

                    }
                }
                else
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
    }
#endif
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        //write contribution of Lagrange Multipliers
        ConstraintBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK, rMatrixKJ, rMatrixKK);
    }

    return errorGlobal;
}

// based on the global dofs build submatrices of the global stiffness matrix
NuTo::Error::eError NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK)
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == false);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    //assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    //assert(rMatrixJK.GetNumEntries() == 0);
    assert(rMatrixKJ.IsSymmetric() == false);
    assert(rMatrixKJ.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKJ.GetNumColumns() == this->mNumActiveDofs);
    //assert(rMatrixKJ.GetNumEntries() == 0);
    assert(rMatrixKK.IsSymmetric() == false);
    assert(rMatrixKK.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    //assert(rMatrixKK.GetNumEntries() == 0);

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

    boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_COLUMN_DOF);

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral] maximum independent set not calculated.");
        if (rMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral] MatrixJK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixKJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral] MatrixKJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        if (rMatrixKK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticStiffnessSubMatricesGeneral] MatrixKK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            //here != had to replaced by < in order to make it compile under openmp
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
                    //bool symmetryFlag = false;

                    Error::eError error;
                    error = elementPtr->Evaluate(elementOutput);

                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::StructureBase::BuildGlobalElasticStiffnessSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& elementMatrix(elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC)->second->GetFullMatrixDouble());

                        std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                        std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                        // write element contribution to global matrix
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                {
                                    int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                    if (globalColumnDof < this->mNumActiveDofs)
                                    {
                                        rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                    }
                                    else
                                    {
                                        rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                    }

                                }
                            }
                            else
                            {
                                for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                {
                                    int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                    if (globalColumnDof < this->mNumActiveDofs)
                                    {
                                        rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                                    }
                                    else
                                    {
                                        rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //dont use MIS
        // loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                //bool symmetryFlag = false;
                Error::eError error;
                error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)
                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::StructureBase::BuildGlobalElasticStiffnessSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& elementMatrix(elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC)->second->GetFullMatrixDouble());

                    std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
                    std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                    assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                    // write element contribution to global matrix
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
                            for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                            {
                                int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                if (globalColumnDof < this->mNumActiveDofs)
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                }
                                else
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                }

                            }
                        }
                        else
                        {
                            for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                            {
                                int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                if (globalColumnDof < this->mNumActiveDofs)
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                                }
                                else
                                {
#pragma omp critical (StructureBuildGlobalCoefficientSubMatrices0General)
                                    rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        //bool symmetryFlag = false;
        Error::eError error;
        error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)
        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalElasticStiffnessSubMatricesGeneral] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& elementMatrix(elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC)->second->GetFullMatrixDouble());

            std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
            std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
            assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

            // write element contribution to global matrix
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }

                    }
                }
                else
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixKJ.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                        }
                        else
                        {
                            rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
    }
#endif

    //write contribution of Lagrange Multipliers
    ConstraintBuildGlobalCoefficientSubMatrices0General(rMatrixJJ, rMatrixJK, rMatrixKJ, rMatrixKK);

    return errorGlobal;
}

// based on the global dofs build submatrices of the global coefficent matrix0
NuTo::Error::eError NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureEnum::eMatrixType rType, SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == true);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_0_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::MASS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_2_TIME_DERIVATIVE);
    }

    else
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] matrix type not implemented (either stiffness or mass.");

    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_COLUMN_DOF);

    // loop over all elements
#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mMIS.size()==0)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] maximum independent set not calculated.");
    if (rMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
    if (rMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");

    for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
    {
        std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
        //here != had to replaced by < in order to make it compile under openmp
        for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = *elementIter;
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
#endif
        // calculate element contribution
        Error::eError error;
        error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)
        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesSymmetric] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::ElementOutputBase* outputPtrHessian = ((rType == NuTo::StructureEnum::eMatrixType::STIFFNESS) ? elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second : elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second);

            if (outputPtrHessian->GetSymmetry() == false)
            {
                throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] element matrix is not symmetric (general sparse matrix required).");
            }

            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& elementMatrix = outputPtrHessian->GetFullMatrixDouble();

            std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
            std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
            assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

            // write element contribution to global matrix
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            if (globalColumnDof >= globalRowDof)
                            {
                                rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                            }
                        }
                        else
                        {
                            rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
    }
}
}
#else
    }
#endif
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        //write contribution of Lagrange Multipliers
        ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(rMatrixJJ, rMatrixJK);
    }

    return errorGlobal;
}

// based on the global dofs build submatrices of the global coefficent matrix0
NuTo::Error::eError NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureEnum::eMatrixType rType, SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKK)
{
    assert(this->mNodeNumberingRequired == false);
    assert(rMatrixJJ.IsSymmetric() == true);
    assert(rMatrixJJ.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumColumns() == this->mNumActiveDofs);
    assert(rMatrixJJ.GetNumEntries() == 0);
    assert(rMatrixJK.IsSymmetric() == false);
    assert(rMatrixJK.GetNumRows() == this->mNumActiveDofs);
    assert(rMatrixJK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixJK.GetNumEntries() == 0);
    assert(rMatrixKK.IsSymmetric() == true);
    assert(rMatrixKK.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumColumns() == this->mNumDofs - this->mNumActiveDofs);
    assert(rMatrixKK.GetNumEntries() == 0);

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_0_TIME_DERIVATIVE);
    }
    else if (rType == NuTo::StructureEnum::eMatrixType::MASS)
    {
        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutput)(Element::HESSIAN_2_TIME_DERIVATIVE);
    }
    else
        throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesGeneral] matrix type not implemented (either stiffness or mass.");

    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_COLUMN_DOF);

    // loop over all elements
#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mMIS.size()==0)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] maximum independent set not calculated.");
    if (rMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] MatrixJJ does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
    if (rMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] MatrixJK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
    if (rMatrixKK.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
    throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] MatrixKK does not allow for parallel assembly, use SparseMatrixCSRVector2 instead.");
    for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
    {
        std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
        for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = *elementIter;
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
#endif
        // calculate element contribution
        Error::eError error;
        error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)
        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatricesSymmetric] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::ElementOutputBase* outputPtrHessian = ((rType == NuTo::StructureEnum::eMatrixType::STIFFNESS) ? elementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE)->second : elementOutput.find(Element::HESSIAN_2_TIME_DERIVATIVE)->second);

            if (outputPtrHessian->GetSymmetry() == false)
            {
                throw MechanicsException("[NuTo::Structure::BuildGlobalCoefficientSubMatricesSymmetric] element matrix is not symmetric (general sparse matrix required).");
            }

            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& elementMatrix = outputPtrHessian->GetFullMatrixDouble();

            std::vector<int>& elementVectorGlobalDofsRow(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());
            std::vector<int>& elementVectorGlobalDofsColumn(elementOutput.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
            assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

            // write element contribution to global matrix
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            if (globalColumnDof >= globalRowDof)
                            {
                                rMatrixJJ.AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                            }
                        }
                        else
                        {
                            rMatrixJK.AddValue(globalRowDof, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                        }

                    }
                }
                else
                {
                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                    {
                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                        if (globalColumnDof >= this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            if (globalColumnDof >= globalRowDof)
                            {
                                rMatrixKK.AddValue(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, elementMatrix(rowCount, colCount));
                            }
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
    }
}
}
#else
    }
#endif
    if (rType == NuTo::StructureEnum::eMatrixType::STIFFNESS)
    {
        //write contribution of Lagrange Multipliers
        ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(rMatrixJJ, rMatrixJK, rMatrixKK);
    }
    return errorGlobal;
}

//! @brief ... based on the global dofs build sub-vectors of the global lumped mass
//! @param rActiveDofVector ... global lumped mass which corresponds to the active dofs
//! @param rDependentDofVector ... global lumped mass which corresponds to the dependent dofs
NuTo::Error::eError NuTo::Structure::BuildGlobalLumpedHession2(NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofVector, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofVector)
{
    // initialize vectors
    assert(rActiveDofVector.GetNumRows() == this->mNumActiveDofs);
    assert(rActiveDofVector.GetNumColumns() == 1);
    for (int row = 0; row < rActiveDofVector.GetNumRows(); row++)
    {
        rActiveDofVector(row, 0) = 0.0;
    }
    assert(rDependentDofVector.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rDependentDofVector.GetNumColumns() == 1);
    for (int row = 0; row < rDependentDofVector.GetNumRows(); row++)
    {
        rDependentDofVector(row, 0) = 0.0;
    }

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

    boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>(elementOutput)(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);

    // loop over all elements
#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalGradientInternalPotentialSubVectors] maximum independent set not calculated.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
                    Error::eError error = elementPtr->Evaluate(elementOutput);
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->second->GetFullVectorDouble());
                        std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                        assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                        // write element contribution to global vectors
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofs[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                rActiveDofVector(globalRowDof) += elementVector(rowCount);
                            }
                            else
                            {
                                globalRowDof -= this->mNumActiveDofs;
                                assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                                rDependentDofVector(globalRowDof) += elementVector(rowCount);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //dont use MIS
        //loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                Error::eError error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)

                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->second->GetFullVectorDouble());
                    std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                    // write element contribution to global vectors
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofs[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
#pragma omp critical (StructureBuildGlobalLumpedHession2_Active)
                            rActiveDofVector(globalRowDof) += elementVector(rowCount);
                        }
                        else
                        {
                            globalRowDof -= this->mNumActiveDofs;
                            assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
#pragma omp critical (StructureBuildGlobalLumpedHession2_Dependent)
                            rDependentDofVector(globalRowDof) += elementVector(rowCount);
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        Error::eError error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)

        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::FullVector<double, Eigen::Dynamic>& elementVector(elementOutput.find(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->second->GetFullVectorDouble());
            std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

            // write element contribution to global vectors
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofs[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    rActiveDofVector(globalRowDof) += elementVector(rowCount);
                }
                else
                {
                    globalRowDof -= this->mNumActiveDofs;
                    assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                    rDependentDofVector(globalRowDof) += elementVector(rowCount);
                }
            }
        }
    }
#endif

    return errorGlobal;
}

NuTo::Error::eError NuTo::Structure::BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofGradientVector, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofGradientVector, bool rUpdateHistoryVariables)
{
    // initialize vectors
    assert(rActiveDofGradientVector.GetNumRows() == this->mNumActiveDofs);
    assert(rActiveDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rActiveDofGradientVector.GetNumRows(); row++)
    {
        rActiveDofGradientVector(row, 0) = 0.0;
    }
    assert(rDependentDofGradientVector.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rDependentDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rDependentDofGradientVector.GetNumRows(); row++)
    {
        rDependentDofGradientVector(row, 0) = 0.0;
    }

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

    boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>(elementOutput)(Element::INTERNAL_GRADIENT);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);
    if (rUpdateHistoryVariables)
    {
        boost::assign::ptr_map_insert<ElementOutputDummy>(elementOutput)(Element::UPDATE_STATIC_DATA);
    }
    // loop over all elements
#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalGradientInternalPotentialSubVectors] maximum independent set not calculated.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
                    Error::eError error = elementPtr->Evaluate(elementOutput);
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble());
                        std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                        assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                        // write element contribution to global vectors
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofs[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                            }
                            else
                            {
                                globalRowDof -= this->mNumActiveDofs;
                                assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                                rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //dont use MIS
        //loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                Error::eError error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)

                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble());
                    std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                    // write element contribution to global vectors
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofs[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
#pragma omp critical (StructureBuildGlobalGradientInternalPotentialSubVectors)
                            rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                        }
                        else
                        {
                            globalRowDof -= this->mNumActiveDofs;
                            assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
#pragma omp critical (StructureBuildGlobalGradientInternalPotentialSubVectors)
                            rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        Error::eError error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)

        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::StructureBase::BuildGlobalCoefficientSubMatrices0Symmetric] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::FullVector<double, Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble());
            std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

            // write element contribution to global vectors
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofs[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                }
                else
                {
                    globalRowDof -= this->mNumActiveDofs;
                    assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                    rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                }
            }
        }
    }
#endif

    //write contribution of Lagrange Multipliers
    ConstraintBuildGlobalGradientInternalPotentialSubVectors(rActiveDofGradientVector, rDependentDofGradientVector);
    return errorGlobal;
}

NuTo::Error::eError NuTo::Structure::BuildGlobalElasticGradientInternalPotentialSubVectors(NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofGradientVector, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofGradientVector)
{
    // initialize vectors
    assert(rActiveDofGradientVector.GetNumRows() == this->mNumActiveDofs);
    assert(rActiveDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rActiveDofGradientVector.GetNumRows(); row++)
    {
        rActiveDofGradientVector(row, 0) = 0.0;
    }
    assert(rDependentDofGradientVector.GetNumRows() == this->mNumDofs - this->mNumActiveDofs);
    assert(rDependentDofGradientVector.GetNumColumns() == 1);
    for (int row = 0; row < rDependentDofGradientVector.GetNumRows(); row++)
    {
        rDependentDofGradientVector(row, 0) = 0.0;
    }

    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

    boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>(elementOutput)(Element::INTERNAL_GRADIENT_ELASTIC);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutput)(Element::GLOBAL_ROW_DOF);

    // loop over all elements
#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticGradientInternalPotentialSubVectors] maximum independent set not calculated.");
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
            std::vector<ElementBase*>::iterator elementIter;
#pragma omp parallel default(shared) private(elementIter) firstprivate(elementOutput)
            for (elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // calculate element contribution
                    Error::eError error = elementPtr->Evaluate(elementOutput);
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (errorGlobal==Error::SUCCESSFUL)
                        errorGlobal = error;
                        else if (errorGlobal!=error)
                        throw MechanicsException("[NuTo::Structure::BuildGlobalElasticGradientInternalPotentialSubVectors] elements have returned multiple different error codes, can't handle that.");
                    }
                    else
                    {
                        NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT_ELASTIC)->second->GetFullVectorDouble());
                        std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                        assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                        // write element contribution to global vectors
                        for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                        {
                            int globalRowDof = elementVectorGlobalDofs[rowCount];
                            if (globalRowDof < this->mNumActiveDofs)
                            {
                                rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                            }
                            else
                            {
                                globalRowDof -= this->mNumActiveDofs;
                                assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                                rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        //do not use MIS
        //loop over all elements
#pragma omp parallel default(shared) firstprivate(elementOutput)
        for (boost::ptr_map<int,ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
        {
#pragma omp single nowait
            {
                ElementBase* elementPtr = elementIter->second;
                // calculate element contribution
                Error::eError error = elementPtr->Evaluate(elementOutput);

                if (error!=Error::SUCCESSFUL)

                {
                    if (errorGlobal==Error::SUCCESSFUL)
                    errorGlobal = error;
                    else if (errorGlobal!=error)
                    throw MechanicsException("[NuTo::Structure::BuildGlobalElasticGradientInternalPotentialSubVectors] elements have returned multiple different error codes, can't handle that.");
                }
                else
                {
                    NuTo::FullVector<double,Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT_ELASTIC)->second->GetFullVectorDouble());
                    std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::INTERNAL_GRADIENT_ELASTIC)->second->GetVectorInt());

                    assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

                    // write element contribution to global vectors
                    for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
                    {
                        int globalRowDof = elementVectorGlobalDofs[rowCount];
                        if (globalRowDof < this->mNumActiveDofs)
                        {
#pragma omp critical (StructureBuildGlobalGradientInternalPotentialSubVectors)
                            rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                        }
                        else
                        {
                            globalRowDof -= this->mNumActiveDofs;
                            assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
#pragma omp critical (StructureBuildGlobalGradientInternalPotentialSubVectors)
                            rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                        }
                    }
                }
            }
        }
    }
#else
    // loop over all elements
    for (boost::ptr_map<int, ElementBase>::iterator elementIter = this->mElementMap.begin(); elementIter != this->mElementMap.end(); elementIter++)
    {
        ElementBase* elementPtr = elementIter->second;
        // calculate element contribution
        Error::eError error = elementPtr->Evaluate(elementOutput);

        if (error != Error::SUCCESSFUL)

        {
            if (errorGlobal == Error::SUCCESSFUL)
                errorGlobal = error;
            else if (errorGlobal != error)
                throw MechanicsException("[NuTo::Structure::BuildGlobalElasticGradientInternalPotentialSubVectors] elements have returned multiple different error codes, can't handle that.");
        }
        else
        {
            NuTo::FullVector<double, Eigen::Dynamic>& elementVector(elementOutput.find(Element::INTERNAL_GRADIENT_ELASTIC)->second->GetFullVectorDouble());
            std::vector<int>& elementVectorGlobalDofs(elementOutput.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

            assert(static_cast<unsigned int>(elementVector.GetNumRows()) == elementVectorGlobalDofs.size());

            // write element contribution to global vectors
            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofs.size(); rowCount++)
            {
                int globalRowDof = elementVectorGlobalDofs[rowCount];
                if (globalRowDof < this->mNumActiveDofs)
                {
                    rActiveDofGradientVector(globalRowDof) += elementVector(rowCount);
                }
                else
                {
                    globalRowDof -= this->mNumActiveDofs;
                    assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                    rDependentDofGradientVector(globalRowDof) += elementVector(rowCount);
                }
            }
        }
    }
#endif

    //write contribution of Lagrange Multipliers
    ConstraintBuildGlobalGradientInternalPotentialSubVectors(rActiveDofGradientVector, rDependentDofGradientVector);
    return errorGlobal;
}



//! @brief ... evaluates the structur
void NuTo::Structure::Evaluate(std::map<StructureEnum::eOutput, StructureOutputBase *> &rStructureOutput)
{
    if (rStructureOutput.empty())
    {
        return;
    }

    // Get number of active and dependend Dofs

    int activeDofs      = this->mNumActiveDofs;
    int dependendDofs   = this->mNumDofs - this->mNumActiveDofs;

    // check dof numbering

    if (this->mNodeNumberingRequired)
    {
        try
        {
            this->NodeBuildGlobalDofs();
        }
        catch (MechanicsException& e)
        {
            e.AddMessage("[NuTo::Structure::Evaluate] error building global dof numbering.");
            throw e;
        }
    }

    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::Structure::Evaluate] First update of tmp static data required.");
    }

    // get dof values stored at the nodes
    FullVector<double,Eigen::Dynamic> activeDofValues;
    FullVector<double,Eigen::Dynamic> dependentDofValues;
    try
    {
        this->NodeExtractDofValues(0,activeDofValues, dependentDofValues);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::Structure::Evaluate] error extracting dof values from node.");
        throw e;
    }

    //std::map<NuTo::Element::eOutput,StructureOutputBase *> AssignmentMap;
    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutputMap;

    // loop over all outputs
    for(auto iteratorOutput : rStructureOutput)
    {
        // first switch: sets data types and adjusts sizes
        switch(iteratorOutput.first)
        {
            case NuTo::StructureEnum::eOutput::DAMPING:
            case NuTo::StructureEnum::eOutput::MASS:
            case NuTo::StructureEnum::eOutput::STIFFNESS:
            {

                switch(iteratorOutput.second->GetNumSubmatrices())
                {
                    case 1:
                    {

                        // Get reference to matrix data
                        SparseMatrix<double>& OutputMatrix =iteratorOutput.second->GetSparseMatrixDouble();


#ifdef _OPENMP
                        if (OutputMatrix.AllowParallelAssemblyUsingMaximumIndependentSets()==false)
                        {
                            throw MechanicsException("[NuTo::Structure::Evaluate] Parallel assembly not possible! Use SparseMatrixCSRVector2 for all your output matrices.");
                        }
#endif

                        // Check matrix dimensions and resize if neccessary
                        if (OutputMatrix.GetNumColumns()!=this->mNumActiveDofs || OutputMatrix.GetNumRows()!=this->mNumActiveDofs)
                        {
                            OutputMatrix.Resize(this->mNumActiveDofs, this->mNumActiveDofs);
                        }
                        else
                        {
                            OutputMatrix.SetZeroEntries();
                        }
                        break;
                    }
                    case 4:
                    {


                        // Get reference to matrix data
                        SparseMatrix<double>& OutputMatrixJJ =iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::JJ);
                        SparseMatrix<double>& OutputMatrixJK =iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::JK);
                        SparseMatrix<double>& OutputMatrixKJ =iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::KJ);
                        SparseMatrix<double>& OutputMatrixKK =iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::KK);

#ifdef _OPENMP
                        if (     OutputMatrixJJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false ||
                                 OutputMatrixJK.AllowParallelAssemblyUsingMaximumIndependentSets()==false ||
                                 OutputMatrixKJ.AllowParallelAssemblyUsingMaximumIndependentSets()==false ||
                                 OutputMatrixKK.AllowParallelAssemblyUsingMaximumIndependentSets()==false )
                        {
                            throw MechanicsException("[NuTo::Structure::Evaluate] Parallel assembly not possible! Use SparseMatrixCSRVector2 for all your output matrices and submatrices.");
                        }
#endif


                        // Check matrix dimensions and resize if neccessary
                        if (OutputMatrixJJ.GetNumRows()!= activeDofs || OutputMatrixJJ.GetNumColumns()!= activeDofs)
                        {
                            OutputMatrixJJ.Resize( activeDofs,    activeDofs);
                            OutputMatrixJK.Resize( activeDofs,    dependendDofs);
                            OutputMatrixKJ.Resize( dependendDofs, activeDofs);
                            OutputMatrixKK.Resize( dependendDofs, dependendDofs);
                        }
                        else
                        {
                            OutputMatrixJJ.SetZeroEntries();
                            OutputMatrixJK.SetZeroEntries();
                            OutputMatrixKJ.SetZeroEntries();
                            OutputMatrixKK.SetZeroEntries();
                        }
                        break;
                    }

                    default:
                    {
                        throw NuTo::MechanicsException( std::string("[NuTo::Structure::Evaluate] StructureOutput for matrices with ") +
                                                        std::to_string(iteratorOutput.second->GetNumSubmatrices()) +
                                                        std::string(" submatrices not supported."));
                    }
                }   // switch number of submatrices


                switch(iteratorOutput.first)
                {
                    case NuTo::StructureEnum::eOutput::STIFFNESS:
                    {
                        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutputMap)(Element::HESSIAN_0_TIME_DERIVATIVE);
                        break;
                    }
                    case NuTo::StructureEnum::eOutput::DAMPING:
                    {
                        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutputMap)(Element::HESSIAN_1_TIME_DERIVATIVE);
                        break;
                    }
                    case NuTo::StructureEnum::eOutput::MASS:
                    {
                        boost::assign::ptr_map_insert<ElementOutputFullMatrixDouble>(elementOutputMap)(Element::HESSIAN_2_TIME_DERIVATIVE);
                        break;
                    }
                    default:
                    {
                        throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Element output for request structure output not implemented.");
                    }
                }


                break;
            }   // case StructureOutput is a matrix



            case NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT:
            case NuTo::StructureEnum::eOutput::RESIDUAL_NORM_FACTOR:
            {
                switch(iteratorOutput.second->GetNumSubvectors())
                {
                    case 1:
                    {
                        FullVector<double, Eigen::Dynamic>& OutputVector = iteratorOutput.second->GetFullVectorDouble();


                        // Check matrix dimensions and resize if neccessary
                        if (OutputVector.GetNumRows()!=activeDofs)
                        {
                            OutputVector.Resize(activeDofs);
                        }
                        else
                        {
                            OutputVector.setZero();
                        }
                    break;
                    }
                    case 2:
                    {
                        FullVector<double, Eigen::Dynamic>& OutputVectorJ = iteratorOutput.second->GetFullVectorDouble(StructureEnum::eSubVector::J);
                        FullVector<double, Eigen::Dynamic>& OutputVectorK = iteratorOutput.second->GetFullVectorDouble(StructureEnum::eSubVector::K);


                        // Check matrix dimensions and resize if neccessary
                        if (OutputVectorJ.GetNumRows()!= activeDofs)
                        {
                            OutputVectorJ.Resize(activeDofs);
                            OutputVectorK.Resize(dependendDofs);
                        }
                        else
                        {
                            OutputVectorJ.setZero();
                            OutputVectorK.setZero();
                        }
                        break;
                    }

                    default:
                    {
                        throw NuTo::MechanicsException( std::string("[NuTo::Structure::Evaluate] StructureOutput for vectors with ") +
                                                        std::to_string(iteratorOutput.second->GetNumSubvectors()) +
                                                        std::string(" subvectors not supported."));
                    }
                }   // switch number of subvectors

                switch(iteratorOutput.first)
                {
                    case NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT:
                    {
                        boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>(elementOutputMap)(Element::INTERNAL_GRADIENT);
                        break;
                    }
                    case NuTo::StructureEnum::eOutput::RESIDUAL_NORM_FACTOR:
                    {
                        boost::assign::ptr_map_insert<ElementOutputFullVectorDouble>(elementOutputMap)(Element::RESIDUAL_NORM_FACTOR);
                        break;
                    }

                    default:
                    {
                        throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Element output for request structure output not implemented.");
                    }
                }
                break;
            }


            default:
            {
                throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Output request not implemented.");
            }
        }
    }
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutputMap)(Element::GLOBAL_ROW_DOF);
    boost::assign::ptr_map_insert<ElementOutputVectorInt>(elementOutputMap)(Element::GLOBAL_COLUMN_DOF);


    // define variables storing the element contribution
    Error::eError errorGlobal(Error::SUCCESSFUL);
#ifdef _OPENMP
    if (mNumProcessors!=0)
    {
        omp_set_num_threads(mNumProcessors);
    }

    // use independents sets
    if (mUseMIS)
    {
        if (mMIS.size()==0)
        {
            throw MechanicsException("[NuTo::Structure::Evaluate] maximum independent set not calculated.");
        }
        for (unsigned int misCounter=0; misCounter<mMIS.size(); misCounter++)
        {
#pragma omp parallel default(shared) firstprivate(elementOutputMap)
            {
                //here != had to replaced by < in order to make it compile under openmp
                for (auto elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end(); elementIter++)
                {
#pragma omp single nowait
                    {
                        ElementBase* elementPtr = *elementIter;

#else
    for(auto elementIter : this->mElementMap)
    {
                        ElementBase* elementPtr = elementIter->second;
#endif

                        // calculate element contribution
                        //bool symmetryFlag = false;
                        Error::eError error;

                        error = elementPtr->Evaluate(elementOutputMap);

                        if (error!=Error::SUCCESSFUL)
                        {
                            if (errorGlobal==Error::SUCCESSFUL)
                            errorGlobal = error;
                            else if (errorGlobal!=error)
                            throw MechanicsException("[NuTo::Structure::Evaluate] elements have returned multiple different error codes, can't handle that.");
                        }
                        else
                        {

                            for(auto iteratorOutput : rStructureOutput)
                            {
                                // Needed by Vectors and matrices as well
                                std::vector<int>& elementVectorGlobalDofsRow(elementOutputMap.find(Element::GLOBAL_ROW_DOF)->second->GetVectorInt());

                                switch(iteratorOutput.first)
                                {
                                    case NuTo::StructureEnum::eOutput::DAMPING:
                                    case NuTo::StructureEnum::eOutput::MASS:
                                    case NuTo::StructureEnum::eOutput::STIFFNESS:
                                    {
                                        std::vector<int>& elementVectorGlobalDofsColumn(elementOutputMap.find(Element::GLOBAL_COLUMN_DOF)->second->GetVectorInt());

                                        NuTo::Element::eOutput ElementEnum;

                                        switch(iteratorOutput.first)
                                        {
                                            case NuTo::StructureEnum::eOutput::STIFFNESS:
                                            {
                                                ElementEnum = Element::HESSIAN_0_TIME_DERIVATIVE;
                                                break;
                                            }
                                            case NuTo::StructureEnum::eOutput::DAMPING:
                                            {
                                                ElementEnum = Element::HESSIAN_1_TIME_DERIVATIVE;
                                                break;
                                            }
                                            case NuTo::StructureEnum::eOutput::MASS:
                                            {
                                                ElementEnum = Element::HESSIAN_2_TIME_DERIVATIVE;
                                                break;
                                            }
                                            default:
                                            {
                                                throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Element output for request structure output not implemented.");
                                            }
                                        }


                                        NuTo::ElementOutputBase* elementOutput = elementOutputMap.find(ElementEnum)->second;
                                        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& elementMatrix = elementOutput->GetFullMatrixDouble();;

                                        assert(static_cast<unsigned int>(elementMatrix.GetNumRows()) == elementVectorGlobalDofsRow.size());
                                        assert(static_cast<unsigned int>(elementMatrix.GetNumColumns()) == elementVectorGlobalDofsColumn.size());

                                        switch(iteratorOutput.second->GetNumSubmatrices())
                                        {
                                            case 1:
                                            {
                                                for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                                                {
                                                    int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                                                    for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                                    {
                                                        int globalColumnDof = elementVectorGlobalDofsColumn[colCount];

                                                        iteratorOutput.second->GetSparseMatrixDouble().AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                                    }
                                                }
                                                break;
                                            }

                                            case 4:
                                            {
                                                // write element contribution to global matrix
                                                for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                                                {
                                                    int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                                                    if (globalRowDof < activeDofs)
                                                    {
                                                        for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                                        {
                                                            int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                                            if (globalColumnDof < activeDofs)
                                                            {
                                                                iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::JJ).AddValue(globalRowDof, globalColumnDof, elementMatrix(rowCount, colCount));
                                                            }
                                                            else
                                                            {
                                                                iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::JK).AddValue(globalRowDof, globalColumnDof - activeDofs, elementMatrix(rowCount, colCount));
                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        for (unsigned int colCount = 0; colCount < elementVectorGlobalDofsColumn.size(); colCount++)
                                                        {
                                                            int globalColumnDof = elementVectorGlobalDofsColumn[colCount];
                                                            if (globalColumnDof < activeDofs)
                                                            {
                                                                iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::KJ).AddValue(globalRowDof - activeDofs, globalColumnDof, elementMatrix(rowCount, colCount));
                                                            }
                                                            else
                                                            {
                                                                iteratorOutput.second->GetSparseMatrixDouble(StructureEnum::eSubMatrix::KK).AddValue(globalRowDof - activeDofs, globalColumnDof - activeDofs, elementMatrix(rowCount, colCount));
                                                            }
                                                        }
                                                    }
                                                }
                                                break;
                                            }

                                            default:
                                            {
                                                throw NuTo::MechanicsException( std::string("[NuTo::Structure::Evaluate] StructureOutput for matrices with ") +
                                                                                std::to_string(iteratorOutput.second->GetNumSubmatrices()) +
                                                                                std::string(" submatrices not supported."));
                                            }
                                        }

                                        if(!elementOutput->GetConstant() && iteratorOutput.second->GetConstant())
                                        {
                                            iteratorOutput.second->SetConstant(false);
                                        }

                                        break;
                                    }


                                    case NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT:
                                    case NuTo::StructureEnum::eOutput::RESIDUAL_NORM_FACTOR:
                                    {
                                        NuTo::FullVector<double, Eigen::Dynamic>* elementVector;
                                        switch(iteratorOutput.first)
                                        {
                                            case NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT:
                                            {
                                                elementVector = &elementOutputMap.find(Element::INTERNAL_GRADIENT)->second->GetFullVectorDouble();
                                                break;
                                            }
                                            case NuTo::StructureEnum::eOutput::RESIDUAL_NORM_FACTOR:
                                            {
                                                elementVector = &elementOutputMap.find(Element::RESIDUAL_NORM_FACTOR)->second->GetFullVectorDouble();
                                                break;
                                            }

                                            default:
                                            {
                                                throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Element output for request structure output not implemented.");
                                            }
                                        }


                                        // TODO: Find a better solution
                                        if(elementPtr->GetEnumType() != Element::BOUNDARYELEMENT2DADDITIONALNODE)
                                        {
                                            assert(static_cast<unsigned int>(elementVector->GetNumRows()) == elementVectorGlobalDofsRow.size());
                                        }
                                        switch(iteratorOutput.second->GetNumSubvectors())
                                        {
                                        case 1:
                                        {
                                            // write element contribution to global vectors
                                            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                                            {
                                                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                                                if (globalRowDof < activeDofs)
                                                {
                                                    if(iteratorOutput.first == NuTo::StructureEnum::eOutput::RESIDUAL_NORM_FACTOR)
                                                    {

                                                        if(elementPtr->GetEnumType() != Element::BOUNDARYELEMENT2DADDITIONALNODE &&
                                                           iteratorOutput.second->GetFullVectorDouble()(globalRowDof) < elementVector->GetValue(rowCount))
                                                        {
                                                            iteratorOutput.second->GetFullVectorDouble()(globalRowDof) = elementVector->GetValue(rowCount);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        iteratorOutput.second->GetFullVectorDouble()(globalRowDof) += elementVector->GetValue(rowCount);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case 2:
                                        {

                                            // write element contribution to global vectors
                                            for (unsigned int rowCount = 0; rowCount < elementVectorGlobalDofsRow.size(); rowCount++)
                                            {
                                                int globalRowDof = elementVectorGlobalDofsRow[rowCount];
                                                if (globalRowDof < activeDofs)
                                                {
                                                    iteratorOutput.second->GetFullVectorDouble(StructureEnum::eSubVector::J)(globalRowDof) += elementVector->GetValue(rowCount);
                                                }
                                                else
                                                {
                                                    globalRowDof -= activeDofs;
                                                    assert(globalRowDof < dependendDofs);
                                                    iteratorOutput.second->GetFullVectorDouble(StructureEnum::eSubVector::K)(globalRowDof) += elementVector->GetValue(rowCount);
                                                }
                                            }
                                            break;
                                        }
                                        default:
                                        {
                                            throw NuTo::MechanicsException( std::string("[NuTo::Structure::Evaluate] StructureOutput for vectors with ") +
                                                                            std::to_string(iteratorOutput.second->GetNumSubvectors()) +
                                                                            std::string(" subvectors not supported."));
                                        }
                                        }   // switch number of subvectors

                                        break;
                                    }

                                    default:
                                    {
                                        throw NuTo::MechanicsException("[NuTo::Structure::Evaluate] Output request not implemented.");
                                    }
                                }
                            }
                        }
#ifdef _OPENMP
                    }
                }   // end loop over elements
            }   // end parallel region
        }   // end loop over independent sets
    }
    // DONT use independents sets
    else
    {
        throw("[NuTo::Structure::Evaluate] The use of OpenMP without using independent sets is deprecated and will be removed soon.");
    }
#else
    }   // end loop over elements
#endif
}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(int rConstitutiveId)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    boost::ptr_map<int, ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveId);
    if (itConstitutive == mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::Structure::BuildNonlocalData] Constitutive law with the given identifier does not exist.");

    try
    {
        BuildNonlocalData(itConstitutive->second);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::Structure::BuildNonlocalData] Error calculating nonlocal data.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::BuildNonlocalData] Error calculating nonlocal data.");
    }
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::BuildNonlocalData] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Builds the nonlocal data for integral type nonlocal constitutive models
//! @param rConstitutiveId constitutive model for which the data is build
void NuTo::Structure::BuildNonlocalData(const ConstitutiveBase* rConstitutive)
{
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
            if (elementPtr->GetConstitutiveLaw(theIp) == rConstitutive)
            {
                elementPtr->DeleteNonlocalElements();
                indexElement.push_back(elementPtr);
                indexIp.push_back(theIp);
                indexIpVolume.push_back(ipVolume[theIp]);
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
    /*
     if(EXIT_SUCCESS != 0)
     {
     INTERPRET_INTERN error_mess("ELEMENT_BUILD_NL_ELEMENTS_ANN: Error using ANN library.");
     return -1;
     }
     */
}

//! @brief import from gmsh, creates groups according to gmsh's physical entities and creates an interpolation types for each group
//! @param rFileName .. file name
//! @param rElementData .. element data for the elements to be created
//! @param rIPData .. ip data for the integration points to be created
//! @return .. Matrix [NumGroups x 2] with [: x 0] group ids and [ : x 1] corresponding interpolation type ids
NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ImportFromGmsh(const std::string& rFileName, ElementData::eElementDataType rElementData, IpData::eIpDataType rIPData)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> ids;

    try
    {
        ids = ImportFromGmshAux(rFileName, rElementData, rIPData);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
    }
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::ImportFromGmsh] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
    return ids;
}

NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ImportFromGmsh(const std::string& rFileName, const std::string& rElementData, const std::string& rIPData)
{
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> ids;

    try
    {
        ElementData::eElementDataType elementDataType = ElementData::ElementDataTypeToEnum(rElementData);
        IpData::eIpDataType ipDataType = NuTo::IpData::IpDataTypeToEnum(rIPData);

        ids = ImportFromGmshAux(rFileName, elementDataType, ipDataType);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::Structure::ImportFromGmsh] Error importing from Gmsh.");
    }
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::ImportFromGmsh] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
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

NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> NuTo::Structure::ImportFromGmshAux(const std::string& rFileName, ElementData::eElementDataType rElementData, IpData::eIpDataType rIPData)
{
    const unsigned int num_elm_nodes[21] =
    { 0, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13 };

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
        throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error opening input file for read access.");
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
        throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error reading file information.");
    }

    //std::cout << "version " << mayor_version <<"," << minor_version  << " binary " << binary << " double_size " << double_size << std::endl;

    if (mayor_version != 2)
    {
        throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid file format.");
    }

    if (binary == 0) // read ASCII - file
    {
        getline(file, line);
        if (line != "$EndMeshFormat")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndMeshFormat not found.");
        }

        // begin node section
        getline(file, line);
        if (line != "$Nodes")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Nodes not found.");
        }

        // read number of nodes
        getline(file, line);
        itFirst = line.begin();
        match = phrase_parse(itFirst, line.end(), uint_[ref(num_nodes) = qi::_1], ascii::space);

        if (!match || itFirst != line.end())
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of nodes.");
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
                throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading node data.");
            }
        }

        // end node section
        getline(file, line);
        if (line != "$EndNodes")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndNodes not found.");
        }

        // begin element section
        getline(file, line);
        if (line != "$Elements")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Elements not found.");
        }

        // read number of elements
        getline(file, line);
        itFirst = line.begin();
        match = phrase_parse(itFirst, line.end(), uint_[ref(num_elements) = qi::_1], ascii::space);

        if (!match || itFirst != line.end())
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading number of elements.");
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
                throw MechanicsException("[NuTo::Structure::ImportFromGmsh] error reading element data.");
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
                throw MechanicsException("[NuTo::Structure::ImportFromGmsh] invalid number of element data.");
            }
            tmp_elem_data.clear();
        }

        // end element section
        getline(file, line);
        if (line != "$EndElements")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndElements not found.");
        }
    }

    else // binary format
    {
        if (double_size != sizeof(double))
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid size of double.");
        }

        // close file and open as binary
        file.close();
        file.open(rFileName.c_str(), std::ios::in | std::ios::binary);
        if (file.is_open() == false)
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Error opening input file for read access.");
        }

        // read the first two lines again
        getline(file, line);
        getline(file, line);

        // check size of integer
        int one;
        file.read((char *) &one, sizeof(int));
        if (one != 1)
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Invalid binary format.");
        }
        file.seekg(1, std::ios::cur);

        getline(file, line);
        if (line != "$EndMeshFormat")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndMeshFormat not found.");
        }

        // begin node section
        getline(file, line);
        if (line != "$Nodes")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Nodes not found.");
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
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndNodes not found.");
        }

        // begin element section
        getline(file, line);
        if (line != "$Elements")
        {
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $Elements not found.");
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
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] $EndElements not found.");
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
        throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Only implemented for 2D and 3D.");
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
            shapeType = Interpolation::TRUSSXD;
            typeOrder = Interpolation::EQUIDISTANT1;
            break;
        case 2: // 3-node triangle.
            shapeType = Interpolation::TRIANGLE2D;
            typeOrder = Interpolation::EQUIDISTANT1;
            break;

        case 3: // 4-node quadrangle.
            shapeType = Interpolation::QUAD2D;
            typeOrder = Interpolation::EQUIDISTANT1;
            break;

        case 4: // 4-node tetrahedron.
            shapeType = Interpolation::TETRAHEDRON3D;
            typeOrder = Interpolation::EQUIDISTANT1;
            break;

        case 5: // 8-node hexahedron.
            shapeType = Interpolation::BRICK3D;
            typeOrder = Interpolation::EQUIDISTANT1;
            break;

//    	case 6: // 6-node prism.

//    	case 7: // 5-node pyramid.

    	case 8: // 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
            shapeType = Interpolation::TRUSSXD;
            typeOrder = Interpolation::EQUIDISTANT2;
            //ordering is different than in gmsh, fix this first
            {
            NuTo::FullVector<int, Eigen::Dynamic> nodeNumbersCopy(nodeNumbers);
            nodeNumbers(0)  = nodeNumbersCopy(0);
            nodeNumbers(1)  = nodeNumbersCopy(2);
            nodeNumbers(2)  = nodeNumbersCopy(1);
            }
            break;
        case 9: // 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
            shapeType = Interpolation::TRIANGLE2D;
            typeOrder = Interpolation::EQUIDISTANT2;
            break;

        case 10: // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
            shapeType = Interpolation::QUAD2D;
            typeOrder = Interpolation::LOBATTO2;
            break;

        case 11: // 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
            shapeType = Interpolation::TETRAHEDRON3D;
            typeOrder = Interpolation::EQUIDISTANT2;
            break;

    	case 12: // 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
    	{
    		shapeType = Interpolation::BRICK3D;
            typeOrder = Interpolation::LOBATTO2;
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
            shapeType = Interpolation::QUAD2D;
            typeOrder = Interpolation::EQUIDISTANT2;
            break;

        case 17: // 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
            shapeType = Interpolation::BRICK3D;
            typeOrder = Interpolation::EQUIDISTANT2;
            break;

//    	case 18: // 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).

//    	case 19: // 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

//    	case 20: // 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)

        case 21: // 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
            shapeType = Interpolation::TRIANGLE2D;
            typeOrder = Interpolation::EQUIDISTANT3;
            break;

//    	case 22: // 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)

        case 23: // 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
            shapeType = Interpolation::TRIANGLE2D;
            typeOrder = Interpolation::EQUIDISTANT4;
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
            throw MechanicsException("[NuTo::Structure::ImportFromGmsh] Element type not implemented in the import routine.");
        }

        // get gmsh group id and create a corresponding nuto group if needed
        int groupId = elements[elementCount].tags[0]; // NuTo groupId == gmsh groupId. // This might cause errors if groups exist before the gmsh import.
        boost::ptr_map<int, GroupBase>::iterator itGroup(mGroupMap.find(groupId));
        if (itGroup == mGroupMap.end())
            GroupCreate(groupId, NuTo::Groups::Elements);

        int interpolationTypeId;

        // find the group in the id map, add it if needed
        auto itGroupInterpolation = groupInterpolationIds.find(groupId);
        if (itGroupInterpolation == groupInterpolationIds.end())
        {
            // create a new interpolation type
            interpolationTypeId = InterpolationTypeCreate(Interpolation::ShapeTypeToString(shapeType));
            InterpolationTypeAdd(interpolationTypeId, Node::COORDINATES, typeOrder);

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
                Interpolation::eTypeOrder groupTypeOrder = interpolationType.Get(Node::COORDINATES).GetTypeOrder();
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
                InterpolationTypeAdd(interpolationTypeId, Node::COORDINATES, typeOrder);

                // add it to the map
                groupInterpolationIds[groupId].insert(interpolationTypeId);
            }



        } // now, a valid group id and a valid interpolation type exists

        int elementId = ElementCreate(interpolationTypeId, nodeNumbers, rElementData, rIPData);
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
#ifdef SHOW_TIME
    std::clock_t start, end;
    start = clock();
#endif
    try
    {
        std::map<NodeBase*, NodeBase*> old2NewNodePointer;
        std::map<ElementBase*, ElementBase*> old2NewElementPointer;
        CopyAndTranslate(rOffset, old2NewNodePointer, old2NewElementPointer);
    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::Structure::CopyAndTranslate] Error translating and copying structure.");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException("[NuTo::Structure::CopyAndTranslate] Error translating and copying structure.");
    }
#ifdef SHOW_TIME
    end = clock();
    if (mShowTime)
        std::cout << "[NuTo::Structure::CopyAndTranslate] " << difftime(end, start) / CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief copy and move the structure
//! most of the data is kept, but e.g. nonlocal data and
//! @param rOffset offset (dimension x 1 has to be identical with structure dimension)
//! @param rOld2NewNodePointer ptrMap showing the new and old node pointers
//! @param rOld2NewElementPointer ptrMap showing the new and old element pointers
void NuTo::Structure::CopyAndTranslate(NuTo::FullVector<double, Eigen::Dynamic>& rOffset, std::map<NodeBase*, NodeBase*>& rOld2NewNodePointer, std::map<ElementBase*, ElementBase*>& rOld2NewElementPointer)
{
    if (rOffset.GetNumRows() != mDimension)
        throw MechanicsException("[NuTo::Structure::CopyAndTranslate] offset has to have the same dimension as the structure.");
    if (rOffset.GetNumColumns() != 1)
        throw MechanicsException("[NuTo::Structure::CopyAndTranslate] offset has to have a single column.");

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

        int numCoordinates = nodeVector[countNode]->GetNumCoordinates();

        switch (numCoordinates)
        {
        case 0:
            break;
        case 1:
            newNode->SetCoordinates1D(nodeVector[countNode]->GetCoordinates1D() + rOffset);
            break;
        case 2:
            newNode->SetCoordinates2D(nodeVector[countNode]->GetCoordinates2D() + rOffset);
            break;
        case 3:
            newNode->SetCoordinates3D(nodeVector[countNode]->GetCoordinates3D() + rOffset);
            break;
        default:
            throw MechanicsException("[NuTo::Structure::CopyAndTranslate] number of nodes not supported.");
        }
    }
    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired = true;

    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    std::set<ConstitutiveBase*> constitutiveWithNonlocalData;
    for (unsigned int countElement = 0; countElement < elements.size(); countElement++)
    {
        ElementBase* oldElementPtr = elements[countElement];
        ElementData::eElementDataType elementDataType = oldElementPtr->GetElementDataType();
        IpData::eIpDataType ipDataType = oldElementPtr->GetIpDataType(0);
        int numNodes = oldElementPtr->GetNumNodes();
        std::vector<NodeBase*> nodeVector(numNodes);
        for (int countNode = 0; countNode < numNodes; countNode++)
        {
            nodeVector[countNode] = rOld2NewNodePointer[oldElementPtr->GetNode(countNode)];
        }
        // find interpolation type
        int interpolationTypeId = 0;
        const InterpolationType* interpolationTypeOld = oldElementPtr->GetInterpolationType();
        for (auto it = mInterpolationTypeMap.begin(); it != mInterpolationTypeMap.end(); it++)
            if ((it->second) == interpolationTypeOld)
            {
                interpolationTypeId = it->first;
                break;
            }

        int newElementId = ElementCreate(interpolationTypeId, nodeVector, elementDataType, ipDataType);
        ElementBase* newElementPtr = ElementGetElementPtr(newElementId);
        rOld2NewElementPointer[oldElementPtr] = newElementPtr;

        //set integration type
        const IntegrationTypeBase* integrationType = oldElementPtr->GetIntegrationType();
        newElementPtr->SetIntegrationType(integrationType, ipDataType);

        //set section
        const SectionBase* section = oldElementPtr->GetSection();
        newElementPtr->SetSection(section);

        //set constitutive model
        ConstitutiveBase* constitutive = oldElementPtr->GetConstitutiveLaw(0);
        newElementPtr->SetConstitutiveLaw(constitutive);

        if (oldElementPtr->GetNumNonlocalElements() != 0)
            constitutiveWithNonlocalData.insert(constitutive);

        //set static data
        for (int countIp = 0; countIp < integrationType->GetNumIntegrationPoints(); countIp++)
        {
            ConstitutiveStaticDataBase* clonedStaticData = oldElementPtr->GetStaticData(countIp)->Clone();
            newElementPtr->SetStaticData(countIp, clonedStaticData);
            //newElementPtr->SetStaticData(countIp,(oldElementPtr->GetStaticData(countIp))->Clone());
        }
    }

#ifdef _OPENMP
    //there seems to be a problem with the nearest neighbor search library
#pragma omp critical
#endif
    {
        //rebuild nonlocal data
        for (auto it = constitutiveWithNonlocalData.begin(); it != constitutiveWithNonlocalData.end(); it++)
            BuildNonlocalData(*it);
    }
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Structure)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
