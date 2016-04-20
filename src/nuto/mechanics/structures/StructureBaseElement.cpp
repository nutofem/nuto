// $Id$

#include <eigen3/Eigen/Eigenvalues>

#include <assert.h>
#include <boost/tokenizer.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

#include "nuto/base/Timer.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputFullVectorDouble.h"

#include "nuto/mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputBlockVectorInt.h"

#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/ElementOutputDummy.h"

#include "nuto/mechanics/elements/IpDataStaticDataBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"



NuTo::BlockFullVector<double> NuTo::StructureBase::ElementBuildInternalGradient(ElementBase* rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());


    try
    {
        rElement->Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] MechanicsError in element " + std::to_string(ElementGetId(rElement)));
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Non-mechanics error in element " + std::to_string(ElementGetId(rElement)));
    }

    return elementOutputMap.at(Element::INTERNAL_GRADIENT)->GetBlockFullVectorDouble();
}

NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian(Element::eOutput rHessianType, ElementBase* rElement)
{
    std::set<Element::eOutput> supportedTypes({Element::HESSIAN_0_TIME_DERIVATIVE, Element::HESSIAN_1_TIME_DERIVATIVE, Element::HESSIAN_2_TIME_DERIVATIVE});
    if (supportedTypes.find(rHessianType) == supportedTypes.end())
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] requested matrix type is not supported or not implemented yet.");


    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[rHessianType] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());

    try
    {
        rElement->Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] MechanicsError in element " + std::to_string(ElementGetId(rElement)));
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Non-mechanics error in element " + std::to_string(ElementGetId(rElement)));
    }

    return elementOutputMap.at(rHessianType)->GetBlockFullMatrixDouble();
}

NuTo::BlockFullVector<int> NuTo::StructureBase::ElementBuildGlobalDofsRow(ElementBase* rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    try
    {
        rElement->Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "MechanicsError in element " + std::to_string(ElementGetId(rElement)));
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(rElement)));
    }

    return elementOutputMap.at(Element::GLOBAL_ROW_DOF)->GetBlockFullVectorInt();
}

NuTo::BlockFullVector<int> NuTo::StructureBase::ElementBuildGlobalDofsColumn(ElementBase* rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::GLOBAL_COLUMN_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    try
    {
        rElement->Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "MechanicsError in element " + std::to_string(ElementGetId(rElement)));
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(rElement)));
    }

    return elementOutputMap.at(Element::GLOBAL_COLUMN_DOF)->GetBlockFullVectorInt();
}


NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian0_CDF(ElementBase* rElement, double rDelta)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    auto internalGradient0 = ElementBuildInternalGradient(rElement);
    auto globalColumnDofs  = ElementBuildGlobalDofsColumn(rElement);
    auto dofs = rElement->GetInterpolationType()->GetActiveDofs();



    NuTo::BlockFullMatrix<double> hessian0_CDF(GetDofStatus());

    auto dofValues = this->NodeExtractDofValues(0);


    for (auto dofRow : dofs)
    {
        for (auto dofCol : dofs)
        {
            auto& dofActDofValues       = dofValues.J[dofCol];
            auto& dofDepDofValues       = dofValues.K[dofCol];
            auto& hessian0_CDF_dof      = hessian0_CDF(dofRow, dofCol);
            auto& globalColumnDofsCol   = globalColumnDofs[dofCol];

            int numCols = globalColumnDofsCol.GetNumRows();
            int numRows = internalGradient0[dofRow].GetNumRows();

            hessian0_CDF_dof.resize(numRows,numCols);
            hessian0_CDF_dof.setZero();
            for (int iCol = 0; iCol < numCols; ++iCol)
            {
                // Apply rDelta to the corresponding Dof
                int index = globalColumnDofsCol(iCol);
                int numActiveDofs = GetNumActiveDofs(dofCol);

                if (index < numActiveDofs)
                    dofActDofValues(index) += rDelta;
                else
                    dofDepDofValues(index-numActiveDofs) += rDelta;

                NodeMergeDofValues(0,dofValues);
//                 this->ElementTotalUpdateTmpStaticData();

                // Calculate CDF
                hessian0_CDF_dof.SetColumn(iCol, (ElementBuildInternalGradient(rElement)[dofRow]-internalGradient0[dofRow])/rDelta);

                // Restore the orignial dof values
                if (index < numActiveDofs)
                    dofActDofValues(index) -= rDelta;
                else
                    dofDepDofValues(index-numActiveDofs) -= rDelta;
            }   // end of loop over all columns
            this->NodeMergeDofValues(0,dofValues);
        }
    }
    return hessian0_CDF;
}

bool NuTo::StructureBase::ElementCheckHessian0(ElementBase* rElement, double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    bool isHessianCorrect = true;

    auto hessianRef = ElementBuildHessian0_CDF(rElement, rDelta);
    auto hessian    = ElementBuildHessian0(rElement);

    auto differenceAbsolute = hessianRef - hessian;

    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");

    for (auto dofRow : GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : GetDofStatus().GetActiveDofTypes())
        {
            double scaling = hessianRef(dofRow, dofCol).cwiseAbs().maxCoeff();
            FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  differenceRelative = differenceAbsolute(dofRow, dofCol) / scaling;
            int row = 0, col = 0;
            double maxRelativeDifference = differenceRelative.cwiseAbs().maxCoeff(&row, &col);
            if (maxRelativeDifference > rRelativeTolerance)
            {
                isHessianCorrect = false;
                GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 at dof combination [" << Node::DofToString(dofRow) << " - " << Node::DofToString(dofCol) << "]\n";
                GetLogger() << "maxRelativeDifference " << maxRelativeDifference << " at entry (" << row << "," << col << ")\n";
                if (rPrintWrongMatrices)
                {
                    differenceRelative.SetSmallEntriesZero(1.e-10);
                    GetLogger() << "####### relative difference\n" << differenceRelative.format(fmt) << "\n";
                    GetLogger() << "####### hessian0\n" << hessian(dofRow, dofCol).format(fmt) << "\n";
                    GetLogger() << "####### hessian0_CDF\n" << hessianRef(dofRow, dofCol).format(fmt) << "\n";
                }
            }
        }
    }

    return isHessianCorrect;
}

bool NuTo::StructureBase::ElementCheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    bool showTime = GetShowTime();
    SetShowTime(false);

    std::vector<std::pair<int, ElementBase*>> elements;
    GetElementsTotal(elements);

    bool areAllElementsCorrect = true;

    for (auto elementIdPair : elements)
    {
        bool isElementCorrect = ElementCheckHessian0(elementIdPair.second, rDelta, rRelativeTolerance, rPrintWrongMatrices);
        areAllElementsCorrect = areAllElementsCorrect && isElementCorrect;

        if (not isElementCorrect)
        {
            GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 in " << Element::ElementTypeToString(elementIdPair.second->GetEnumType()) << " " << elementIdPair.first << "\n";
            GetLogger() << "################################################################################################\n";
        }

    }
    SetShowTime(showTime);

    return areAllElementsCorrect;
}

//! @brief sets the constitutive law of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    try
    {
    	ElementSetConstitutiveLaw(elementPtr,itConstitutive->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element " + ss.str() + ".");
    }
}

//! @brief sets the constitutive law of a single element
//! @param rElementIdent identifier for the element
//! @param rIp  id of integration point
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId,int rIp, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    try
    {
    	ElementSetConstitutiveLaw(elementPtr,rIp,itConstitutive->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetConstitutiveLaw] Error setting constitutive law  for element " + ss.str() + ".");
    }
}
//! @brief sets the constitutive law of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetConstitutiveLaw(int rGroupIdent, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(itElement->second,itConstitutive->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law  for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting constitutive law for element " + ss.str() + ".");
        }
    }
}

//! @brief sets the constitutive law of a all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Constitutive law with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetConstitutiveLaw(elementVector[countElement],itConstitutive->second);
        	if (elementVector[countElement]->GetNumNonlocalElements()>0)
        	    elementVector[countElement]->DeleteNonlocalElements();
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetConstitutiveLaw] Error setting constitutive law for element "
        			   + ss.str() + ".");
        }
    }
}

//! @brief sets the constitutive law of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive)
{
	rElement->SetConstitutiveLaw(rConstitutive);
}

//! @brief sets the constitutive law of a single ip at an element
//! @param rElement element pointer
//! @param rIp number of integration point
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement,int rIp, ConstitutiveBase* rConstitutive)
{
    //std::cout<< "[NuTo::StructureBase::ElementSetConstitutiveLaw]" << "\n";
	rElement->SetConstitutiveLaw(rIp,rConstitutive);
}


//! @brief sets the section of a single element
//! @param rElementIdent identifier for the element
//! @param rConstitutiveLawIdent identifier for the section
void NuTo::StructureBase::ElementSetSection(int rElementId, int rSectionId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetSection] Section with the given identifier does not exist.");

    try
    {
    	ElementSetSection(elementPtr,itSection->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetSection] Error setting section for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementSetSection] Error setting section for element " + ss.str() + ".");
    }
}

//! @brief sets the section of a group of elements
//! @param rGroupIdent identifier for the group of elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementGroupSetSection(int rGroupIdent, int rSectionId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Section with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	ElementSetSection(itElement->second,itSection->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Error setting section for element " + ss.str() + ".");
        }
    }
}

//! @brief sets the section for all elements
//! @param rConstitutiveLawIdent identifier for the material
void NuTo::StructureBase::ElementTotalSetSection(int rSectionId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,SectionBase>::iterator itSection = mSectionMap.find(rSectionId);
    if (itSection==mSectionMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementTotalSetSection] Section with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        try
        {
        	ElementSetSection(elementVector[countElement],itSection->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalSetSection] Error setting section  for element "
            		+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementTotalSetSection] Error setting section for element "
        			   + ss.str() + ".");
        }
    }
}

//! @brief modifies the material of a single element
//! @param rElement element pointer
//! @param rConstitutive material pointer
void NuTo::StructureBase::ElementSetSection(ElementBase* rElement, SectionBase* rSection)
{
    rElement->SetSection(rSection);
}

//! @brief modifies the interpolation type of a single element
//! @param rElementId ... element number
//! @param rInterpolationTypeId ... interpolation type id
void NuTo::StructureBase::ElementSetInterpolationType(int rElementId, int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementSetInterpolationType] Interpolation type with the given identifier does not exist.");

    try
    {
        ElementSetInterpolationType(elementPtr, itInterpolationType->second);
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementSetInterpolationType] Error setting interpolation type for element "
            + ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
        throw NuTo::MechanicsException
           ("[NuTo::StructureBase::ElementSetInterpolationType] Error setting interpolation type for element " + ss.str() + ".");
    }
}

//! @brief modifies the interpolation type of a group of elements
//! @param rGroupId ... identifier for the group of elements
//! @param rInterpolationTypeId ... interpolation type id
void NuTo::StructureBase::ElementGroupSetInterpolationType(int rGroupId, int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Interpolation type with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
            ElementSetInterpolationType(itElement->second, itInterpolationType->second);
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupSetInterpolationType] Error setting interpolation type for element "
                + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupSetInterpolationType] Error setting interpolation type for element " + ss.str() + ".");
        }
    }
}

//! @brief modifies the interpolation type of a single element
//! @param rElement element pointer
//! @param rInterpolationType interpolation type
void NuTo::StructureBase::ElementSetInterpolationType(ElementBase* rElement, InterpolationType* rInterpolationType)
{
    rElement->SetInterpolationType(rInterpolationType);
}


//! @brief returns the enum of string identifier for an integration type
//! @param rIpDataTypeStr string
//! @return enum
NuTo::IpData::eIpDataType NuTo::StructureBase::ElementGetEnumIntegrationType(const std::string& rIpDataTypeStr)
{
    // get ip data type
    std::string upperCaseIpDataTypeStr;
    std::transform(rIpDataTypeStr.begin(), rIpDataTypeStr.end(), std::back_inserter(upperCaseIpDataTypeStr), (int(*)(int)) toupper);

    NuTo::IpData::eIpDataType ipDataType;
    if (upperCaseIpDataTypeStr=="NOIPDATA")
    {
    	ipDataType = NuTo::IpData::NOIPDATA;
    }
    else if (upperCaseIpDataTypeStr=="STATICDATA")
	{
    	ipDataType = NuTo::IpData::STATICDATA;
	}
    else if (upperCaseIpDataTypeStr=="STATICDATANONLOCAL")
    {
    	ipDataType = NuTo::IpData::STATICDATANONLOCAL;
    }
    else
    {
    	throw MechanicsException("[NuTo::Structure::ElementGetEnumIntegrationType] Ip data type "+upperCaseIpDataTypeStr +" does not exist.");
    }
    return ipDataType;
}


//! @brief calculates static ip data
//! @param rElemIdent  element number
//! @param rType static ip data type
//! @param rIPData matrix with (... x numIP), x varies depending on IPData type
NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> NuTo::StructureBase::ElementGetStaticIPData(int rElementId, IpData::eIpStaticDataType rType)
{
    Timer timer(std::string(__FUNCTION__) + ":" + IpData::IpStaticDataTypeToString(rType), GetShowTime(), GetLogger());

    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] First update of tmp static data required.");

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
        elementOutputMap[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(rType);

        elementPtr->Evaluate(elementOutputMap);

        return elementOutputMap.at(Element::IP_DATA)->GetIpData().GetIpDataMap().at(rType);

    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting " + IpData::IpStaticDataTypeToString(rType) + " for element " + std::to_string(rElementId) + ".");
        throw e;
    } catch (...)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting " + IpData::IpStaticDataTypeToString(rType) + " for element "  + std::to_string(rElementId)  + ".");
    }
}


//! @brief calculates the global integration point coordinates
//! @param rElemIdent  identifier for the element
//! @param rCoordinates integration point coordinates (return value, always 3xnumIp matrix)
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::StructureBase::ElementGetIntegrationPointCoordinates(int rElementId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] First update of tmp static data required.");
    }



    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
		//evaluate the coordinates
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> coordinates (3,elementPtr->GetNumIntegrationPoints());
    	for (int count=0; count<elementPtr->GetNumIntegrationPoints(); count++)
    	{
    	    Eigen::Vector3d coords = elementPtr->GetGlobalIntegrationPointCoordinates(count);
    	    coordinates.SetBlock(0, count, coords);
    	}
    	return coordinates;
    }
    catch(NuTo::MechanicsException &e)
    {
        std::stringstream ss;
        ss << rElementId;
        e.AddMessage("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] Error getting integration point coordinates for element "
        	+ ss.str() + ".");
        throw e;
    }
    catch(...)
    {
        std::stringstream ss;
        ss << rElementId;
    	throw NuTo::MechanicsException
    	   ("[NuTo::StructureBase::ElementGetIntegrationPointCoordinates] Error getting integration point coordinates for element " + ss.str() + ".");
    }
}

//! @brief calculates the maximum damage in all elements
//! @param rElemIdent  identifier for the element
//! @return max damage value
double NuTo::StructureBase::ElementTotalGetMaxDamage()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] First update of tmp static data required.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::DAMAGE);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    double maxDamage = 0;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rIPDamage;

    for (auto element : elementVector)
    {
        try
        {
            element->Evaluate(elementOutputMap);
            rIPDamage = elementOutputMap.at(Element::IP_DATA)->GetIpData().GetIpDataMap()[IpData::DAMAGE];
        } catch (NuTo::MechanicsException &e)
        {
            e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting damage for element " + std::to_string(ElementGetId(element)) + ".");
            throw e;
        } catch (...)
        {
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting damage for element " +std::to_string(ElementGetId(element)) + ".");
        }
        maxDamage = std::max(maxDamage, rIPDamage.Max());

    }

    return maxDamage;
}

double NuTo::StructureBase::ElementTotalGetStaticDataExtrapolationError()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] First update of tmp static data required.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::EXTRAPOLATION_ERROR);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    double max = 0;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ipValues;

    for (auto element : elementVector)
    {
        try
        {
            element->Evaluate(elementOutputMap);
            ipValues = elementOutputMap.at(Element::IP_DATA)->GetIpData().GetIpDataMap()[IpData::EXTRAPOLATION_ERROR];
        } catch (NuTo::MechanicsException &e)
        {
            e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting EXTRAPOLATION_ERROR for element " + std::to_string(ElementGetId(element)) + ".");
            throw e;
        } catch (...)
        {
            throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Error getting EXTRAPOLATION_ERROR for element " +std::to_string(ElementGetId(element)) + ".");
        }
//        if (ipValues.Max() > max)
//        {
//            std::cout << "local max error in element " << ElementGetId(element) << std::endl;
//            ipValues.Info(10,5,true);
//        }
        max = std::max(max, ipValues.Max());
    }

    return max;
}


//! @brief allocates additional static data for an element group
//! @param rElementGroupId ... element group id
//! @param rNumAdditionalStaticData ... number of addidional static data objects
void NuTo::StructureBase::ElementGroupAllocateAdditionalStaticData(int rElementGroupId, int rNumAdditionalStaticData)
{
    if (GroupGetGroupPtr(rElementGroupId)->GetType() != Groups::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Element group required.");

    auto elementIds = GroupGetMemberIds(rElementGroupId);
    for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
    {
        ElementBase* element = ElementGetElementPtr(elementIds[iElement]);
        for (int iIp = 0; iIp < element->GetNumIntegrationPoints(); ++iIp)
        {
            element->GetStaticDataBase(iIp).AllocateAdditionalStaticData(rNumAdditionalStaticData);
        }
    }

}

//! @brief updates the history data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementTotalUpdateStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        try
        {
            Error::eError error = this->ElementTotalUpdateTmpStaticData();
            if (error!=Error::SUCCESSFUL)
            	return error;
        }
        catch (NuTo::Exception& e)
        {
        	e.AddMessage("[NuTo::StructureBase::ElementTotalUpdateStaticData] error building tmp static data.");
        	throw e;
        }
        catch(...)
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateStaticData] error building tmp static data.");
        }
    }

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
	Error::eError errorGlobal (Error::SUCCESSFUL);
    int exception(0);
    std::string exceptionStringTotal;

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();

#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
	{
		try
		{
            Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
            if (error!=Error::SUCCESSFUL)
            {
            	if (errorGlobal==Error::SUCCESSFUL)
            		errorGlobal = error;
            	else if (errorGlobal!=error)
            		throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateStaticData] elements have returned multiple different error codes, can't handle that.");
            }
		}
		catch(NuTo::Exception& e)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal(e.ErrorMessage()
					+"[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
		catch(...)
		{
			std::stringstream ss;
			ss << ElementGetId(elementVector[countElement]);
			std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() + ".\n");
			exception+=1;
			exceptionStringTotal+=exceptionStringLocal;
		}
	}
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
    return errorGlobal;
}

//! @brief updates the history data of a all elements
NuTo::Error::eError NuTo::StructureBase::ElementTotalUpdateTmpStaticData()
{
    Error::eError errorGlobal (Error::SUCCESSFUL);

    //std::cout << "do we really have tmp static data " << mHaveTmpStaticData << std::endl;
    if (mHaveTmpStaticData)
	{
		std::vector<ElementBase*> elementVector;
		GetElementsTotal(elementVector);
	    int exception(0);
	    std::string exceptionStringTotal;

	    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
	    elementOutput[Element::UPDATE_TMP_STATIC_DATA] = std::make_shared<ElementOutputDummy>();


#ifdef _OPENMP
		if (mNumProcessors!=0)
			omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
    #pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
	    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
		{
			try
			{
				Error::eError error = elementVector[countElement]->Evaluate(elementOutput);
				if (error!=Error::SUCCESSFUL)
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] elements have returned multiple different error codes, can't handle that.");
				}
			}
			catch(NuTo::Exception& e)
			{
				std::stringstream ss;
				ss << ElementGetId(elementVector[countElement]);
				std::string exceptionStringLocal(e.ErrorMessage()
						+"[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
			catch(...)
			{
				std::stringstream ss;
				ss << ElementGetId(elementVector[countElement]);
				std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
		}
		if(exception>0)
		{
			throw MechanicsException(exceptionStringTotal);
		}
	}
	//std::cout << "NuTo::StructureBase::ElementTotalUpdateTmpStaticData " << mUpdateTmpStaticDataRequired << "\n";
	mUpdateTmpStaticDataRequired = false;
    return errorGlobal;
}

//! @brief saves static data of all elements
void NuTo::StructureBase::ElementTotalSaveStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        try
        {
            ElementBase* element = elementVector[iElement];
            for (int iIP = 0; iIP < element->GetNumIntegrationPoints(); ++iIP)
            {
                auto& ipDataBase = elementVector[iElement]->GetStaticDataBase(iIP);
                ipDataBase.SaveStaticData();
            }
        } catch (NuTo::Exception& e)
        {
            exception = e;
            exception.AddMessage(__PRETTY_FUNCTION__, "Error saving static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        } catch (...)
        {
            exception.AddMessage(__PRETTY_FUNCTION__, "Non-NuTo Error saving static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        }
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}

//! @brief restores static data of a all elements
void NuTo::StructureBase::ElementTotalRestoreStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        try
        {
            ElementBase* element = elementVector[iElement];
            for (int iIP = 0; iIP < element->GetNumIntegrationPoints(); ++iIP)
            {
                auto& ipDataBase = elementVector[iElement]->GetStaticDataBase(iIP);
                ipDataBase.RestoreStaticData();
            }
        } catch (NuTo::Exception& e)
        {
            exception = e;
            exception.AddMessage(__PRETTY_FUNCTION__, "Error restoring static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        } catch (...)
        {
            exception.AddMessage(__PRETTY_FUNCTION__, "Non-NuTo Error restoring static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        }
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}

//! @brief saves static data of a all elements
void NuTo::StructureBase::ElementTotalExtrapolateStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();

    ConstitutiveCalculateStaticData calculateStaticData(CalculateStaticData::EULER_FORWARD);
    ConstitutiveInputMap input;
    input[Constitutive::Input::CALCULATE_STATIC_DATA] = &calculateStaticData;

#ifdef _OPENMP
    if (mNumProcessors!=0)
    omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        try
        {
            ElementBase* element = elementVector[iElement];
            element->Evaluate(input, elementOutput);

        } catch (NuTo::Exception& e)
        {
            exception = e;
            exception.AddMessage(__PRETTY_FUNCTION__, "Error restoring static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        } catch (...)
        {
            exception.AddMessage(__PRETTY_FUNCTION__, "Non-NuTo Error restoring static data for element " + std::to_string(ElementGetId(elementVector[iElement])));
        }
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}


//! @brief calculates the average stress
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
//! @param rEngineeringStress  average stress (return value)
void NuTo::StructureBase::ElementTotalGetAverageStress(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStress;
    rEngineeringStress.Resize(6,1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
            elementVector[elementCount]->GetIntegratedStress(elementEngineeringStress);
            rEngineeringStress+=elementEngineeringStress;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStress] Error calculating integrated stress  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStress] Error calculating integrated stress  for element " + ss.str() + ".");
        }
    }
    rEngineeringStress*=1./rVolume;

}

//! @brief calculates the average stress
//! @param rGroupId  group number
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
//! @param rEngineeringStress  average stress (return value)
void NuTo::StructureBase::ElementGroupGetAverageStress(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStress] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStress] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStress;
    rEngineeringStress.Resize(6,1);

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	itElement->second->GetIntegratedStress(elementEngineeringStress);
            rEngineeringStress+=elementEngineeringStress;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetAverageStress] Error calculating integrated stress  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetAverageStress] Error calculating integrated stress  for element " + ss.str() + ".");
        }
    }
    rEngineeringStress*=(1./rVolume);

}


//! @brief calculates the average strain
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero strain, but still nonzero area)
//! @param rEngineeringStraiu  average strain (return value)
void NuTo::StructureBase::ElementTotalGetAverageStrain(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStrain;
    rEngineeringStrain.Resize(6,1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
            elementVector[elementCount]->GetIntegratedStrain(elementEngineeringStrain);
            rEngineeringStrain+=elementEngineeringStrain;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=1./rVolume;
}

//! @brief calculates the average strain
//! @param rGroupId  group number
//! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
//! this is a parameter of the model, since holes have to be considered (zero strain, but still nonzero area)
//! @param rEngineeringStrain  average strain (return value)
void NuTo::StructureBase::ElementGroupGetAverageStrain(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStrain] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetAverageStrain] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> elementEngineeringStrain;
    rEngineeringStrain.Resize(6,1);

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	itElement->second->GetIntegratedStrain(elementEngineeringStrain);
            rEngineeringStrain+=elementEngineeringStrain;
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=(1./rVolume);

}


void NuTo::StructureBase::ElementGroupGetMembers(int rGroupId, NuTo::FullVector<int,Eigen::Dynamic>& rMembers)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetMembers] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetMembers] Group is not an element group.");
    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);

    rMembers.Resize(elementGroup->GetNumMembers());
    int countElement(0);
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++,countElement++)
    {
       	rMembers[countElement] = itElement->first;
    }
}


//! @brief calculates the total energy of the system
//! @return total energy
double NuTo::StructureBase::ElementTotalGetInternalEnergy()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    double totalEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    std::vector<const ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalGetInternalEnerg] not yet implemented on ip level.");
//            elementVector[elementCount]->GetIpData(NuTo::IpData::INTERNAL_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                totalEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetInternalEnergy] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetInternalEnergy] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    return totalEnergy;
}

//! @brief calculates the total energy of a group of elements
//! @return total energy
double NuTo::StructureBase::ElementGroupGetTotalEnergy(int rGroupId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Group is not an element group.");
    const Group<ElementBase> *elementGroup = dynamic_cast<const Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    double totalEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementGroupGetTotalEnergy] not yet implemented on ip level.");
        	//itElement->second->GetIpData(NuTo::IpData::INTERNAL_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                totalEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetTotalEnergy] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    return totalEnergy;
}

//! @brief calculates the elastic energy of the system
//! @return elastic energy
double NuTo::StructureBase::ElementTotalGetElasticEnergy()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    double elasticEnergy(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ipEnergy;

    std::vector<const ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount=0; elementCount<elementVector.size();elementCount++)
    {
        try
        {
        	throw MechanicsException("[NuTo::StructureBase::ElementTotalGetElasticEnergy] not yet implemented on ip level.");
        	//elementVector[elementCount]->GetIpData(NuTo::IpData::ELASTIC_ENERGY,ipEnergy);
            for (int theIP=0; theIP<ipEnergy.GetNumColumns(); theIP++)
            {
                elasticEnergy+=ipEnergy(0,theIP)*ipEnergy(1,theIP);
            }
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            e.AddMessage("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementTotalGetAverageStrain] Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    return elasticEnergy;
}

//! @brief calculate the largest element eigenvalue for a group of elements solving the generalized eigenvalue problem Ku=lambda Mu
//! this is used for the estimation of the critical time step
double NuTo::StructureBase::ElementGroupCalculateLargestElementEigenvalue(int rGroupId)
{
    std::vector< NuTo::ElementBase*> elementVector;
    NuTo::GroupBase* grp_PtrBase = this->GroupGetGroupPtr(rGroupId);
    Group<NuTo::ElementBase> *grp_Ptr = grp_PtrBase->AsGroupElement();
    assert(grp_Ptr!=0);
	this->GetElementsByGroup(grp_Ptr,elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}

//! @brief calculate the critical time step for all elements solving the generalized eigenvalue problem Ku=lambda Mu
double NuTo::StructureBase::ElementTotalCalculateLargestElementEigenvalue()
{
    std::vector< ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}


//! @brief calculate the critical time step for a vector of elements solving the generalized eigenvalue problem Ku=lambda Mu
double NuTo::StructureBase::ElementCalculateLargestElementEigenvalue(const std::vector< ElementBase*>& rElementVector)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Error::eError errorGlobal (Error::SUCCESSFUL);

    int exception(0);
    std::string exceptionStringTotal;

	double maxGlobalEigenValue(0);


#ifdef _OPENMP
	if (mNumProcessors!=0)
		omp_set_num_threads(mNumProcessors);
    #pragma omp parallel default(shared)
#endif //_OPENMP
	{

	    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
        elementOutput[Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockVectorDouble>(mDofStatus);
        elementOutput[Element::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(mDofStatus);

		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;


#ifdef _OPENMP
    	#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
		for (unsigned int countElement=0;  countElement<rElementVector.size();countElement++)
		{
			try
			{
				Error::eError error = rElementVector[countElement]->Evaluate(elementOutput);
				if (error!=Error::SUCCESSFUL)
#ifdef _OPENMP
#pragma omp critical
				{
					if (errorGlobal==Error::SUCCESSFUL)
						errorGlobal = error;
					else if (errorGlobal!=error)
						throw MechanicsException("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] multiple elements have returned different error codes, can't handle that.");
				}
#else //_OPENMP
				errorGlobal = error;
#endif //_OPENMP

				auto lumpedMass = elementOutput.at(Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->GetBlockFullVectorDouble().Export();
				auto stiffness = elementOutput.at(Element::HESSIAN_0_TIME_DERIVATIVE)->GetBlockFullMatrixDouble().Export();

				//assuming the stiffness matrix is symmetric
				//std::cout << "lumped mass in element routine\n" <<  lumpedMass << std::endl;
				//std::cout << "element stiffness\n" << stiffness << std::endl;

				//eigenSolver.compute(stiffness);
				//std::cout << "eigenvalues element stiffness\n" << eigenSolver.eigenvalues() << std::endl;

				//std::cout << "eigenmatrix element\n" << lumpedMass.cwiseInverse().asDiagonal()*stiffness << std::endl;

				//invert the lumped mass matrix
				eigenSolver.compute(stiffness,lumpedMass.asDiagonal());
				//std::cout << "eigenvalues in element routine\n" <<  eigenSolver.eigenvalues() << std::endl;

				double maxElementEigenValue = eigenSolver.eigenvalues().maxCoeff();
				if (maxElementEigenValue>maxGlobalEigenValue)
				{
#ifdef _OPENMP
#pragma omp critical
#endif //_OPENMP
					maxGlobalEigenValue = maxElementEigenValue;
				}
			}
			catch(NuTo::Exception& e)
			{
				std::stringstream ss;
				ss << ElementGetId(rElementVector[countElement]);
				std::string exceptionStringLocal(e.ErrorMessage()
						+"[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] Error calculating critical time step for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
			catch(...)
			{
				std::stringstream ss;
				ss << ElementGetId(rElementVector[countElement]);
				std::string exceptionStringLocal("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] Error calculating critical time step for element " + ss.str() + ".\n");
				exception+=1;
				exceptionStringTotal+=exceptionStringLocal;
			}
		}
	} //end of parallel region
    if(exception>0)
    {
	    throw MechanicsException(exceptionStringTotal);
    }
    if (errorGlobal!=Error::SUCCESSFUL)
    {
	    throw MechanicsException("[NuTo::StructureBase::ElementTotalCalculateCriticalTimeStep] error calculating critical time step.");
    }

    return maxGlobalEigenValue;

}

//! @brief calculates the volume of the elements
//! @param rGroupId  group number
//! @return volume of the structure in 3D /area in 2D/ length in 1D
double NuTo::StructureBase::ElementGroupGetVolume(int rGroupId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetVolume] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupGetVolume] Group is not an element group.");
    const Group<ElementBase> *elementGroup = dynamic_cast<const Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    double totalVolume(0);

    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        try
        {
        	totalVolume += itElement->second->GetIntegrationPointVolume().sum();
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupGetVolume] Error calculating volume for element "  + ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               ("[NuTo::StructureBase::ElementGroupGetVolume] Error calculating volume for element " + ss.str() + ".");
        }
    }
    return totalVolume;
}


#ifdef ENABLE_VISUALIZE
//! @brief ... adds all the elements in the vector to the data structure that is finally visualized
void NuTo::StructureBase::ElementTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    std::vector< ElementBase*> elementVec;
    this->GetElementsTotal(elementVec);
    ElementVectorAddToVisualize(rVisualize,rVisualizationList,elementVec);
}

//! @brief ... adds all the elements in a group to the data structure that is finally visualized
void NuTo::StructureBase::ElementGroupAddToVisualize(int rGroupId, VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
	// find group by name
	Group<ElementBase>* elementGroup = this->GroupGetGroupPtr(rGroupId)->AsGroupElement();
	std::vector< ElementBase*> elementVec;
	this->GetElementsByGroup(elementGroup,elementVec);
    ElementVectorAddToVisualize(rVisualize,rVisualizationList,elementVec,mGroupVisualizationType.at(rGroupId));

}

//! @brief ... adds all the elements in the vector to the data structure that is finally visualized
void NuTo::StructureBase::ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements, const VisualizeBase::eVisualizationType rVisualizationType)
{
    // build global tmp static data
    if (mHaveTmpStaticData and mUpdateTmpStaticDataRequired)
    	throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Update of tmpStaticData required first.");

    switch (rVisualizationType)
    {
        case VisualizeBase::VORONOI_CELL:
            for (auto const & iElePtr : rElements)
                iElePtr->Visualize(rVisualize, rVisualizationList);
            break;
        case VisualizeBase::EXTRAPOLATION_TO_NODES:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeExtrapolateToNodes(rVisualize, rVisualizationList);
            break;
        case VisualizeBase::POINTS:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeIntegrationPointData(rVisualize, rVisualizationList);
            break;
        default:
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization type not implemented.");
    }
}





#endif //VISUALIZE
