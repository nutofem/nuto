#include <eigen3/Eigen/Eigenvalues>

#include <cassert>
#include <boost/tokenizer.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

#include "base/Timer.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "mechanics/elements/ElementOutputFullVectorDouble.h"

#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"

#include "mechanics/elements/ElementOutputVectorInt.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/ElementOutputDummy.h"

#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

NuTo::BlockFullVector<double> NuTo::StructureBase::ElementBuildInternalGradient(ElementBase& rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());


    try
    {
        rElement.Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "MechanicsError in element " + std::to_string(ElementGetId(&rElement)));
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(&rElement)));
    }

    return elementOutputMap.at(Element::eOutput::INTERNAL_GRADIENT)->GetBlockFullVectorDouble();
}

NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian(Element::eOutput rHessianType, ElementBase& rElement)
{
    std::set<Element::eOutput> supportedTypes({Element::eOutput::HESSIAN_0_TIME_DERIVATIVE, Element::eOutput::HESSIAN_1_TIME_DERIVATIVE, Element::eOutput::HESSIAN_2_TIME_DERIVATIVE});
    if (supportedTypes.find(rHessianType) == supportedTypes.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "requested matrix type is not supported or not implemented yet.");


    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[rHessianType] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());

    try
    {
        rElement.Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__ , "MechanicsError in element " + std::to_string(ElementGetId(&rElement)));
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(&rElement)));
    }

    return elementOutputMap.at(rHessianType)->GetBlockFullMatrixDouble();
}


NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian0(int rElementId)
{
    return ElementBuildHessian0 (*ElementGetElementPtr(rElementId));
}
NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian1(int rElementId)
{
    return ElementBuildHessian1 (*ElementGetElementPtr(rElementId));
}
NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian2(int rElementId)
{
    return ElementBuildHessian2 (*ElementGetElementPtr(rElementId));
}


NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian0(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE, rElement);
}

NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian1(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_1_TIME_DERIVATIVE, rElement);
}

NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian2(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_2_TIME_DERIVATIVE, rElement);
}

NuTo::BlockFullVector<int> NuTo::StructureBase::ElementBuildGlobalDofsRow(ElementBase& rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    try
    {
        rElement.Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "MechanicsError in element " + std::to_string(ElementGetId(&rElement)));
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(&rElement)));
    }

    return elementOutputMap.at(Element::eOutput::GLOBAL_ROW_DOF)->GetBlockFullVectorInt();
}

NuTo::BlockFullVector<int> NuTo::StructureBase::ElementBuildGlobalDofsColumn(ElementBase& rElement)
{
    std::map<Element::eOutput,std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::GLOBAL_COLUMN_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    try
    {
        rElement.Evaluate(elementOutputMap);

    } catch (MechanicsException& e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "MechanicsError in element " + std::to_string(ElementGetId(&rElement)));
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Non-mechanics error in element " + std::to_string(ElementGetId(&rElement)));
    }

    return elementOutputMap.at(Element::eOutput::GLOBAL_COLUMN_DOF)->GetBlockFullVectorInt();
}


NuTo::BlockFullMatrix<double> NuTo::StructureBase::ElementBuildHessian0_CDF(ElementBase& rElement, double rDelta)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    auto internalGradient0 = ElementBuildInternalGradient(rElement);
    auto globalColumnDofs  = ElementBuildGlobalDofsColumn(rElement);
    auto dofs = rElement.GetInterpolationType().GetActiveDofs();



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

            int numCols = globalColumnDofsCol.rows();
            int numRows = internalGradient0[dofRow].rows();

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
                hessian0_CDF_dof.col(iCol) = (ElementBuildInternalGradient(rElement)[dofRow]-internalGradient0[dofRow])/rDelta;

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

bool NuTo::StructureBase::ElementCheckHessian0(ElementBase& rElement, double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
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
            // TODO: Do not loop over all possible combinations of DOFs but over a list of combinations created by the constitutive law of the corresponding element. What if an element has multiple constitutive laws assigned?
            if(not rElement.GetConstitutiveLaw(0).CheckDofCombinationComputable(dofRow,dofCol,0))
                continue;

            double scaling = hessianRef(dofRow, dofCol).cwiseAbs().maxCoeff();
            Eigen::MatrixXd  differenceRelative = differenceAbsolute(dofRow, dofCol) / scaling;
            int row = 0, col = 0;
            double maxRelativeDifference = differenceRelative.cwiseAbs().maxCoeff(&row, &col);
            if (maxRelativeDifference > rRelativeTolerance)
            {
                isHessianCorrect = false;
                GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 at dof combination [" << Node::DofToString(dofRow) << " - " << Node::DofToString(dofCol) << "]\n";
                GetLogger() << "maxRelativeDifference " << maxRelativeDifference << " at entry (" << row << "," << col << ")\n";
                if (rPrintWrongMatrices)
                {
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
        bool isElementCorrect = ElementCheckHessian0(*elementIdPair.second, rDelta, rRelativeTolerance, rPrintWrongMatrices);
        areAllElementsCorrect = areAllElementsCorrect && isElementCorrect;

        if (not isElementCorrect)
        {
            GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 in Element " << elementIdPair.first << "\n";
            GetLogger() << "################################################################################################\n";
        }

    }
    SetShowTime(showTime);

    return areAllElementsCorrect;
}

void NuTo::StructureBase::ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    ElementSetConstitutiveLaw(elementPtr,itConstitutive->second);
}

void NuTo::StructureBase::ElementGroupSetConstitutiveLaw(int rGroupIdent, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
    	throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

	boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        ElementSetConstitutiveLaw(itElement->second,itConstitutive->second);
    }
}

void NuTo::StructureBase::ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive==mConstitutiveLawMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        ElementSetConstitutiveLaw(elementVector[countElement],itConstitutive->second);
    }
}
void NuTo::StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive)
{
    try
    {
        rElement->SetConstitutiveLaw(*rConstitutive);
    }
    catch(MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error setting constitutive law  for element " + std::to_string(ElementGetId(rElement)) + ".");
        throw;
    }
    catch(...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Error setting constitutive law  for element " + std::to_string(ElementGetId(rElement)) + ".");
    }
}


void NuTo::StructureBase::ElementGroupSetSection(int rGroupIdent, std::shared_ptr<Section> section)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
    	throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        ElementSetSection(itElement->second, section);
    }
}

void NuTo::StructureBase::ElementTotalSetSection(std::shared_ptr<Section> section)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement=0;  countElement<elementVector.size();countElement++)
    {
        ElementSetSection(elementVector[countElement], section);
    }
}


void StructureBase::ElementSetSection(int rElementId, std::shared_ptr<Section> section)
{
    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    ElementSetSection(elementPtr, section);
}


void StructureBase::ElementSetSection(ElementBase* rElement, std::shared_ptr<Section> section)
{
    rElement->SetSection(section);
}


void NuTo::StructureBase::ElementSetInterpolationType(int rElementId, int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type with the given identifier does not exist.");

    ElementSetInterpolationType(elementPtr, itInterpolationType->second);
}

void NuTo::StructureBase::ElementGroupSetInterpolationType(int rGroupId, int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetInterpolationType] Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    boost::ptr_map<int,InterpolationType>::iterator itInterpolationType = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType==mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupSetConstitutiveLaw] Interpolation type with the given identifier does not exist.");

    for (Group<ElementBase>::iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++)
    {
        ElementSetInterpolationType(itElement->second, itInterpolationType->second);
    }
}



void NuTo::StructureBase::ElementSetInterpolationType(ElementBase* rElement, InterpolationType* rInterpolationType)
{
    try
    {
        rElement->SetInterpolationType(*rInterpolationType);
    }
    catch(NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error setting interpolation type for element " + std::to_string(ElementGetId(rElement)) + ".");
        throw;
    }
    catch(...)
    {
        throw NuTo::MechanicsException
            (__PRETTY_FUNCTION__, "Error setting interpolation type for element " + std::to_string(ElementGetId(rElement)) + ".");
    }
}


Eigen::MatrixXd NuTo::StructureBase::ElementGetStaticIPData(int rElementId, std::string rType)
{
    return ElementGetStaticIPData(rElementId, IpData::IpStaticDataTypeToEnum(rType));
}

Eigen::MatrixXd NuTo::StructureBase::ElementGetEngineeringStrain(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_STRAIN);
}

Eigen::MatrixXd NuTo::StructureBase::ElementGetEngineeringPlasticStrain(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN);
}

Eigen::MatrixXd NuTo::StructureBase::ElementGetEngineeringStress(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_STRESS);
}

Eigen::MatrixXd NuTo::StructureBase::ElementGetDamage(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::DAMAGE);
}

Eigen::MatrixXd NuTo::StructureBase::ElementGetStaticIPData(int rElementId, IpData::eIpStaticDataType rType)
{
    Timer timer(std::string(__FUNCTION__) + ":" + IpData::IpStaticDataTypeToString(rType), GetShowTime(), GetLogger());

    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
        elementOutputMap[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(rType);

        elementPtr->Evaluate(elementOutputMap);

        return elementOutputMap.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap().at(rType);

    } catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error getting " + IpData::IpStaticDataTypeToString(rType) + " for element " + std::to_string(rElementId) + ".");
        throw;
    } catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error getting " + IpData::IpStaticDataTypeToString(rType) + " for element "  + std::to_string(rElementId)  + ".");
    }
}


Eigen::MatrixXd NuTo::StructureBase::ElementGetIntegrationPointCoordinates(int rElementId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");
    }



    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    try
    {
		//evaluate the coordinates
        Eigen::MatrixXd coordinates (3,elementPtr->GetNumIntegrationPoints());
    	for (int count=0; count<elementPtr->GetNumIntegrationPoints(); count++)
    	{
    	    Eigen::Vector3d coords = elementPtr->GetGlobalIntegrationPointCoordinates(count);
    	    coordinates.block<3,1>(0, count) = coords;
    	}
    	return coordinates;
    }
    catch(NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__, "Error getting integration point coordinates for for element " + std::to_string(rElementId)  + ".");
        throw;
    }
    catch(...)
    {
    	throw NuTo::MechanicsException
    	   (__PRETTY_FUNCTION__, "Error getting integration point coordinates for for element " + std::to_string(rElementId)  + ".");
    }
}

double NuTo::StructureBase::ElementTotalGetMaxDamage()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::DAMAGE);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    double maxDamage = 0;

    Eigen::MatrixXd rIPDamage;

    for (auto element : elementVector)
    {
        try
        {
            element->Evaluate(elementOutputMap);
            rIPDamage = elementOutputMap.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap()[IpData::eIpStaticDataType::DAMAGE];
        } catch (NuTo::MechanicsException &e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, "Error getting damage for element " + std::to_string(ElementGetId(element)) + ".");
            throw;
        } catch (...)
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Error getting damage for element " +std::to_string(ElementGetId(element)) + ".");
        }
        maxDamage = std::max(maxDamage, rIPDamage.maxCoeff());
    }
    return maxDamage;
}


double StructureBase::ElementTotalGetStaticDataExtrapolationError()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::IP_DATA] =
            std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::EXTRAPOLATION_ERROR);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    double maxError = 0;
    Eigen::MatrixXd ipValues;
    for (auto element : elementVector)
    {
        element->Evaluate(elementOutputMap);
        ipValues = elementOutputMap.at(Element::eOutput::IP_DATA)
                           ->GetIpData()
                           .GetIpDataMap()[IpData::eIpStaticDataType::EXTRAPOLATION_ERROR];
        maxError = std::max(maxError, ipValues.maxCoeff());
    }
    return maxError;
}


void NuTo::StructureBase::ElementGroupAllocateAdditionalStaticData(int rElementGroupId, int rNumAdditionalStaticData)
{
    if (GroupGetGroupPtr(rElementGroupId)->GetType() != eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Element group required.");

    for (int elementId : GroupGetMemberIds(rElementGroupId))
    {
        ElementBase* element = ElementGetElementPtr(elementId);
        IPData& ipdata = element->GetIPData();
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            ipdata.GetIPConstitutiveLaw(i).AllocateAdditional(rNumAdditionalStaticData);
    }

}

void NuTo::StructureBase::ElementTotalUpdateStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        this->ElementTotalUpdateTmpStaticData();
    }

	std::vector<ElementBase*> elementVector;
	GetElementsTotal(elementVector);
    int exception(0);
    std::string exceptionStringTotal;

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();

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
            elementVector[countElement]->Evaluate(elementOutput);
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
}

void NuTo::StructureBase::ElementTotalUpdateTmpStaticData()
{
    if (mHaveTmpStaticData)
	{
		std::vector<ElementBase*> elementVector;
		GetElementsTotal(elementVector);
	    int exception(0);
	    std::string exceptionStringTotal;

	    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
        elementOutput[Element::eOutput::UPDATE_TMP_STATIC_DATA] = std::make_shared<ElementOutputDummy>();


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
				elementVector[countElement]->Evaluate(elementOutput);
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
}

void NuTo::StructureBase::ElementTotalShiftStaticDataToPast()
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
            IPData& ipdata = element->GetIPData();
            for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
                ipdata.GetIPConstitutiveLaw(i).ShiftToPast();
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


void NuTo::StructureBase::ElementTotalShiftStaticDataToFuture()
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
            IPData& ipdata = element->GetIPData();
            for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
                ipdata.GetIPConstitutiveLaw(i).ShiftToFuture();
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

void NuTo::StructureBase::ElementTotalExtrapolateStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();

    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

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


void NuTo::StructureBase::ElementTotalGetAverageStress(double rVolume, Eigen::MatrixXd& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Eigen::MatrixXd elementEngineeringStress;
    rEngineeringStress.resize(6,1);

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
            throw;
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

void NuTo::StructureBase::ElementGroupGetAverageStress(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    Eigen::MatrixXd elementEngineeringStress;
    rEngineeringStress.resize(6,1);

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
            e.AddMessage(__PRETTY_FUNCTION__, "Error calculating integrated stress  for element "  + ss.str() + ".");
            throw;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               (__PRETTY_FUNCTION__, "Error calculating integrated stress  for element " + ss.str() + ".");
        }
    }
    rEngineeringStress*=(1./rVolume);

}


void NuTo::StructureBase::ElementTotalGetAverageStrain(double rVolume, Eigen::MatrixXd& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Eigen::MatrixXd elementEngineeringStrain;
    rEngineeringStrain.resize(6,1);

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
            e.AddMessage(__PRETTY_FUNCTION__, "Error calculating integrated strain  for element "  + ss.str() + ".");
            throw;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[elementCount]);
            throw NuTo::MechanicsException
               (__PRETTY_FUNCTION__, "Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=1./rVolume;
}

void NuTo::StructureBase::ElementGroupGetAverageStrain(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
    Group<ElementBase> *elementGroup = dynamic_cast<Group<ElementBase>*>(itGroup->second);
    assert(elementGroup!=0);

    Eigen::MatrixXd elementEngineeringStrain;
    rEngineeringStrain.resize(6,1);

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
            e.AddMessage(__PRETTY_FUNCTION__, "Error calculating integrated strain  for element "  + ss.str() + ".");
            throw;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               (__PRETTY_FUNCTION__, "Error calculating integrated strain  for element " + ss.str() + ".");
        }
    }
    rEngineeringStrain*=(1./rVolume);

}


void NuTo::StructureBase::ElementGroupGetMembers(int rGroupId, std::vector<int>& rMembers)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
    Group<ElementBase> *elementGroup = itGroup->second->AsGroupElement();
    assert(elementGroup!=0);

    rMembers.resize(elementGroup->GetNumMembers());
    int countElement(0);
    for (Group<ElementBase>::const_iterator itElement=elementGroup->begin(); itElement!=elementGroup->end();itElement++,countElement++)
    {
       	rMembers[countElement] = itElement->first;
    }
}


double NuTo::StructureBase::ElementGroupCalculateLargestElementEigenvalue(int rGroupId)
{
    std::vector< NuTo::ElementBase*> elementVector;
    NuTo::GroupBase* grp_PtrBase = this->GroupGetGroupPtr(rGroupId);
    Group<NuTo::ElementBase> *grp_Ptr = grp_PtrBase->AsGroupElement();
    assert(grp_Ptr!=0);
	this->GetElementsByGroup(grp_Ptr,elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}

double NuTo::StructureBase::ElementTotalCalculateLargestElementEigenvalue()
{
    std::vector< ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}


double NuTo::StructureBase::ElementCalculateLargestElementEigenvalue(const std::vector< ElementBase*>& rElementVector)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

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
        elementOutput[Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
        elementOutput[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());

		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;


#ifdef _OPENMP
    	#pragma omp for schedule(dynamic,1) nowait
#endif //_OPENMP
		for (unsigned int countElement=0;  countElement<rElementVector.size();countElement++)
		{
			try
			{
				rElementVector[countElement]->Evaluate(elementOutput);

                auto lumpedMass = elementOutput.at(Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE)->GetBlockFullVectorDouble().Export();
                auto stiffness = elementOutput.at(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE)->GetBlockFullMatrixDouble().Export();

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
    return maxGlobalEigenValue;
}

double NuTo::StructureBase::ElementGroupGetVolume(int rGroupId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Elements)
        throw MechanicsException(__PRETTY_FUNCTION__, "Group is not an element group.");
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
            e.AddMessage(__PRETTY_FUNCTION__, "Error calculating volume for element "  + ss.str() + ".");
            throw;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            throw NuTo::MechanicsException
               (__PRETTY_FUNCTION__, "Error calculating volume for element " + ss.str() + ".");
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
void NuTo::StructureBase::ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements)
{
    ElementVectorAddToVisualize(rVisualize,rVisualizationList,rElements,eVisualizationType::VORONOI_CELL);
}

//! @brief ... adds all the elements in the vector to the data structure that is finally visualized
void NuTo::StructureBase::ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements, const eVisualizationType rVisualizationType)
{
    // build global tmp static data
    if (mHaveTmpStaticData and mUpdateTmpStaticDataRequired)
    	throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Update of tmpStaticData required first.");

    switch (rVisualizationType)
    {
        case eVisualizationType::VORONOI_CELL:
            for (auto const & iElePtr : rElements)
                iElePtr->Visualize(rVisualize, rVisualizationList);
            break;
        case eVisualizationType::EXTRAPOLATION_TO_NODES:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeExtrapolateToNodes(rVisualize, rVisualizationList);
            break;
        case eVisualizationType::POINTS:
            for (auto const & iElePtr : rElements)
                iElePtr->VisualizeIntegrationPointData(rVisualize, rVisualizationList);
            break;
        default:
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization type not implemented.");
    }
}





#endif //VISUALIZE
