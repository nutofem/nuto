#include <cassert>
#include <eigen3/Eigen/Eigenvalues>

#include "mechanics/structures/StructureBase.h"

#include "base/Timer.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/ElementOutputDummy.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/nodes/NodeEnum.h"

#include "visualize/VisualizeEnum.h"

using namespace NuTo;

BlockFullVector<double> StructureBase::ElementBuildInternalGradient(ElementBase& rElement)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] =
            std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());

    rElement.Evaluate(elementOutputMap);
    return elementOutputMap.at(Element::eOutput::INTERNAL_GRADIENT)->GetBlockFullVectorDouble();
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian(Element::eOutput rHessianType, ElementBase& rElement)
{
    std::set<Element::eOutput> supportedTypes({Element::eOutput::HESSIAN_0_TIME_DERIVATIVE,
                                               Element::eOutput::HESSIAN_1_TIME_DERIVATIVE,
                                               Element::eOutput::HESSIAN_2_TIME_DERIVATIVE});
    if (supportedTypes.find(rHessianType) == supportedTypes.end())
        throw Exception(__PRETTY_FUNCTION__, "requested matrix type is not supported or not implemented yet.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[rHessianType] = std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());

    rElement.Evaluate(elementOutputMap);

    return elementOutputMap.at(rHessianType)->GetBlockFullMatrixDouble();
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian0(int rElementId)
{
    return ElementBuildHessian0(*ElementGetElementPtr(rElementId));
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian1(int rElementId)
{
    return ElementBuildHessian1(*ElementGetElementPtr(rElementId));
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian2(int rElementId)
{
    return ElementBuildHessian2(*ElementGetElementPtr(rElementId));
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian0(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE, rElement);
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian1(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_1_TIME_DERIVATIVE, rElement);
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian2(ElementBase& rElement)
{
    return ElementBuildHessian(Element::eOutput::HESSIAN_2_TIME_DERIVATIVE, rElement);
}


BlockFullVector<int> StructureBase::ElementBuildGlobalDofsRow(ElementBase& rElement)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] = std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    rElement.Evaluate(elementOutputMap);

    return elementOutputMap.at(Element::eOutput::GLOBAL_ROW_DOF)->GetBlockFullVectorInt();
}


BlockFullVector<int> StructureBase::ElementBuildGlobalDofsColumn(ElementBase& rElement)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::GLOBAL_COLUMN_DOF] =
            std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());

    rElement.Evaluate(elementOutputMap);

    return elementOutputMap.at(Element::eOutput::GLOBAL_COLUMN_DOF)->GetBlockFullVectorInt();
}


BlockFullMatrix<double> StructureBase::ElementBuildHessian0_CDF(ElementBase& rElement, double rDelta)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    auto internalGradient0 = ElementBuildInternalGradient(rElement);
    auto globalColumnDofs = ElementBuildGlobalDofsColumn(rElement);
    auto dofs = rElement.GetInterpolationType().GetActiveDofs();


    BlockFullMatrix<double> hessian0_CDF(GetDofStatus());

    auto dofValues = this->NodeExtractDofValues(0);

    for (auto dofRow : dofs)
    {
        for (auto dofCol : dofs)
        {
            auto& dofActDofValues = dofValues.J[dofCol];
            auto& dofDepDofValues = dofValues.K[dofCol];
            auto& hessian0_CDF_dof = hessian0_CDF(dofRow, dofCol);
            auto& globalColumnDofsCol = globalColumnDofs[dofCol];

            int numCols = globalColumnDofsCol.rows();
            int numRows = internalGradient0[dofRow].rows();

            hessian0_CDF_dof.resize(numRows, numCols);
            hessian0_CDF_dof.setZero();
            for (int iCol = 0; iCol < numCols; ++iCol)
            {
                // Apply rDelta to the corresponding Dof
                int index = globalColumnDofsCol(iCol);
                int numActiveDofs = GetNumActiveDofs(dofCol);

                if (index < numActiveDofs)
                    dofActDofValues(index) += rDelta;
                else
                    dofDepDofValues(index - numActiveDofs) += rDelta;

                NodeMergeDofValues(0, dofValues);

                // Calculate CDF
                hessian0_CDF_dof.col(iCol) =
                        (ElementBuildInternalGradient(rElement)[dofRow] - internalGradient0[dofRow]) / rDelta;

                // Restore the orignial dof values
                if (index < numActiveDofs)
                    dofActDofValues(index) -= rDelta;
                else
                    dofDepDofValues(index - numActiveDofs) -= rDelta;
            }
            this->NodeMergeDofValues(0, dofValues);
        }
    }
    return hessian0_CDF;
}


bool StructureBase::ElementCheckHessian0(ElementBase& rElement, double rDelta, double rRelativeTolerance,
                                         bool rPrintWrongMatrices)
{
    bool isHessianCorrect = true;

    auto hessianRef = ElementBuildHessian0_CDF(rElement, rDelta);
    auto hessian = ElementBuildHessian0(rElement);

    auto differenceAbsolute = hessianRef - hessian;

    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");

    for (auto dofRow : GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : GetDofStatus().GetActiveDofTypes())
        {
            // TODO: Do not loop over all possible combinations of DOFs but over a list of combinations created by the
            // constitutive law of the corresponding element. What if an element has multiple constitutive laws
            // assigned?
            if (not rElement.GetConstitutiveLaw(0).CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;

            double scaling = hessianRef(dofRow, dofCol).cwiseAbs().maxCoeff();
            Eigen::MatrixXd differenceRelative = differenceAbsolute(dofRow, dofCol) / scaling;
            int row = 0, col = 0;
            double maxRelativeDifference = differenceRelative.cwiseAbs().maxCoeff(&row, &col);
            if (maxRelativeDifference > rRelativeTolerance)
            {
                isHessianCorrect = false;
                GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 at dof combination ["
                            << Node::DofToString(dofRow) << " - " << Node::DofToString(dofCol) << "]\n";
                GetLogger() << "maxRelativeDifference " << maxRelativeDifference << " at entry (" << row << "," << col
                            << ")\n";
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


bool StructureBase::ElementCheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    bool showTime = GetShowTime();
    SetShowTime(false);

    std::vector<std::pair<int, ElementBase*>> elements;
    GetElementsTotal(elements);

    bool areAllElementsCorrect = true;

    for (auto elementIdPair : elements)
    {
        bool isElementCorrect =
                ElementCheckHessian0(*elementIdPair.second, rDelta, rRelativeTolerance, rPrintWrongMatrices);
        areAllElementsCorrect = areAllElementsCorrect && isElementCorrect;

        if (not isElementCorrect)
        {
            GetLogger() << "[" << __FUNCTION__ << "] wrong hessian0 in Element " << elementIdPair.first << "\n";
            GetLogger() << "#######################################################################################\n";
        }
    }
    SetShowTime(showTime);

    return areAllElementsCorrect;
}


void StructureBase::ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);

    boost::ptr_map<int, ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive == mConstitutiveLawMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    ElementSetConstitutiveLaw(elementPtr, itConstitutive->second);
}


void StructureBase::ElementGroupSetConstitutiveLaw(int rGroupIdent, int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    boost::ptr_map<int, ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive == mConstitutiveLawMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    for (auto& element : elementGroup)
    {
        ElementSetConstitutiveLaw(element.second, itConstitutive->second);
    }
}


void StructureBase::ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, ConstitutiveBase>::iterator itConstitutive = mConstitutiveLawMap.find(rConstitutiveLawIdent);
    if (itConstitutive == mConstitutiveLawMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law with the given identifier does not exist.");

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
    {
        ElementSetConstitutiveLaw(elementVector[countElement], itConstitutive->second);
    }
}


void StructureBase::ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive)
{
    rElement->SetConstitutiveLaw(*rConstitutive);
}


void StructureBase::ElementGroupSetSection(int rGroupIdent, std::shared_ptr<Section> section)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    for (auto& element : elementGroup)
    {
        ElementSetSection(element.second, section);
    }
}


void StructureBase::ElementTotalSetSection(std::shared_ptr<Section> section)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
    {
        ElementSetSection(elementVector[countElement], section);
    }
}


void StructureBase::ElementSetSection(ElementBase* rElement, std::shared_ptr<Section> section)
{
    rElement->SetSection(section);
}


void StructureBase::ElementGroupSetInterpolationType(int rGroupId, int rInterpolationTypeId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup == mGroupMap.end())
        throw Exception(
                "[StructureBase::ElementGroupSetInterpolationType] Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception("[StructureBase::ElementGroupSetInterpolationType] Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    boost::ptr_map<int, InterpolationType>::iterator itInterpolationType =
            mInterpolationTypeMap.find(rInterpolationTypeId);
    if (itInterpolationType == mInterpolationTypeMap.end())
        throw Exception("[StructureBase::ElementGroupSetConstitutiveLaw] Interpolation type with the given "
                                 "identifier does not exist.");

    for (auto& element : elementGroup)
    {
        ElementSetInterpolationType(element.second, itInterpolationType->second);
    }
}


void StructureBase::ElementSetInterpolationType(ElementBase* rElement, InterpolationType* rInterpolationType)
{
    rElement->SetInterpolationType(*rInterpolationType);
}


Eigen::MatrixXd StructureBase::ElementGetStaticIPData(int rElementId, std::string rType)
{
    return ElementGetStaticIPData(rElementId, IpData::IpStaticDataTypeToEnum(rType));
}


Eigen::MatrixXd StructureBase::ElementGetEngineeringStrain(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_STRAIN);
}


Eigen::MatrixXd StructureBase::ElementGetEngineeringPlasticStrain(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN);
}


Eigen::MatrixXd StructureBase::ElementGetEngineeringStress(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::ENGINEERING_STRESS);
}


Eigen::MatrixXd StructureBase::ElementGetDamage(int rElementId)
{
    return ElementGetStaticIPData(rElementId, IpData::eIpStaticDataType::DAMAGE);
}


Eigen::MatrixXd StructureBase::ElementGetStaticIPData(int rElementId, IpData::eIpStaticDataType rType)
{
    Timer timer(std::string(__FUNCTION__) + ":" + IpData::IpStaticDataTypeToString(rType), GetShowTime(), GetLogger());

    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw Exception(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(rType);

    elementPtr->Evaluate(elementOutputMap);

    return elementOutputMap.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap().at(rType);
}


Eigen::MatrixXd StructureBase::ElementGetIntegrationPointCoordinates(int rElementId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
    {
        throw Exception(__PRETTY_FUNCTION__, "First update of tmp static data required.");
    }

    ElementBase* elementPtr = ElementGetElementPtr(rElementId);
    Eigen::MatrixXd coordinates(3, elementPtr->GetNumIntegrationPoints());
    for (int count = 0; count < elementPtr->GetNumIntegrationPoints(); count++)
    {
        Eigen::Vector3d coords = elementPtr->GetGlobalIntegrationPointCoordinates(count);
        coordinates.block<3, 1>(0, count) = coords;
    }
    return coordinates;
}


double StructureBase::ElementTotalGetMaxDamage()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw Exception(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;
    elementOutputMap[Element::eOutput::IP_DATA] =
            std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::DAMAGE);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    double maxDamage = 0;
    Eigen::MatrixXd rIPDamage;
    for (auto element : elementVector)
    {
        element->Evaluate(elementOutputMap);
        rIPDamage = elementOutputMap.at(Element::eOutput::IP_DATA)
                            ->GetIpData()
                            .GetIpDataMap()[IpData::eIpStaticDataType::DAMAGE];
        maxDamage = std::max(maxDamage, rIPDamage.maxCoeff());
    }
    return maxDamage;
}


double StructureBase::ElementTotalGetStaticDataExtrapolationError()
{
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw Exception(__PRETTY_FUNCTION__, "First update of tmp static data required.");

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


void StructureBase::ElementGroupAllocateAdditionalStaticData(int rElementGroupId, int rNumAdditionalStaticData)
{
    if (GroupGetGroupPtr(rElementGroupId)->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Element group required.");

    for (int elementId : GroupGetMemberIds(rElementGroupId))
    {
        ElementBase* element = ElementGetElementPtr(elementId);
        IPData& ipdata = element->GetIPData();
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            ipdata.GetIPConstitutiveLaw(i).AllocateAdditional(rNumAdditionalStaticData);
    }
}


void StructureBase::ElementTotalUpdateStaticData()
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
    if (mNumProcessors != 0)
        omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
    for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
    {
        try
        {
            elementVector[countElement]->Evaluate(elementOutput);
        }
        catch (Exception& e)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            std::string exceptionStringLocal(
                    e.ErrorMessage() +
                    "[StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() +
                    ".\n");
            exception += 1;
            exceptionStringTotal += exceptionStringLocal;
        }
        catch (...)
        {
            std::stringstream ss;
            ss << ElementGetId(elementVector[countElement]);
            std::string exceptionStringLocal(
                    "[StructureBase::ElementTotalUpdateStaticData] Error updating static data for element " + ss.str() +
                    ".\n");
            exception += 1;
            exceptionStringTotal += exceptionStringLocal;
        }
    }
    if (exception > 0)
    {
        throw Exception(exceptionStringTotal);
    }
}


void StructureBase::ElementTotalUpdateTmpStaticData()
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
        if (mNumProcessors != 0)
            omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
        for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
        {
            try
            {
                elementVector[countElement]->Evaluate(elementOutput);
            }
            catch (Exception& e)
            {
                std::stringstream ss;
                ss << ElementGetId(elementVector[countElement]);
                std::string exceptionStringLocal(
                        e.ErrorMessage() +
                        "[StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " +
                        ss.str() + ".\n");
                exception += 1;
                exceptionStringTotal += exceptionStringLocal;
            }
            catch (...)
            {
                std::stringstream ss;
                ss << ElementGetId(elementVector[countElement]);
                std::string exceptionStringLocal(
                        "[StructureBase::ElementTotalUpdateTmpStaticData] exception updating static data for element " +
                        ss.str() + ".\n");
                exception += 1;
                exceptionStringTotal += exceptionStringLocal;
            }
        }
        if (exception > 0)
        {
            throw Exception(exceptionStringTotal);
        }
    }
    mUpdateTmpStaticDataRequired = false;
}


void StructureBase::ElementTotalShiftStaticDataToPast()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

#ifdef _OPENMP
    if (mNumProcessors != 0)
        omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        ElementBase* element = elementVector[iElement];
        IPData& ipdata = element->GetIPData();
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            ipdata.GetIPConstitutiveLaw(i).ShiftToPast();
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}


void StructureBase::ElementTotalShiftStaticDataToFuture()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

#ifdef _OPENMP
    if (mNumProcessors != 0)
        omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        ElementBase* element = elementVector[iElement];
        IPData& ipdata = element->GetIPData();
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            ipdata.GetIPConstitutiveLaw(i).ShiftToFuture();
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}


void StructureBase::ElementTotalExtrapolateStaticData()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    Exception exception("");

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();

    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

#ifdef _OPENMP
    if (mNumProcessors != 0)
        omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
    for (unsigned int iElement = 0; iElement < elementVector.size(); iElement++)
    {
        ElementBase* element = elementVector[iElement];
        element->Evaluate(input, elementOutput);
    }
    if (not exception.ErrorMessage().empty())
        throw exception;
}


void StructureBase::ElementTotalGetAverageStress(double rVolume, Eigen::MatrixXd& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Eigen::MatrixXd elementEngineeringStress;
    rEngineeringStress.resize(6, 1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount = 0; elementCount < elementVector.size(); elementCount++)
    {
        elementVector[elementCount]->GetIntegratedStress(elementEngineeringStress);
        rEngineeringStress += elementEngineeringStress;
    }
    rEngineeringStress *= 1. / rVolume;
}


void StructureBase::ElementGroupGetAverageStress(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStress)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    Eigen::MatrixXd elementEngineeringStress;
    rEngineeringStress.resize(6, 1);

    for (auto& element : elementGroup)
    {
        element.second->GetIntegratedStress(elementEngineeringStress);
        rEngineeringStress += elementEngineeringStress;
    }
    rEngineeringStress *= (1. / rVolume);
}


void StructureBase::ElementTotalGetAverageStrain(double rVolume, Eigen::MatrixXd& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Eigen::MatrixXd elementEngineeringStrain;
    rEngineeringStrain.resize(6, 1);

    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    for (unsigned int elementCount = 0; elementCount < elementVector.size(); elementCount++)
    {
        elementVector[elementCount]->GetIntegratedStrain(elementEngineeringStrain);
        rEngineeringStrain += elementEngineeringStrain;
    }
    rEngineeringStrain *= 1. / rVolume;
}


void StructureBase::ElementGroupGetAverageStrain(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStrain)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    Eigen::MatrixXd elementEngineeringStrain;
    rEngineeringStrain.resize(6, 1);

    for (const auto& element : elementGroup)
    {
        element.second->GetIntegratedStrain(elementEngineeringStrain);
        rEngineeringStrain += elementEngineeringStrain;
    }
    rEngineeringStrain *= (1. / rVolume);
}


void StructureBase::ElementGroupGetMembers(int groupId, std::vector<int>& members)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::iterator itGroup = mGroupMap.find(groupId);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    auto& elementGroup = dynamic_cast<Group<ElementBase>&>(*itGroup->second);

    members.clear();
    for (const auto& element : elementGroup)
    {
        members.push_back(element.first);
    }
}


double StructureBase::ElementGroupCalculateLargestElementEigenvalue(int rGroupId)
{
    std::vector<ElementBase*> elementVector;
    GroupBase* grp_PtrBase = this->GroupGetGroupPtr(rGroupId);
    auto grp_Ptr = dynamic_cast<Group<ElementBase>*>(grp_PtrBase);
    assert(grp_Ptr != nullptr);
    this->GetElementsByGroup(grp_Ptr, elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}


double StructureBase::ElementTotalCalculateLargestElementEigenvalue()
{
    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);
    return this->ElementCalculateLargestElementEigenvalue(elementVector);
}


double StructureBase::ElementCalculateLargestElementEigenvalue(const std::vector<ElementBase*>& rElementVector)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    int exception(0);
    std::string exceptionStringTotal;

    double maxGlobalEigenValue(0);


#ifdef _OPENMP
    if (mNumProcessors != 0)
        omp_set_num_threads(mNumProcessors);
#pragma omp parallel default(shared)
#endif //_OPENMP
    {

        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
        elementOutput[Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE] =
                std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
        elementOutput[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] =
                std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> eigenSolver;


#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1) nowait
#endif //_OPENMP
        for (unsigned int countElement = 0; countElement < rElementVector.size(); countElement++)
        {
            try
            {
                rElementVector[countElement]->Evaluate(elementOutput);

                auto lumpedMass = elementOutput.at(Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE)
                                          ->GetBlockFullVectorDouble()
                                          .Export();
                auto stiffness = elementOutput.at(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE)
                                         ->GetBlockFullMatrixDouble()
                                         .Export();

                // invert the lumped mass matrix
                eigenSolver.compute(stiffness, lumpedMass.asDiagonal());

                double maxElementEigenValue = eigenSolver.eigenvalues().maxCoeff();
                if (maxElementEigenValue > maxGlobalEigenValue)
                {
#ifdef _OPENMP
#pragma omp critical
#endif //_OPENMP
                    maxGlobalEigenValue = maxElementEigenValue;
                }
            }
            catch (Exception& e)
            {
                std::stringstream ss;
                ss << ElementGetId(rElementVector[countElement]);
                std::string exceptionStringLocal(e.ErrorMessage() + "[StructureBase::"
                                                                    "ElementTotalCalculateCriticalTimeStep] Error "
                                                                    "calculating critical time step for element " +
                                                 ss.str() + ".\n");
                exception += 1;
                exceptionStringTotal += exceptionStringLocal;
            }
            catch (...)
            {
                std::stringstream ss;
                ss << ElementGetId(rElementVector[countElement]);
                std::string exceptionStringLocal("[StructureBase::ElementTotalCalculateCriticalTimeStep] Error "
                                                 "calculating critical time step for element " +
                                                 ss.str() + ".\n");
                exception += 1;
                exceptionStringTotal += exceptionStringLocal;
            }
        }
    } // end of parallel region
    if (exception > 0)
    {
        throw Exception(exceptionStringTotal);
    }
    return maxGlobalEigenValue;
}

double StructureBase::ElementGroupGetVolume(int rGroupId)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    boost::ptr_map<int, GroupBase>::const_iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Group with the given identifier does not exist.");
    if (itGroup->second->GetType() != eGroupId::Elements)
        throw Exception(__PRETTY_FUNCTION__, "Group is not an element group.");
    const auto& elementGroup = dynamic_cast<const Group<ElementBase>&>(*itGroup->second);

    double totalVolume(0);

    for (auto& element : elementGroup)
    {
        totalVolume += element.second->GetIntegrationPointVolume().sum();
    }
    return totalVolume;
}


#ifdef ENABLE_VISUALIZE
void StructureBase::ElementTotalAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                               const std::vector<eVisualizeWhat>& visualizeComponents)
{
    std::vector<ElementBase*> elementVec;
    this->GetElementsTotal(elementVec);
    ElementVectorAddToVisualize(visualizer, visualizeComponents, elementVec);
}


void NuTo::StructureBase::ElementGroupAddToVisualize(int rGroupId, Visualize::UnstructuredGrid& visualizer,
                                                     const std::vector<eVisualizeWhat>& visualizeComponents)
{
    // find group by name
    Group<ElementBase>* elementGroup = this->GroupGetGroupPtr(rGroupId)->AsGroupElement();
    std::vector<ElementBase*> elementVec;
    this->GetElementsByGroup(elementGroup, elementVec);
    ElementVectorAddToVisualize(visualizer, visualizeComponents, elementVec, mGroupVisualizationType.at(rGroupId));
}


void StructureBase::ElementVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                                const std::vector<eVisualizeWhat>& visualizeComponents,
                                                const std::vector<ElementBase*>& elements)
{
    ElementVectorAddToVisualize(visualizer, visualizeComponents, elements, eVisualizationType::VORONOI_CELL);
}


void StructureBase::ElementVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                                const std::vector<eVisualizeWhat>& visualizeComponents,
                                                const std::vector<ElementBase*>& elements,
                                                const eVisualizationType rVisualizationType)
{
    if (mHaveTmpStaticData and mUpdateTmpStaticDataRequired)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Update of tmpStaticData required first.");

    switch (rVisualizationType)
    {
    case eVisualizationType::VORONOI_CELL:
        for (auto const& iElePtr : elements)
            iElePtr->Visualize(visualizer, visualizeComponents);
        break;
    case eVisualizationType::EXTRAPOLATION_TO_NODES:
        for (auto const& iElePtr : elements)
            iElePtr->VisualizeExtrapolateToNodes(visualizer, visualizeComponents);
        break;
    case eVisualizationType::POINTS:
        for (auto const& iElePtr : elements)
            iElePtr->VisualizeIntegrationPointData(visualizer, visualizeComponents);
        break;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Visualization type not implemented.");
    }
}
#endif // VISUALIZE
