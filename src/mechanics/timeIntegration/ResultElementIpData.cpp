
#include <boost/assign/ptr_map_inserter.hpp>

#include "mechanics/timeIntegration/ResultElementIpData.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultElementIpData::ResultElementIpData(const std::string& rIdent, int rElementId, NuTo::IpData::eIpStaticDataType rIpDataType) :
        ResultBase(rIdent),
        mElementId(rElementId),
        mIpDataType(rIpDataType)
{
}

void NuTo::ResultElementIpData::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
    assert(rTimeStepPlot>=0);
    Eigen::Matrix<double, 1, Eigen::Dynamic> ipValues(1,this->GetNumData(rStructure));
    this->CalculateValues(rStructure,ipValues);
    if (rTimeStepPlot>=mData.rows())
    {
        this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
    }
    if (ipValues.cols()!=mData.cols())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: The allocated number of columns is wrong.");
    mData.row(rTimeStepPlot) = ipValues;
}

void NuTo::ResultElementIpData::CalculateValues(const StructureBase& rStructure, Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues)const
{
    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));

    std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(mIpDataType);

    const_cast<ElementBase*>(element)->Evaluate(elementOutput);

    const auto& ipDataResult = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap()[mIpDataType];

    // iterate over all ips
    assert(ipDataResult.cols() == element->GetNumIntegrationPoints());
    unsigned int numComponents = ipDataResult.rows();
    for (int iCol = 0; iCol < ipDataResult.cols(); ++iCol)
    {
        rValues.block(0, iCol * numComponents, 1, ipDataResult.rows()) = ipDataResult.col(iCol).transpose();
    }
}

int NuTo::ResultElementIpData::GetNumData(const StructureBase& rStructure) const
{
    unsigned int numComponents;

    switch (mIpDataType)
    {
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS:
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN:
    case NuTo::IpData::eIpStaticDataType::SHRINKAGE_STRAIN:
        numComponents = 6;
        break;
    case NuTo::IpData::eIpStaticDataType::DAMAGE:
        numComponents = 1;
        break;
    case NuTo::IpData::eIpStaticDataType::SLIP:
    case NuTo::IpData::eIpStaticDataType::BOND_STRESS:
        numComponents = 3;
            break;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Ip data type not supported yet.");
    } // switch

    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));
    return element->GetNumIntegrationPoints() * numComponents;
}


NuTo::eTimeIntegrationResultType NuTo::ResultElementIpData::GetResultType() const
{
    switch (mIpDataType)
    {
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS:
        return NuTo::eTimeIntegrationResultType::ELEMENT_IP_STRESS;
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN:
        return NuTo::eTimeIntegrationResultType::ELEMENT_IP_STRAIN;
    case NuTo::IpData::eIpStaticDataType::DAMAGE:
        return NuTo::eTimeIntegrationResultType::ELEMENT_IP_DAMAGE;
    case NuTo::IpData::eIpStaticDataType::BOND_STRESS:
        return NuTo::eTimeIntegrationResultType::ELEMENT_IP_BOND_STRESS;
    case NuTo::IpData::eIpStaticDataType::SLIP:
        return NuTo::eTimeIntegrationResultType::ELEMENT_IP_SLIP;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Ip data type not supported yet.");
    } // switch

}

void NuTo::ResultElementIpData::Info() const
{
    std::cout << "ResultElementIpData Info:      " << std::endl;
    std::cout << "Integration point data type:   " << IpData::IpStaticDataTypeToString(mIpDataType) << std::endl;
    std::cout << "Element id:                    " << mElementId << std::endl;
}


