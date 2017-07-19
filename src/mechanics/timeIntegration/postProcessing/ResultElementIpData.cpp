#include "mechanics/timeIntegration/postProcessing/ResultElementIpData.h"

#include <boost/assign/ptr_map_inserter.hpp>

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/structures/StructureBase.h"

using namespace NuTo;

ResultElementIpData::ResultElementIpData(const std::string& rIdent, int rElementId,
                                         IpData::eIpStaticDataType rIpDataType)
    : ResultBase(rIdent)
    , mElementId(rElementId)
    , mIpDataType(rIpDataType)
{
}


void ResultElementIpData::CalculateAndAddValues(const StructureBase& structure, int timeStep,
                                                const StructureOutputBlockVector& residual, double currentTime)
{
    assert(timeStep >= 0);
    Eigen::Matrix<double, 1, Eigen::Dynamic> ipValues(1, GetNumData(structure));
    CalculateValues(structure, ipValues);
    if (timeStep >= mData.rows())
    {
        Resize(structure, 2 * (timeStep + 1), false);
    }
    if (ipValues.cols() != mData.cols())
        throw Exception(__PRETTY_FUNCTION__, "The allocated number of columns is wrong.");
    mData.row(timeStep) = ipValues;
}

void ResultElementIpData::CalculateValues(const StructureBase& rStructure,
                                          Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues) const
{
    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));

    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
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


int ResultElementIpData::GetNumData(const StructureBase& rStructure) const
{
    unsigned int numComponents;

    switch (mIpDataType)
    {
    case IpData::eIpStaticDataType::ENGINEERING_STRESS:
    case IpData::eIpStaticDataType::ENGINEERING_STRAIN:
    case IpData::eIpStaticDataType::SHRINKAGE_STRAIN:
        numComponents = 6;
        break;
    case IpData::eIpStaticDataType::DAMAGE:
        numComponents = 1;
        break;
    case IpData::eIpStaticDataType::SLIP:
    case IpData::eIpStaticDataType::BOND_STRESS:
        numComponents = 3;
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Ip data type not supported yet.");
    } // switch

    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));
    return element->GetNumIntegrationPoints() * numComponents;
}
