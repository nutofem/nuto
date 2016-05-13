
#include <boost/assign/ptr_map_inserter.hpp>
#include "nuto/mechanics/timeIntegration/ResultElementIpData.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"

#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ResultElementIpData::ResultElementIpData(const std::string& rIdent, int rElementId, NuTo::IpData::eIpStaticDataType rIpDataType) :
        ResultBase(rIdent),
        mElementId(rElementId),
        mIpDataType(rIpDataType)
{
}

void NuTo::ResultElementIpData::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
    assert(rTimeStepPlot>=0);
    FullMatrix<double,1,Eigen::Dynamic> ipValues(1,this->GetNumData(rStructure));
    this->CalculateValues(rStructure,ipValues);
    if (rTimeStepPlot>=mData.GetNumRows())
    {
        this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
    }
    if (ipValues.GetNumColumns()!=mData.GetNumColumns())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: The allocated number of columns is wrong.");
    mData.SetRow(rTimeStepPlot,ipValues);
}

void NuTo::ResultElementIpData::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const
{
    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));

    std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(mIpDataType);

    const_cast<ElementBase*>(element)->Evaluate(elementOutput);

    const auto& ipDataResult = elementOutput.at(Element::IP_DATA)->GetIpData().GetIpDataMap()[mIpDataType];

    // iterate over all ips
    assert(ipDataResult.GetNumColumns() == element->GetNumIntegrationPoints());
    unsigned int numComponents = ipDataResult.GetNumRows();
    for (int iCol = 0; iCol < ipDataResult.GetNumColumns(); ++iCol)
    {
        rValues.SetBlock(0, iCol * numComponents, ipDataResult.GetColumn(iCol).Trans());
    }
}

int NuTo::ResultElementIpData::GetNumData(const StructureBase& rStructure) const
{
    unsigned int numComponents;

    switch (mIpDataType)
    {
    case NuTo::IpData::ENGINEERING_STRESS:
    case NuTo::IpData::ENGINEERING_STRAIN:
    case NuTo::IpData::SHRINKAGE_STRAIN:
        numComponents = 6;
        break;
    case NuTo::IpData::DAMAGE:
        numComponents = 1;
        break;
    case NuTo::IpData::SLIP:
    case NuTo::IpData::BOND_STRESS:
        numComponents = 3;
            break;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Ip data type not supported yet.");
    } // switch

    const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));
    return element->GetNumIntegrationPoints() * numComponents;
}


NuTo::TimeIntegration::eResultType NuTo::ResultElementIpData::GetResultType() const
{
    switch (mIpDataType)
    {
    case NuTo::IpData::ENGINEERING_STRESS:
        return NuTo::TimeIntegration::ELEMENT_IP_STRESS;
    case NuTo::IpData::ENGINEERING_STRAIN:
        return NuTo::TimeIntegration::ELEMENT_IP_STRAIN;
    case NuTo::IpData::DAMAGE:
        return NuTo::TimeIntegration::ELEMENT_IP_DAMAGE;
    case NuTo::IpData::BOND_STRESS:
        return NuTo::TimeIntegration::ELEMENT_IP_BOND_STRESS;
    case NuTo::IpData::SLIP:
        return NuTo::TimeIntegration::ELEMENT_IP_SLIP;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Ip data type not supported yet.");
    } // switch

}


