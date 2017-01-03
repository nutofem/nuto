
#include <boost/assign/ptr_map_inserter.hpp>
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/timeIntegration/ResultElementGroupIpData.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"

NuTo::ResultElementGroupIpData::ResultElementGroupIpData(const std::string& rIdent, int rElementGroupId, int rComponent, const std::vector<int> &rIP, NuTo::IpData::eIpStaticDataType rIpDataType) :
        ResultBase(rIdent),
        mElementGroupId(rElementGroupId),
        mComponent(rComponent),
        mIP(rIP),
        mIpDataType(rIpDataType)
{}

void NuTo::ResultElementGroupIpData::CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot)
{
    assert(rTimeStepPlot>=0);

    FullVector<int, Eigen::Dynamic> goupMemberIds = rStructure.GroupGetMemberIds(mElementGroupId);

    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ipValues(goupMemberIds.rows()*mIP.size(),4);
    this->CalculateValues(rStructure,ipValues);

    if (rTimeStepPlot*goupMemberIds.rows()*mIP.size() >= (size_t)mData.GetNumRows())
    {
        this->Resize(rStructure, 2*(rTimeStepPlot+1),false);
    }
    if (ipValues.GetNumColumns()!=mData.GetNumColumns())
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) +"\t: The allocated number of columns is wrong.");
    mData.SetBlock(rTimeStepPlot*goupMemberIds.rows()*mIP.size(), 0, ipValues);
}

void NuTo::ResultElementGroupIpData::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rValues) const
{
    FullMatrix<int, Eigen::Dynamic> goupMemberIds = rStructure.GroupGetMemberIds(mElementGroupId);

    int count = 0;
    for(int i = 0; i < goupMemberIds.rows(); i++)
    {
        const ElementBase* element(rStructure.ElementGetElementPtr(goupMemberIds.GetValue(i)));

        std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
        elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(mIpDataType);

        const_cast<ElementBase*>(element)->Evaluate(elementOutput);

        const auto& ipDataResult = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap()[mIpDataType];

        for(auto &it : mIP)
        {
            assert(mComponent < 0 || ipDataResult.GetNumRows() > mComponent || ipDataResult.GetNumColumns() > it);

            rValues.SetBlock(count, 0, element->GetGlobalIntegrationPointCoordinates(it).transpose());
            rValues(count, 3) = ipDataResult(mComponent, it);
            count++;
        }
    }
}

void NuTo::ResultElementGroupIpData::GetNumData(const StructureBase& rStructure, int &rows, int &cols) const
{
    FullVector<int, Eigen::Dynamic> goupMemberIds = rStructure.GroupGetMemberIds(mElementGroupId);

    rows = goupMemberIds.rows()*mIP.size();
    cols = 4 ; // 3 for the x,y, and z coordinates and 1 for the ip value
}


NuTo::eTimeIntegrationResultType NuTo::ResultElementGroupIpData::GetResultType() const
{
    switch (mIpDataType)
    {
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS:
        return NuTo::eTimeIntegrationResultType::ELEMENTGROUP_IP_STRESS;
    case NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN:
        return NuTo::eTimeIntegrationResultType::ELEMENTGROUP_IP_STRAIN;
    case NuTo::IpData::eIpStaticDataType::DAMAGE:
        return NuTo::eTimeIntegrationResultType::ELEMENTGROUP_IP_DAMAGE;
    case NuTo::IpData::eIpStaticDataType::BOND_STRESS:
        return NuTo::eTimeIntegrationResultType::ELEMENTGROUP_IP_BOND_STRESS;
    case NuTo::IpData::eIpStaticDataType::SLIP:
        return NuTo::eTimeIntegrationResultType::ELEMENTGROUP_IP_SLIP;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Ip data type not supported yet.");
    } // switch

}

void NuTo::ResultElementGroupIpData::Info() const
{
    std::cout << "ResultElementGroupIpData Info:      " << std::endl;
    std::cout << "Integration point data type:   " << IpData::IpStaticDataTypeToString(mIpDataType) << std::endl;
    std::cout << "ElementGroup id:                    " << mElementGroupId << std::endl;
}


