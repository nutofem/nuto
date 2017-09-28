#pragma once

#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "mechanics/groups/Group.h"

namespace NuTo
{
template <int TDim>
class ContinuumContactElement : public ContinuumBoundaryElement<TDim>
{
public:
    ContinuumContactElement(const ContinuumElement<TDim>& rSlaveElement, int rSurfaceId,
                            const Group<ElementBase>* elementGroup, const Group<NodeBase>* nodeGroup,
                            const IntegrationTypeBase& integrationType);

    virtual ~ContinuumContactElement() = default;

    void ProjectIntegrationPointOnMaster();

    void CalculateElementOutputs(std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
                                 EvaluateDataContinuumBoundary<TDim>& rData, int rTheIP,
                                 const ConstitutiveInputMap& constitutiveInput,
                                 const ConstitutiveOutputMap& constitutiveOutput) const;

    void CalculateElementOutputGapMatrixMortar(BlockFullMatrix<double>& rGapMatrix,
                                               EvaluateDataContinuumBoundary<TDim>& rData,
                                               const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const;

protected:
    std::vector<std::pair<const ContinuumElement<TDim>*, int>> mElementsMaster;

    const IntegrationTypeBase* mIntegrationType;
};
} /* namespace NuTo */
