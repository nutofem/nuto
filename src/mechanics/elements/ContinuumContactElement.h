#pragma once

#include "mechanics/elements/ContinuumBoundaryElement.h"

namespace NuTo
{
template<int TDim>
class ContinuumContactElement: public ContinuumBoundaryElement<TDim>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    ContinuumContactElement() = default;
#endif // ENABLE_SERIALIZATION
public:
    ContinuumContactElement(const ContinuumElement<TDim>& rSlaveElement, int rSurfaceId, int rElementGroupId, int rNodeGroupId, const IntegrationTypeBase *rIntegrationType);

    virtual ~ContinuumContactElement() = default;

    void ProjectIntegrationPointOnMaster();

    void CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
                                 EvaluateDataContinuumBoundary<TDim> &rData, int rTheIP,
                                 const ConstitutiveInputMap& constitutiveInput,
                                 const ConstitutiveOutputMap& constitutiveOutput) const;

    void CalculateElementOutputGapMatrixMortar(BlockFullMatrix<double>& rGapMatrix,
                                               EvaluateDataContinuumBoundary<TDim> &rData,
                                               const ConstitutiveOutputMap& constitutiveOutput,
                                               int rTheIP) const;

protected:

    std::vector<std::pair<const ContinuumElement<TDim>*, int> > mElementsMaster;

    const IntegrationTypeBase *mIntegrationType;

};
} /* namespace NuTo */

