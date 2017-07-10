#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/NewmarkBase.h"

#include "mechanics/structures/StructureOutputBlockVector.h"

NuTo::NewmarkBase::NewmarkBase(StructureBase* rStructure)
    : TimeIntegrationBase(rStructure), mDampingMatrix(rStructure->GetDofStatus())
{
}

void NuTo::NewmarkBase::MergeDofValues(const StructureOutputBlockVector& rDof_dt0,
                                       const StructureOutputBlockVector& rDof_dt1,
                                       const StructureOutputBlockVector& rDof_dt2, bool rMergeAll)
{
    mStructure->NodeMergeDofValues(0, rDof_dt0.J, rDof_dt0.K);

    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        if (rMergeAll or mMuDampingMass == 0)
        {
            mStructure->NodeMergeDofValues(1, rDof_dt1.J, rDof_dt1.K);
        }
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        if (rMergeAll)
        {
            mStructure->NodeMergeDofValues(2, rDof_dt2.J, rDof_dt2.K);
        }
    }
    mStructure->ElementTotalUpdateTmpStaticData();
}
