#include "boost/filesystem.hpp"

#include "base/Timer.h"

#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/timeIntegration/postProcessing/ResultElementIpData.h"
#include "mechanics/timeIntegration/postProcessing/ResultElementIpData.h"
#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeForce.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeDisp.h"
#include "mechanics/timeIntegration/postProcessing/ResultNodeAcceleration.h"
#include "mechanics/timeIntegration/postProcessing/ResultTime.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/nodes/NodeEnum.h"

#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/structures/Assembler.h"

using namespace NuTo;

NuTo::TimeIntegrationBase::TimeIntegrationBase(StructureBase* rStructure)
    : mStructure(rStructure)
    , mSolver(std::make_unique<SolverMUMPS>(false))
    , mLoadVectorStatic(rStructure->GetDofStatus())
    , mLoadVectorTimeDependent(rStructure->GetDofStatus())
    , mToleranceResidual(rStructure->GetDofStatus())
    , mCallback(nullptr)
    , mPostProcessor(std::make_unique<PostProcessor>(*rStructure, mTimeControl))
{
    ResetForNextLoad();
}

NuTo::TimeIntegrationBase::~TimeIntegrationBase()
{
}

void NuTo::TimeIntegrationBase::ResetForNextLoad()
{
    mTimeDependentLoadCase = -1;
    mTimeDependentLoadFactor.resize(0, 0);
}

const NuTo::BlockScalar& NuTo::TimeIntegrationBase::GetToleranceResidual() const
{
    return mToleranceResidual;
}

void NuTo::TimeIntegrationBase::UpdateConstraints(double rCurrentTime)
{
    mStructure->GetAssembler().ConstraintUpdateRhs(rCurrentTime);
}

void NuTo::TimeIntegrationBase::SetTimeDependentLoadCase(int rTimeDependentLoadCase,
                                                         const Eigen::MatrixXd& rTimeDependentLoadFactor)
{
    if (rTimeDependentLoadFactor.cols() != 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of columns must be 2, first column contains the time, "
                                                      "second column contains the corresponding value.");
    if (rTimeDependentLoadFactor.rows() < 2)
        throw MechanicsException(__PRETTY_FUNCTION__, "number of rows must be at least 2.");
    if (rTimeDependentLoadFactor(0, 0) != 0)
        throw MechanicsException(__PRETTY_FUNCTION__, "the first time should always be zero.");

    // check, if the time is monotonically increasing
    for (int count = 0; count < rTimeDependentLoadFactor.rows() - 1; count++)
    {
        if (rTimeDependentLoadFactor(count, 0) >= rTimeDependentLoadFactor(count + 1, 0))
            throw MechanicsException(__PRETTY_FUNCTION__, "time has to increase monotonically.");
    }

    mTimeDependentLoadFactor = rTimeDependentLoadFactor;
    mTimeDependentLoadCase = rTimeDependentLoadCase;
}

void NuTo::TimeIntegrationBase::SetToleranceResidual(NuTo::Node::eDof rDof, double rTolerance)
{
    mToleranceResidual[rDof] = rTolerance;
}




void NuTo::TimeIntegrationBase::CalculateStaticAndTimeDependentExternalLoad()
{
    mLoadVectorStatic = StructureOutputBlockVector(mStructure->GetDofStatus(), true);
    mLoadVectorTimeDependent = StructureOutputBlockVector(mStructure->GetDofStatus(), true);

    auto tmp = mStructure->BuildGlobalExternalLoadVector();
    if (mTimeDependentLoadCase == 0)
    {
        mLoadVectorTimeDependent += tmp;
    }
    else
    {
        mLoadVectorStatic += tmp;
    }
    mStructure->GetLogger() << "Sum of loads is " << tmp.J.Export().colwise().sum() + tmp.K.Export().colwise().sum()
                            << "\n";
}




NuTo::StructureOutputBlockVector NuTo::TimeIntegrationBase::CalculateCurrentExternalLoad(double curTime)
{
    if (mTimeDependentLoadCase != -1)
    {
        if (mTimeDependentLoadFactor.rows() == 0)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "TimeDependentLoadFactor not set.");
        }
        int curStep(0);
        while (mTimeDependentLoadFactor(curStep, 0) < curTime && curStep < mTimeDependentLoadFactor.rows() - 1)
            curStep++;
        if (curStep == 0)
            curStep++;

        // extract the two data points
        double s1 = mTimeDependentLoadFactor(curStep - 1, 1);
        double s2 = mTimeDependentLoadFactor(curStep, 1);
        double t1 = mTimeDependentLoadFactor(curStep - 1, 0);
        double t2 = mTimeDependentLoadFactor(curStep, 0);

        double s = s1 + (s2 - s1) / (t2 - t1) * (curTime - t1);

        return mLoadVectorStatic + mLoadVectorTimeDependent * s;
    }
    else
    {
        return mLoadVectorStatic;
    }
}

const NuTo::BlockFullVector<double>& NuTo::TimeIntegrationBase::UpdateAndGetConstraintRHS(double rCurrentTime)
{
    UpdateConstraints(rCurrentTime);
    return mStructure->GetAssembler().GetConstraintRhs();
}

const NuTo::BlockFullVector<double>&
NuTo::TimeIntegrationBase::UpdateAndGetAndMergeConstraintRHS(double rCurrentTime, StructureOutputBlockVector& rDof_dt0)
{
    UpdateConstraints(rCurrentTime);


    rDof_dt0.K = mStructure->NodeCalculateDependentDofValues(rDof_dt0.J);
    mStructure->NodeMergeDofValues(0, rDof_dt0);

    mStructure->ElementTotalUpdateTmpStaticData();
    return mStructure->GetAssembler().GetConstraintRhs();
}


std::array<StructureOutputBlockVector, 3> NuTo::TimeIntegrationBase::ExtractDofValues() const
{
    const auto& dofStatus = mStructure->GetDofStatus();
    std::array<StructureOutputBlockVector, 3> dofValues = {StructureOutputBlockVector(dofStatus, true),
                                                           StructureOutputBlockVector(dofStatus, true),
                                                           StructureOutputBlockVector(dofStatus, true)};
    dofValues[0] = mStructure->NodeExtractDofValues(0);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        dofValues[1] = mStructure->NodeExtractDofValues(1);

    if (mStructure->GetNumTimeDerivatives() >= 2)
        dofValues[2] = mStructure->NodeExtractDofValues(2);
    return dofValues;
}

double NuTo::TimeIntegrationBase::CalculateNorm(const BlockFullVector<double>& rResidual) const
{
    double norm = 0;
    for (auto rDofType : rResidual.GetDofStatus().GetActiveDofTypes())
        norm += rResidual[rDofType].norm();

    return norm;
}


void NuTo::TimeIntegrationBase::AddCalculationStep(const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs.push_back(rActiveDofs);
}


void NuTo::TimeIntegrationBase::SetNumCalculationSteps(int rNumSteps)
{
    mStepActiveDofs.resize(rNumSteps);
}

void NuTo::TimeIntegrationBase::SetActiveDofsCalculationStep(int rStepNum,
                                                             const std::set<NuTo::Node::eDof>& rActiveDofs)
{
    mStepActiveDofs[rStepNum] = rActiveDofs;
}

void NuTo::TimeIntegrationBase::Info() const
{
}


bool TimeIntegrationBase::GetShowTime() const
{
    return mShowTime;
}

void TimeIntegrationBase::SetShowTime(bool showTime)
{
    mShowTime = showTime;
}
