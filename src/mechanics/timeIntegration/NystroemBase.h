
#pragma once

#include "mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for explicit time integration using Nyström methods
//! see for example Hairer, Lubich and Wanner ("Geometrical numerical integration")
class NystroemBase : public TimeIntegrationBase
{

public:
    //! @brief constructor
    NystroemBase(StructureBase* rStructure);

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    void Solve(double rTimeDelta) override;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

    //! @brief ... return number of intermediate stages
    virtual int GetNumStages() const = 0;

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    virtual double GetStageTimeFactor(int rStage) const = 0;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    virtual void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const = 0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    virtual double GetStageWeights1(int rStage) const = 0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    //! for standard Runge Kutta-Methods, this is identical
    //! for Nyström methods, this migh be different
    virtual double GetStageWeights2(int rStage) const = 0;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    virtual bool HasTimeChanged(int rStage) const = 0;

    //! @brief ... set whether to use a diagonal mass matrix (standard) or a full mass matrix
    void SetUseLumpedMass(bool rUseDiagonalMassMatrix)
    {
        mUseDiagonalMassMatrix = rUseDiagonalMassMatrix;
    }


protected:
    double mTime = 0.;
private:
    // use diagonal mass matrix (standard is true, only for test cases use false)
    bool mUseDiagonalMassMatrix;
};
} // namespace NuTo
