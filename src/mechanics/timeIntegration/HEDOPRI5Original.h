

#include "mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
class HEDOPRI5Original : public TimeIntegrationBase
{

public:
    //! @brief constructor
    HEDOPRI5Original(StructureBase* rStructure);

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep() const
    {
        return true;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep() const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in
    //! the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId() const;

    //! @brief ... return number of intermediate stages
    int GetNumStages() const
    {
        return 8;
    }

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    double GetStageTimeFactor(int rStage) const;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    double GetStageWeights(int rStage) const;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    bool HasTimeChanged(int rStage) const;

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta);

protected:
    // empty private construct required for serialization
    HEDOPRI5Original()
    {
    }
};
} // namespace NuTo
