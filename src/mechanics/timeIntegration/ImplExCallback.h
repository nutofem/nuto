//
// Created by Thomas Titscher on 1/17/17.
//

#pragma once

namespace NuTo
{

class ImplExCallback
{
public:

    virtual ~ImplExCallback() = default;

    //! @brief is called after an ImplEx step
    //! @param rMaxError maximum extrapolation error [from GetError(...)]
    //! @param rThreshold error threshold (user defined number)
    //! @return true: continue with next time step [from GetNewTimeStep(...)], false: revert time step and start again
    virtual bool AcceptSolution(double rMaxError, double rThreshold);

    //! @brief caclulates the length of the next time step based on the extrapolation error
    //! @param rMaxError maximum extrapolation error [from GetError(...)]
    //! @param rThreshold error threshold (user defined number)
    //! @param mOldTimeStep old time step
    //! @return length of the next time step
    virtual double GetNewTimeStep(double rMaxError, double rThreshold, double mOldTimeStep);

    //! @brief calculates the error
    //! @param k_tilde_n implex extrapolation at n
    //! @param k_n implicit value at n
    //! @param k_nn implicit value at n-1
    //! @return
    virtual double GetError(double k_tilde_n, double k_n, double k_nn) const;

private:
    bool mForceAcceptOfNextSolution;
};

} // namespace NuTo