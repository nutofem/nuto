#pragma once

#include <eigen3/Eigen/Core>

namespace NuTo
{
class SerializeStreamIn;
class SerializeStreamOut;
namespace Constitutive
{
namespace StaticData
{
//! @brief Storing creep static data.
class DataCreep
{
public:
    // ctor
    DataCreep()
    {
        int a = 1;
    }


    //! @brief sets all the previous data to the value of the current data
    void ProceedToNextTimestep();


    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream)
    {
        SerializeDataCreep(rStream);
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream)
    {
        SerializeDataCreep(rStream);
    }


    //! @brief gets the current time
    //! @return current time
    double GetCurrentTime() const
    {
        return mCurrentTime;
    }

    //! @brief gets the previous time
    //! @return previous time
    double GetPreviousTime() const
    {
        return mPreviousTime;
    }


    //! @brief sets the current time
    //! @param currentTime:  current time
    void SetCurrentTime(double currentTime)
    {
        mCurrentTime = currentTime;
    }

private:
    //! @brief defines the serialization of this class
    //! @param rStream serialize input/output stream
    template <typename TStream>
    void SerializeDataCreep(TStream& rStream);


public:
    Eigen::MatrixXd mHistoryData = Eigen::MatrixXd(0, 0);
    Eigen::VectorXd mHistoryStrain = Eigen::VectorXd(0);
    Eigen::VectorXd mHistoryStress = Eigen::VectorXd(0);
    Eigen::VectorXd mDeltaStress = Eigen::VectorXd(0);
    Eigen::VectorXd mDeltaStrain = Eigen::VectorXd(0);
    Eigen::VectorXd mDeltaCreepStrain = Eigen::VectorXd(0);

protected:
    double mCurrentTime = 0.0;
    double mPreviousTime = 0.0;
};
} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo
