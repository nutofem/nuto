#pragma once

#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

namespace NuTo
{
class SerializeStreamIn;
class SerializeStreamOut;
class MisesPlasticityEngineeringStress;

namespace Constitutive
{
namespace StaticData
{
//! @brief Storing Mises plasticity static data.
template <int TDim>
class DataMisesPlasticity
{
    friend class MisesPlasticityEngineeringStress;
public:
    DataMisesPlasticity();

    //! @brief Get the equivalent plastic strain.
    double GetEquivalentPlasticStrain() const;

    //! @brief Get the plastic strain.
    EngineeringStrain<TDim> GetPlasticStrain() const;

    //! @brief Get the back stress.
    EngineeringStress<TDim> GetBackStress() const;

    //! @brief Set the equivalent plastic strain.
    void SetEquivalentPlasticStrain(double newEpsilonPEq);

    //! @brief Set the plastic strain.
    void SetPlasticStrain(EngineeringStrain<TDim> newEpsilonP);

    //! @brief Set the back stress.
    void SetBackStress(EngineeringStress<TDim> newBackStress);

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream)
    {
        SerializeDataMisesPlasticity(rStream);
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream)
    {
        SerializeDataMisesPlasticity(rStream);
    }

protected:
    //! @brief Accumulated plastic strain (is not always equivalent to epsilon_p)
    double mEpsilonPEq;

    //! @brief plastic strain
    EngineeringStrain<TDim> mEpsilonP;

    //! @brief back stress
    EngineeringStress<TDim> mSigmaB;

private:
    //! @brief defines the serialization of this class
    //! @param rStream serialize input/output stream
    template <typename TStream>
    void SerializeDataMisesPlasticity(TStream &rStream);
};

} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo
