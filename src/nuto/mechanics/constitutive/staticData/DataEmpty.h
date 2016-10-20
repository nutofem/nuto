//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once

namespace NuTo
{
class SerializeStreamIn;
class SerializeStreamOut;

namespace Constitutive
{
namespace StaticData
{
//! @brief empty static data class.
class DataEmpty
{
public:
    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    void NuToSerializeSave(SerializeStreamOut& rStream) {}

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) {}
};
} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo
