// $Id: $
#ifndef NodeTemperaturesBase_H
#define NodeTemperaturesBase_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having rotational degrees of freedom
class NodeTemperaturesBase : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeTemperaturesBase() : NodeBase ()
    {
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {}
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of temperatures of the node
    //! @return number of temperatures
    virtual int GetNumTemperatures()const=0;

    //! @brief gives the global DOF of a temperature component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofTemperature(int rComponent)const=0;

    //! @brief set the rotations
    //! @param rRotations  given rotations
    virtual void SetTemperatures(const double *rTemperatures)=0;

    //! @brief writes the temperature of a node to the prescribed pointer
    //! @param rTemperatur temperature
    virtual void GetTemperatur(double *rTemperatures)const=0;

protected:
};

}

#endif //NodeTemperaturesBase_H
