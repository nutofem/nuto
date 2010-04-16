#ifndef NODE_DISPLACEMENTSBASE_H
#define NODE_DISPLACEMENTSBASE_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having displacement degrees of freedom
class NodeDisplacementsBase : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeDisplacementsBase() : NodeBase ()
    {
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of displacements of the node
    //! @return number of displacements
    virtual int GetNumDisplacements()const=0;

    //! @brief gives the global DOF of a displacement component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofDisplacement(int rComponent)const=0;

    //! @brief set the displacements
    //! @param rDisplacements  given displacements
    virtual void SetDisplacements(const double *rDisplacements)=0;

    //! @brief writes the displacements of a node to the prescribed pointer
    //! @param rDisplacements displacements
    virtual void GetDisplacements(double *rDisplacements)const=0;

protected:
};
}//namespace NuTo
#endif //NODE_DISPLACEMENTSBASE_H
