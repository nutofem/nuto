// $Id: $
#ifndef NodeRadius_H
#define NodeRadius_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having a radius (required for particle models)
class NodeRadius : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeRadius();

    //! @brief constructor
    NodeRadius (double rRadius);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the number of radii
    //! @return number of radii (1)
    int GetNumRadius()const;

    //! @brief returns the radius of the node
    //! @param rRadius ... radius
    void GetRadius(double rRadius[1])const;

    //! @brief set the radius
    //! @param rRadius  given radius
    void SetRadius(const double rRadius[1]);

protected:
    double mRadius;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeRadius)
#endif // ENABLE_SERIALIZATION
#endif //NodeRadius_H
