// $Id$
#ifndef NODECOORDINATESDOF_DEF_H
#define NODECOORDINATESDOF_DEF_H

#ifdef ENABLE_SERIALIZATION
// serialize
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/utility/identity_type.hpp>
#include <boost/serialization/array.hpp>
#else
#include <boost/array.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

namespace NuTo
{
template <class T> class FullMatrix;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for all nodes
template <int TNumCoordinates, int TNumTimeDerivatives, int TNumDisplacements, int TNumRotations, int TNumTemperatures>
class NodeCoordinatesDof: public NodeCoordinates<TNumCoordinates>, public NodeDof<TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeCoordinatesDof();


    //! @brief destructor
    ~NodeCoordinatesDof();

    //! @brief assignment operator
    void operator=(NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures> const& rOther);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp("NodeDof",boost::serialization::base_object< NodeBase >(*this));
        ar & boost::serialization::make_nvp("NodeCoordinates",boost::serialization::base_object< NodeBase >(*this));
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr()const;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    //Node::eNodeType GetNodeType()const;

#ifdef ENABLE_VISUALIZE
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;
#endif // ENABLE_VISUALIZE

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeCoordinatesDof<TNumCoordinates,TNumTimeDerivatives,TNumDisplacements,TNumRotations,TNumTemperatures>* Clone()const;

protected:

};

}//namespace NuTo


#endif //NODECOORDINATESDOF_DEF_H

