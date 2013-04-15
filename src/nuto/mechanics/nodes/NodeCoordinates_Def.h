// $Id$
#ifndef NODECOORDINATES_DEF_H
#define NODECOORDINATES_DEF_H

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
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for all nodes
template <int TNumCoordinates>
class NodeCoordinates: public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeCoordinates();


    //! @brief destructor
    ~NodeCoordinates();

    //! @brief assignment operator
    void operator=(NodeCoordinates const& rOther);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp("NodeCoordinates",boost::serialization::base_object< NodeBase >(*this));
    	ar & BOOST_SERIALIZATION_NVP(mCoordinates);
    }
#endif // ENABLE_SERIALIZATION
    void GetCoordinates1D(double rCoordinates[1] )const override;

    void GetCoordinates2D(double rCoordinates[2] )const override;

    void GetCoordinates3D(double rCoordinates[3] )const override;

    void SetCoordinates1D(const double rCoordinates[1] ) override;

    void SetCoordinates2D(const double rCoordinates[2] ) override;

    void SetCoordinates3D(const double rCoordinates[3] ) override;

    double GetCoordinate(short rComponent)const override;

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates (=dimension)
    int GetNumCoordinates()const override;

    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    std::string GetNodeTypeStr()const override;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    //Node::eNodeType GetNodeType()const;

#ifdef ENABLE_VISUALIZE
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;
#endif // ENABLE_VISUALIZE

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeCoordinates* Clone()const;

protected:
    boost::array<double, TNumCoordinates> mCoordinates;

};

}//namespace NuTo


#endif //NODEDOF_DEF_H

