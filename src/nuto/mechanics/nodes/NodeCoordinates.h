// $Id$

#ifndef NODE_COORDINATES_H
#define NODE_COORDINATES_H
#include "nuto/mechanics/nodes/NodeCoordinates_Def.h"
#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumCoordinates >
NuTo::NodeCoordinates<TNumCoordinates>::NodeCoordinates() : NuTo::NodeBase()
{
}

//! @brief ... destructor
template <int TNumCoordinates >
NuTo::NodeCoordinates<TNumCoordinates>::~NodeCoordinates()
{
}

//! @brief assignment operator
template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::
     operator=(NuTo::NodeCoordinates<TNumCoordinates> const& rOther)
{
    mCoordinates = rOther.mCoordinates;
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::GetCoordinates1D(double rCoordinates[1] )const
{
    assert(TNumCoordinates==1);
    rCoordinates[0] = mCoordinates[0];
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::GetCoordinates2D(double rCoordinates[2] )const
{
    assert(TNumCoordinates==2);
    rCoordinates[0] = mCoordinates[0];
    rCoordinates[1] = mCoordinates[1];
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::GetCoordinates3D(double rCoordinates[3] )const
{
    assert(TNumCoordinates==3);
    rCoordinates[0] = mCoordinates[0];
    rCoordinates[1] = mCoordinates[1];
    rCoordinates[2] = mCoordinates[2];
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::SetCoordinates1D(const double rCoordinates[1] )
{
    assert(TNumCoordinates==1);
    mCoordinates[0] = rCoordinates[0];
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::SetCoordinates2D(const double rCoordinates[2] )
{
    assert(TNumCoordinates==2);
    mCoordinates[0] = rCoordinates[0];
    mCoordinates[1] = rCoordinates[1];
}

template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::SetCoordinates3D(const double rCoordinates[3] )
{
    assert(TNumCoordinates==3);
    mCoordinates[0] = rCoordinates[0];
    mCoordinates[1] = rCoordinates[1];
    mCoordinates[2] = rCoordinates[2];
}

template <int TNumCoordinates >
double NuTo::NodeCoordinates<TNumCoordinates>::GetCoordinate(short rComponent) const
{
    assert(TNumCoordinates>rComponent);
    return mCoordinates[rComponent];
}


//! @brief returns the type of node as a string (all the data stored at the node)
//! @return string
template <int TNumCoordinates >
std::string NuTo::NodeCoordinates<TNumCoordinates>::GetNodeTypeStr()const
{
	throw MechanicsException("");
    return std::string("[NuTo::NodeDof::GetNodeTypeStr]to be done");
}


//! @brief returns the number of coordinates of the node
//! @return number of coordinates (=dimension)
template <int TNumCoordinates >
int NuTo::NodeCoordinates<TNumCoordinates>::GetNumCoordinates()const
{
	return TNumCoordinates;
}

#ifdef ENABLE_VISUALIZE
template <int TNumCoordinates >
void NuTo::NodeCoordinates<TNumCoordinates>::Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{

}
#endif // ENABLE_VISUALIZE

//! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
template <int TNumCoordinates >
NuTo::NodeCoordinates<TNumCoordinates>* NuTo::NodeCoordinates<TNumCoordinates>::
Clone()const
{
    return new NuTo::NodeCoordinates<TNumCoordinates>(*this);
}

#endif //NODE_COORDINATES_H

