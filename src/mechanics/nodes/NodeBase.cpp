#include <iostream>
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "visualize/VisualizeComponent.h"
#include "visualize/VisualizeEnum.h"
#include "visualize/VisualizeUnstructuredGrid.h"
#include <cassert>

using namespace NuTo;

#ifdef ENABLE_VISUALIZE

Eigen::Vector3d NodeBase::GetPointVectorData(Node::eDof rDofType, int rTimeDerivative) const
{
    int dim = this->GetNum(rDofType);
    assert(dim != 0);
    Eigen::Vector3d data = Eigen::Vector3d::Zero();
    data.block(0, 0, dim, 1) = this->Get(rDofType, rTimeDerivative);
    return data;
}

void NodeBase::Visualize(VisualizeUnstructuredGrid& rVisualize,
                         const std::list<std::shared_ptr<VisualizeComponent>>& rVisualizationList) const
{
    Eigen::Matrix<double, 3, 1> coordinates = Eigen::Matrix<double, 3, 1>::Zero();
    int dim = this->GetNum(Node::eDof::COORDINATES);
    coordinates.block(0, 0, dim, 1) = this->Get(Node::eDof::COORDINATES);

    unsigned int PointId = rVisualize.AddPoint(coordinates.data());
    //	std::cout << "add point " << PointId << std::endl;


    // store data
    for (auto const& it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case eVisualizeWhat::DISPLACEMENTS:
        {
            if (GetNum(Node::eDof::DISPLACEMENTS) == 0)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::DISPLACEMENTS, 0).data());
        }
        break;
        case eVisualizeWhat::VELOCITY:
        {
            if (GetNum(Node::eDof::DISPLACEMENTS) == 0 && GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < 1)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::DISPLACEMENTS, 1).data());
        }
        break;
        case eVisualizeWhat::ACCELERATION:
        {
            if (GetNum(Node::eDof::DISPLACEMENTS) == 0 && GetNumTimeDerivatives(Node::eDof::DISPLACEMENTS) < 2)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::DISPLACEMENTS, 2).data());
        }
        break;
        case eVisualizeWhat::ROTATION:
        {
            if (GetNum(Node::eDof::ROTATIONS) == 0)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::ROTATIONS, 0).data());
        }
        break;
        case eVisualizeWhat::ANGULAR_VELOCITY:
        {
            if (GetNum(Node::eDof::ROTATIONS) == 0 && GetNumTimeDerivatives(Node::eDof::ROTATIONS) < 1)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::ROTATIONS, 1).data());
        }
        break;
        case eVisualizeWhat::ANGULAR_ACCELERATION:
        {
            if (GetNum(Node::eDof::ROTATIONS) == 0 && GetNumTimeDerivatives(Node::eDof::ROTATIONS) < 2)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(),
                                          GetPointVectorData(Node::eDof::ROTATIONS, 2).data());
        }
        break;
        default:
            break;
        }
    }
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NodeBase::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NodeBase::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase"
              << "\n";
#endif
// nothing to do here, no members...
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NodeBase)
#endif // ENABLE_SERIALIZATION

int NodeBase::GetDof(Node::eDof rDof) const
{
    assert(GetNum(rDof) == 1);
    return GetDof(rDof, 0);
}

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const NodeBase& node)
{
    node.Info(out);
    return out;
}
} /* NuTo */

