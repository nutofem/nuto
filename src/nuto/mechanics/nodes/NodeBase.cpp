// $Id$
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>

#ifdef ENABLE_VISUALIZE

Eigen::Vector3d NuTo::NodeBase::GetPointVectorData(Node::eDof rDofType, int rTimeDerivative) const
{
    int dim = this->GetNum(rDofType);
    assert(dim != 0);
    Eigen::Vector3d data = Eigen::Vector3d::Zero();
    data.block(0,0,dim, 1) = this->Get(rDofType, rTimeDerivative);
    return data;
}

void NuTo::NodeBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const
{
    Eigen::Matrix<double, 3, 1> coordinates = Eigen::Matrix<double, 3, 1>::Zero();
    int dim = this->GetNum(Node::COORDINATES);
    coordinates.block(0,0,dim, 1) = this->Get(Node::COORDINATES);

    unsigned int PointId = rVisualize.AddPoint(coordinates.data());
    //	std::cout << "add point " << PointId << std::endl;


    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
        {
            if (GetNum(Node::DISPLACEMENTS) == 0)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::DISPLACEMENTS, 0).data());
        }
        break;
        case NuTo::VisualizeBase::VELOCITY:
        {
            if (GetNum(Node::DISPLACEMENTS) == 0 && GetNumTimeDerivatives() < 1)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::DISPLACEMENTS, 1).data());
        }
        break;
        case NuTo::VisualizeBase::ACCELERATION:
        {
            if (GetNum(Node::DISPLACEMENTS) == 0 && GetNumTimeDerivatives() < 2)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::DISPLACEMENTS, 2).data());
        }
        break;
        case NuTo::VisualizeBase::ROTATION:
        {
            if (GetNum(Node::ROTATIONS) == 0)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::ROTATIONS, 0).data());
        }
        break;
        case NuTo::VisualizeBase::ANGULAR_VELOCITY:
        {
            if (GetNum(Node::ROTATIONS) == 0 && GetNumTimeDerivatives() < 1)
                break;
            rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::ROTATIONS, 1).data());
        }
        break;
        case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
        {
            if (GetNum(Node::ROTATIONS) == 0 && GetNumTimeDerivatives() < 2)
                 break;
             rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GetPointVectorData(Node::ROTATIONS, 2).data());
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
template void NuTo::NodeBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase" << "\n";
#endif
    // nothing to do here, no members...
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeBase)
#endif // ENABLE_SERIALIZATION
