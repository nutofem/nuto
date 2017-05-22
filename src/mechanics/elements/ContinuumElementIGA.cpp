#include "mechanics/elements/ContinuumElementIGA.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

#include "mechanics/sections/SectionTruss.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/elements/ElementOutputBase.h"
#include "mechanics/elements/ElementOutputIpData.h"

#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"

#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/EvaluateDataContinuum.h"
#include "mechanics/elements/IpDataEnum.h"

#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

template <int TDim>
NuTo::ContinuumElementIGA<TDim>::ContinuumElementIGA(const std::vector<NuTo::NodeBase*>& rNodes,
                                                     const Eigen::MatrixXd& rKnots, const Eigen::VectorXi& rKnotIDs,
                                                     const InterpolationType& rInterpolationType,
                                                     const DofStatus& dofStatus)
    : ContinuumElement<TDim>(rNodes, rInterpolationType, dofStatus)
    , mKnots(rKnots)
    , mKnotIDs(rKnotIDs)
{
}

template<int TDim>
Eigen::Matrix<double, TDim, TDim> NuTo::ContinuumElementIGA<TDim>::CalculateJacobianParametricSpaceIGA() const
{
    Eigen::Matrix<double, TDim, TDim> jac;
    jac.setZero(TDim, TDim);
    for(int i = 0; i < TDim; i++) jac(i,i) = 0.5*(mKnots(i,1) - mKnots(i,0));

    return jac;
}

template<int TDim>
const Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::GetIntegrationPointVolume() const
{
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::VectorXd volume(this->GetNumIntegrationPoints());
    for (int theIP = 0; theIP < this->GetNumIntegrationPoints(); theIP++)
    {
        const auto ipCoords = this->GetIntegrationType().GetLocalIntegrationPointCoordinates(theIP);
        Eigen::MatrixXd derivativeShapeFunctionsNatural = this->mInterpolationType->Get(Node::eDof::COORDINATES).DerivativeShapeFunctionsNaturalIGA(ipCoords, mKnotIDs);
        double detJacobian = this->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates).determinant();
        volume[theIP] = detJacobian * this->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::CheckElement()
{
    int numIntegrationPoints = this->GetNumIntegrationPoints();

    if (numIntegrationPoints < 1)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] invalid integration type.");
    }

    int theIP = 0;

    const auto ipCoords = this->GetIntegrationType().GetLocalIntegrationPointCoordinates(theIP);
    Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::eDof::COORDINATES).DerivativeShapeFunctionsNaturalIGA(ipCoords, mKnotIDs);
    Eigen::MatrixXd nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);
    double detJacobian = this->CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);
    }

    double size = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const auto ipCoords = this->GetIntegrationType().GetLocalIntegrationPointCoordinates(theIP);
        Eigen::MatrixXd derivativeShapeFunctions = this->mInterpolationType->Get(Node::eDof::COORDINATES).DerivativeShapeFunctionsNaturalIGA(ipCoords, mKnotIDs);
        detJacobian = this->CalculateJacobian(derivativeShapeFunctions, nodeCoordinates).determinant();
        if (detJacobian <= 0)
        {
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian <= zero, no inversion possible.");
        }
        size += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    assert(std::abs(GetIntegrationPointVolume().sum() / size - 1) < 1.e-6); // just to be sure ...

    // check element length
    if (size < 1e-14)
    {
        MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] element with zero size (check nodes).");
    }

}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuum<TDim> &rData, int rTheIP) const
{
    // calculate Jacobian
    const auto ipCoords = this->GetIntegrationType().GetLocalIntegrationPointCoordinates(rTheIP);
    const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = this->mInterpolationType->Get(Node::eDof::COORDINATES).DerivativeShapeFunctionsNaturalIGA(ipCoords, mKnotIDs);

    Eigen::Matrix<double, TDim, TDim> jacobian = this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, rData.mNodalValues[Node::eDof::COORDINATES]);

    // there are two mappings in IGA: 1) reference element <=> parametric space (knots) 2) parametric space <=> physical space
    // the B-matrix only the invJac of mapping 2) is needed
    rData.mDetJacobian = jacobian.determinant();
    for(int i = 0; i < TDim; i++) rData.mDetJacobian *= 0.5*(mKnots(i,1) - mKnots(i,0));

    if (rData.mDetJacobian == 0)
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Determinant of the Jacobian is zero, no inversion possible.");

    Eigen::Matrix<double, TDim, TDim> invJacobian = jacobian.inverse();

    // calculate shape functions and their derivatives
    for (auto dof : this->mInterpolationType->GetDofs())
    {
        if (dof == Node::eDof::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = this->mInterpolationType->Get(dof);
        rData.mNIGA[dof] = interpolationType.MatrixNIGA(ipCoords, mKnotIDs);

        rData.mB[dof] = this->CalculateMatrixB(dof, interpolationType.DerivativeShapeFunctionsNaturalIGA(ipCoords, mKnotIDs), invJacobian);
    }
}

template<int TDim>
NuTo::NodeBase* NuTo::ContinuumElementIGA<TDim>::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
const NuTo::NodeBase* NuTo::ContinuumElementIGA<TDim>::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    return (this->mNodes)[rLocalNodeNumber];
}

template<int TDim>
void NuTo::ContinuumElementIGA<TDim>::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )((this->mNodes).size()));
    assert(rNode != nullptr);
    (this->mNodes)[rLocalNodeNumber] = rNode;
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationType = this->mInterpolationType->Get(rDofType);
    Eigen::MatrixXd nodalValues = this->ExtractNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd matrixN = interpolationType.MatrixNIGA(rNaturalCoordinates, mKnotIDs);

//    std::cout << matrixN * nodalValues << std::endl;

    return matrixN * nodalValues;
}

template<int TDim>
Eigen::VectorXd NuTo::ContinuumElementIGA<TDim>::InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const
{
//    const InterpolationBase& interpolationTypeCoords = this->GetInterpolationType()->Get(Node::eDof::COORDINATES);
//    Eigen::VectorXd parameter = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(rNaturalCoordinates, rSurfaceId, mKnots);

    Eigen::VectorXd nodalInitial       = this->ExtractNodeValues(rTimeDerivative, Node::eDof::COORDINATES);
    Eigen::VectorXd nodalDisplacements = this->ExtractNodeValues(rTimeDerivative, Node::eDof::DISPLACEMENTS);

    Eigen::MatrixXd matrixNDerivative = this->mInterpolationType->Get(Node::eDof::COORDINATES).MatrixNDerivativeIGA(rParameter, mKnotIDs, rDerivative, rDirection);

    return matrixNDerivative * (nodalInitial + nodalDisplacements);
}

namespace NuTo // template specialization in *.cpp somehow requires the definition to be in the namespace...
{

template<>
double NuTo::ContinuumElementIGA<1>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    const auto ipCoords = this->GetIntegrationType().GetLocalIntegrationPointCoordinates(rTheIP);
    const auto globalIPCoordinate = this->InterpolateDofGlobal(0, ipCoords, Node::eDof::COORDINATES);
    return rDetJacobian * GetIntegrationPointWeight(rTheIP) * mSection->GetArea(globalIPCoordinate(0, 0));
}

template<>
double NuTo::ContinuumElementIGA<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return rDetJacobian * GetIntegrationPointWeight(rTheIP) * mSection->GetThickness();
}

}  // namespace NuTo

template class NuTo::ContinuumElementIGA<1>;
template class NuTo::ContinuumElementIGA<2>;


#ifdef ENABLE_SERIALIZATION
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElement",boost::serialization::base_object<ContinuumElement<TDim> >(*this));
    ar & boost::serialization::make_nvp("mKnots_size", mKnots);
    ar & boost::serialization::make_nvp("mKnotIDs_size", mKnotIDs);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ContinuumElementIGA" << std::endl;
#endif
}

template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ContinuumElementIGA<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive>
void NuTo::ContinuumElementIGA<TDim>::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize ContinuumElementIGA " << std::endl;
#endif
    ar & boost::serialization::make_nvp("ContinuumElementIGA_ElementBase",boost::serialization::base_object<ContinuumElement<TDim> >(*this));
    ar & BOOST_SERIALIZATION_NVP(mKnots);
    ar & BOOST_SERIALIZATION_NVP(mKnotIDs);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize ContinuumElementIGA" << std::endl;
#endif
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<1>)), "ContinuumElementIGA_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ContinuumElementIGA<2>)), "ContinuumElementIGA_2")
#endif // ENABLE_SERIALIZATION
