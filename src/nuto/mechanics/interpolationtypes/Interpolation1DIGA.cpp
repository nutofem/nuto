#include "nuto/mechanics/interpolationtypes/Interpolation1DIGA.h"

NuTo::Interpolation1DIGA::Interpolation1DIGA(NuTo::Node::eDof rDofType,  NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension, int rDegree, const Eigen::VectorXd &rKnots)
    : InterpolationBaseIGA::InterpolationBaseIGA(rDofType, rTypeOrder, rDimension),
      mKnots(rKnots),
      mDegree(rDegree)
{
    Initialize();
}

NuTo::IntegrationType::eIntegrationType NuTo::Interpolation1DIGA::GetStandardIntegrationType() const
{
    switch (mDegree)
    {
    case 0: // (0+0+1)/2 = 0.5 ips
        return NuTo::IntegrationType::IntegrationType1D2NGauss1Ip;
    case 1: // (1+1+1)/2 = 1.5 ips or (1+1+3)/2 = 2.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
    case 2: // (2+2+1)/2 = 2.5 ips or (2+2+3)/2 = 3.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType1D2NGauss3Ip;
    case 3: // (3+3+1)/2 = 3.5 ips or (3+3+3)/2 = 4.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType1D2NGauss4Ip;
    case 4: // (4+4+1)/2 = 4.5 ips or (4+4+3)/2 = 5.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType1D2NGauss5Ip;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + std::to_string(mDegree) + " IGA not implemented");
    }
}

void NuTo::Interpolation1DIGA::UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType)
{
    int numIPs = rIntegrationType.GetNumIntegrationPoints();

    mIPCoordinates.resize(numIPs);

    assert(rIntegrationType.GetCoordinateDimension() == 1);

    for (int iIP = 0; iIP < numIPs; ++iIP)
    {
        double coordinate1D;
        rIntegrationType.GetLocalIntegrationPointCoordinates1D(iIP, coordinate1D);
        mIPCoordinates(iIP) = coordinate1D;
    }

    mUpdateRequired = false;
}


int NuTo::Interpolation1DIGA::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return mDimension;
    case NuTo::Node::DISPLACEMENTS:
        return mDimension;
    case NuTo::Node::TEMPERATURE:
        return 1;
    case NuTo::Node::NONLOCALEQSTRAIN:
        return 1;
    case NuTo::Node::NONLOCALEQPLASTICSTRAIN:
        return 2;
    case NuTo::Node::RELATIVEHUMIDITY:
        return 1;
    case NuTo::Node::WATERVOLUMEFRACTION:
        return 1;
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "dof type not found.");
    }
}


// --- shape functions --- //

Eigen::VectorXd NuTo::Interpolation1DIGA::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rCoordinates(0,0), mDegree, mKnots);
    return ShapeFunctionsIGA::BasisFunctions(rCoordinates(0,0), spanIdx, mDegree, mKnots);
}

Eigen::VectorXd NuTo::Interpolation1DIGA::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates, int rKnotID) const
{
    return ShapeFunctionsIGA::BasisFunctions(rCoordinates(0,0), rKnotID, mDegree, mKnots);
}

Eigen::VectorXd NuTo::Interpolation1DIGA::CalculateShapeFunctions(int rIP, const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());
    assert(rKnotIDs.rows() == 1 );

    Eigen::VectorXd IPcoordinates(1);

    IPcoordinates(0) = transformation(mIPCoordinates(rIP), mKnots(rKnotIDs(0)), mKnots(rKnotIDs(0) + 1));

    return CalculateShapeFunctions(IPcoordinates, rKnotIDs(0));
}

// --- derivatives shape functions --- //

Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd &rCoordinates) const
{
    int spanIdx = ShapeFunctionsIGA::FindSpan(rCoordinates(0,0), mDegree, mKnots);
    Eigen::MatrixXd shapeFunctionsAndDerivatives = ShapeFunctionsIGA::BasisFunctionsAndDerivatives(rCoordinates(0,0), spanIdx, 1, mDegree, mKnots);
    return (shapeFunctionsAndDerivatives.row(1)).transpose();
}


Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd &rCoordinates, const Eigen::VectorXi &rKnotIDs) const
{
    Eigen::MatrixXd shapeFunctionsAndDerivatives = ShapeFunctionsIGA::BasisFunctionsAndDerivatives(rCoordinates(0,0), rKnotIDs(0), 1, mDegree, mKnots);
    return (shapeFunctionsAndDerivatives.row(1)).transpose();
}

Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());
    assert(rKnotIDs.rows() == 1 );

    Eigen::VectorXd IPcoordinates(1);

    IPcoordinates(0) = transformation(mIPCoordinates(rIP), mKnots(rKnotIDs(0)), mKnots(rKnotIDs(0) + 1));

    return CalculateDerivativeShapeFunctionsNatural(IPcoordinates, rKnotIDs);
}

// --- N-matrix --- //

Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const
{
    auto shapeFunctions = CalculateShapeFunctions(rCoordinates);
    return shapeFunctions.transpose();
}


Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateMatrixN(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const
{
    auto shapeFunctions = CalculateShapeFunctions(rCoordinates, rKnotIDs(0));
    return shapeFunctions.transpose();
}


Eigen::MatrixXd NuTo::Interpolation1DIGA::CalculateMatrixN(int rIP,  const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());
    assert(rKnotIDs.rows() == 1 );

    Eigen::VectorXd IPcoordinates(1);

    IPcoordinates(0) = transformation(mIPCoordinates(rIP), mKnots(rKnotIDs(0)), mKnots(rKnotIDs(0) + 1));

    return CalculateMatrixN(IPcoordinates, rKnotIDs);
}


int NuTo::Interpolation1DIGA::CalculateNumNodes() const
{
    return mDegree+1;
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Interpolation1DIGA::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1DIGA::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation1DIGA::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1DIGA::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation1DIGA::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1DIGA::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Interpolation1DIGA::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation1DIGA" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBaseIGA);
    ar & boost::serialization::make_nvp("mKnots", mKnots);
    ar & boost::serialization::make_nvp("mDegree", mDegree);
    ar & boost::serialization::make_nvp("mIPCoordinates", mIPCoordinates);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation1DIGA" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation1DIGA)
#endif  // ENABLE_SERIALIZATION
