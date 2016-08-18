#include "nuto/mechanics/interpolationtypes/Interpolation2DIGA.h"

NuTo::Interpolation2DIGA::Interpolation2DIGA(NuTo::Node::eDof rDofType,
                                             NuTo::Interpolation::eTypeOrder rTypeOrder,
                                             int rDimension,
                                             const Eigen::Vector2i &rDegree,
                                             const Eigen::VectorXd &rKnotsX,
                                             const Eigen::VectorXd &rKnotsY)
                        :
                                  InterpolationBaseIGA::InterpolationBaseIGA(rDofType, rTypeOrder, rDimension),
                                  mDegree(rDegree),
                                  mKnotsX(rKnotsX),
                                  mKnotsY(rKnotsY)
{
    Initialize();
}

NuTo::IntegrationType::eIntegrationType NuTo::Interpolation2DIGA::GetStandardIntegrationType() const
{
    // its tensor product
    switch (std::max(mDegree(0), mDegree(1)))
    {
    case 0: // (0+0+1)/2 = 0.5 ips
        return NuTo::IntegrationType::IntegrationType2D4NGauss1Ip;
    case 1: // (1+1+1)/2 = 1.5 ips or (1+1+3)/2 = 2.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType2D4NGauss4Ip;
    case 2: // (2+2+1)/2 = 2.5 ips or (2+2+3)/2 = 3.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType2D4NGauss9Ip;
    case 3: // (3+3+1)/2 = 3.5 ips or (3+3+3)/2 = 4.5 ips lobatto
        return NuTo::IntegrationType::IntegrationType2D4NLobatto25Ip;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + std::to_string(std::max(mDegree(0), mDegree(1))) + " IGA not implemented");
    }
}

std::vector<Eigen::VectorXd> NuTo::Interpolation2DIGA::GetSurfaceEdgesCoordinates(int rSurface) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented yet");
}

int NuTo::Interpolation2DIGA::GetNumDofsPerNode() const
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

/*******************/
// Shape functions //
/*******************/

Eigen::VectorXd NuTo::Interpolation2DIGA::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    Eigen::Vector2i spanIdx;
    spanIdx(0) = ShapeFunctionsIGA::FindSpan(rCoordinates(0), mDegree(0), mKnotsX);
    spanIdx(1) = ShapeFunctionsIGA::FindSpan(rCoordinates(1), mDegree(1), mKnotsY);

    return ShapeFunctionsIGA::BasisFunctions2D(rCoordinates, spanIdx, mDegree, mKnotsX, mKnotsY);
}


Eigen::VectorXd NuTo::Interpolation2DIGA::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates, const Eigen::Vector2i &rKnotIDs) const
{
    assert(rKnotIDs.rows() == 2 );
    assert(rKnotIDs(0)-1 < mKnotsX.rows());
    assert(rKnotIDs(1)-1 < mKnotsX.rows());

    Eigen::Vector2d parameter;

    parameter(0) = transformation(rCoordinates(0), mKnotsX(rKnotIDs(0)), mKnotsX(rKnotIDs(0) + 1));
    parameter(1) = transformation(rCoordinates(1), mKnotsY(rKnotIDs(1)), mKnotsY(rKnotIDs(1) + 1));

    return ShapeFunctionsIGA::BasisFunctions2D(parameter, rKnotIDs, mDegree, mKnotsX, mKnotsY);
}

Eigen::VectorXd NuTo::Interpolation2DIGA::CalculateShapeFunctions(int rIP, const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());

    Eigen::Vector2d IPcoordinates;

    IPcoordinates << mIPCoordinates(rIP, 0), mIPCoordinates(rIP, 1);

    return CalculateShapeFunctions(IPcoordinates, rKnotIDs);
}

/***************/
// Derivatives //
/***************/
Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd &rCoordinates) const
{
    Eigen::Vector2i spanIdx;
    spanIdx(0) = ShapeFunctionsIGA::FindSpan(rCoordinates(0), mDegree(0), mKnotsX);
    spanIdx(1) = ShapeFunctionsIGA::FindSpan(rCoordinates(1), mDegree(1), mKnotsY);

    Eigen::Matrix<double, Eigen::Dynamic, 2> shapeFunctionsAndDerivatives = ShapeFunctionsIGA::BasisFunctionsAndDerivatives2D(rCoordinates, spanIdx, mDegree, mKnotsX, mKnotsY);
    return shapeFunctionsAndDerivatives;
}

Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd &rCoordinates, const Eigen::VectorXi &rKnotIDs) const
{
    assert(rKnotIDs.rows() == 2);
    assert(rKnotIDs(0)-1 < mKnotsX.rows());
    assert(rKnotIDs(1)-1 < mKnotsX.rows());

    Eigen::Vector2d parameter;

    parameter(0) = transformation(rCoordinates(0), mKnotsX(rKnotIDs(0)), mKnotsX(rKnotIDs(0) + 1));
    parameter(1) = transformation(rCoordinates(1), mKnotsY(rKnotIDs(1)), mKnotsY(rKnotIDs(1) + 1));

    Eigen::Matrix<double, Eigen::Dynamic, 2> shapeFunctionsAndDerivatives = ShapeFunctionsIGA::BasisFunctionsAndDerivatives2D(parameter, rKnotIDs, mDegree, mKnotsX, mKnotsY);
    return shapeFunctionsAndDerivatives;
}

Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());

    Eigen::Vector2d IPcoordinates;
    IPcoordinates << mIPCoordinates(rIP, 0), mIPCoordinates(rIP, 1);

    return CalculateDerivativeShapeFunctionsNatural(IPcoordinates, rKnotIDs);
}

/***************/
// N-Matrix    //
/***************/
Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const
{
    auto shapeFunctions = CalculateShapeFunctions(rCoordinates);

    return ConstructMatrixN(shapeFunctions);
}

Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateMatrixN(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const
{
    assert(rKnotIDs.rows() == 2 );
    assert(rKnotIDs(0)-1 < mKnotsX.rows());
    assert(rKnotIDs(1)-1 < mKnotsX.rows());

    Eigen::Vector2d parameter;

    parameter(0) = transformation(rCoordinates(0), mKnotsX(rKnotIDs(0)), mKnotsX(rKnotIDs(0) + 1));
    parameter(1) = transformation(rCoordinates(1), mKnotsY(rKnotIDs(1)), mKnotsY(rKnotIDs(1) + 1));

    auto shapeFunctions = ShapeFunctionsIGA::BasisFunctions2D(parameter, rKnotIDs, mDegree, mKnotsX, mKnotsY);

    return ConstructMatrixN(shapeFunctions);
}

Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateMatrixN(int rIP, const Eigen::VectorXi &rKnotIDs) const
{
    assert(!mUpdateRequired);
    assert(rIP >=0 && rIP < mIPCoordinates.rows());

    Eigen::Vector2d IPcoordinates;
    IPcoordinates << mIPCoordinates(rIP, 0), mIPCoordinates(rIP, 1);

    return CalculateMatrixN(IPcoordinates, rKnotIDs);
}

Eigen::MatrixXd  NuTo::Interpolation2DIGA::ConstructMatrixN(Eigen::VectorXd rShapeFunctions) const
{
    int numNodes = GetNumNodes();
    int dimBlock = GetNumDofsPerNode();

    assert (rShapeFunctions.rows() == (mDegree(0)+1)*(mDegree(1)+1));

    Eigen::MatrixXd matrixN(dimBlock, numNodes * dimBlock);
    for (int iNode = 0, iBlock = 0; iNode < numNodes; ++iNode, iBlock += dimBlock)
    {
        matrixN.block(0, iBlock, dimBlock, dimBlock) = Eigen::MatrixXd::Identity(dimBlock, dimBlock) * rShapeFunctions(iNode);
    }

    return matrixN;
}

Eigen::VectorXd NuTo::Interpolation2DIGA::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);

    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(transformation(rNaturalSurfaceCoordinates(0), rKnots(0,0), rKnots(0,1)), rKnots(1,0));
    case 1:
        return Eigen::Vector2d(rKnots(0,1), transformation(rNaturalSurfaceCoordinates(0), rKnots(1,0), rKnots(1,1)));
    case 2:
        return Eigen::Vector2d(rKnots(0,1) - transformation(rNaturalSurfaceCoordinates(0), rKnots(0,0), rKnots(0,1)), rKnots(1,1));
    case 3:
        return Eigen::Vector2d(rKnots(0,0), rKnots(1,1) - transformation(rNaturalSurfaceCoordinates(0), rKnots(1,0), rKnots(1,1)));
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA2D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation2DIGA::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);
    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(1, 0);
    case 1:
        return Eigen::Vector2d(0, 1);
    case 2:
        return Eigen::Vector2d(-1, 0);
    case 3:
        return Eigen::Vector2d(0, -1);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA2D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
}

//       2
//    -------
//  3 |     | 1
//    |     |
//    -------
//       0
Eigen::VectorXi NuTo::Interpolation2DIGA::GetSurfaceNodeIndices(int rSurface) const
{
    Eigen::VectorXi indices;

    switch (rSurface)
    {
    case 0:
    {
        indices.resize(mDegree(0) + 1);
        for(int i = 0; i <= mDegree(0); i++) indices(i) = i;
        return indices;
    }
    case 1:
    {
        indices.resize(mDegree(1) + 1);
        for(int i = 1; i <= mDegree(1)+1; i++) indices(i-1) = (mDegree(0)+1) * i - 1;
        return indices;
    }
    case 2:
    {
        indices.resize(mDegree(0) + 1);
        int count = 0;
        for(int i = CalculateNumNodes() - mDegree(0) - 1 ; i < CalculateNumNodes(); i++, count++) indices(count) = i;
        return indices;
    }
    case 3:
    {
        indices.resize(mDegree(1) + 1);
        for(int i = 0; i <= mDegree(1); i++) indices(i) = (mDegree(0)+1) * i;
        return indices;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA2D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
}

int NuTo::Interpolation2DIGA::GetSurfaceDegree(int rSurface) const
{
    switch (rSurface)
    {
    case 0:
        return mDegree(0);
    case 1:
        return mDegree(1);
    case 2:
        return mDegree(0);
    case 3:
        return mDegree(1);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA2D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
}

void NuTo::Interpolation2DIGA::UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType)
{
    int dim = rIntegrationType.GetCoordinateDimension();
    assert(dim == 2);

    int numIPs = rIntegrationType.GetNumIntegrationPoints();

    mIPCoordinates.resize(numIPs, 2);

    for (int iIP = 0; iIP < numIPs; ++iIP)
    {
        double coordinate2D[2];
        rIntegrationType.GetLocalIntegrationPointCoordinates2D(iIP, coordinate2D);
        mIPCoordinates(iIP, 0) = coordinate2D[0];
        mIPCoordinates(iIP, 1) = coordinate2D[1];
    }

    mUpdateRequired = false;
}

int NuTo::Interpolation2DIGA::CalculateNumNodes() const
{
    return (mDegree(0)+1)*(mDegree(1)+1);
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Interpolation2DIGA::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation2DIGA::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation2DIGA::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation2DIGA::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation2DIGA::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation2DIGA::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Interpolation2DIGA::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation1DIGA" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBaseIGA);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation1DIGA" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation2DIGA)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Interpolation2DIGA)
#endif  // ENABLE_SERIALIZATION
