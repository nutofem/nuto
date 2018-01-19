#include "mechanics/mesh/MeshFemDofConvert.h"

//! @brief container that performs O(1) addition and O(1) lookup based on subboxes. These subboxes divide the domain
//! with a given number of subdivisions.
//! @tparam value type of the container that provides coordinate access `double operator[](int dim)`
template <typename T>
class SubBoxes
{
public:
    struct Domain
    {
        //! @brief lowest coordinate values of the domain
        Eigen::Vector3d start;

        //! @brief highest coordinate values of the domain
        Eigen::Vector3d end;

        //! @brief number of subdivisions
        Eigen::Vector3i subs;

        //! @brief radius should be in the range of the smallest distance beween two nodes
        double searchRadiusSquared;
    };

    SubBoxes(Domain d)
        : mDomain(d)
        , mData(d.subs.prod())
        // strange syntax to do (end - start) / subs
        , mDelta((mDomain.end - mDomain.start).cwiseQuotient(mDomain.subs.template cast<double>()))
    {
        if (d.searchRadiusSquared > mDelta.minCoeff() * mDelta.minCoeff())
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "Search radius too small. Internal error. Choose a smaller eps.");
    }

    //! @brief adds a T to the container
    void Add(T t)
    {
        mData.at(Index(t)).push_back(t);
    }

    //! @brief returns _some_ entry at `coords` within a search radius (provided by the domain)
    //! @remark This is one of the reasons why this class should not be public. This behavior is very specific to the
    //! algorithm below that only adds a T if nothing is found within the search radius.
    //! @param coords query point
    //! @return pointer to T if something is found, nullptr otherwise
    T* FindAt(Eigen::Vector3d coords)
    {
        int i = (coords[0] - mDomain.start[0]) / mDelta[0];
        int j = (coords[1] - mDomain.start[1]) / mDelta[1];
        int k = (coords[2] - mDomain.start[2]) / mDelta[2];

        // The search radius is smaller than the subbox size. So we have to search the 9 surrounding subboxes. A newly
        // added entry would go to (i,j,k) which makes this subbox most likely for a match.
        for (int kk : {k, k - 1, k + 1})
            for (int jj : {j, j - 1, j + 1})
                for (int ii : {i, i - 1, i + 1})
                {
                    // skip the invalid indices at the boundary of the domain
                    if (ii < 0 || ii >= mDomain.subs[0])
                        continue;
                    if (jj < 0 || jj >= mDomain.subs[1])
                        continue;
                    if (kk < 0 || kk >= mDomain.subs[2])
                        continue;

                    // search the subbox (ii,jj,kk) for matches
                    for (T& entry : mData[Index(ii, jj, kk)])
                        if ((entry - coords).squaredNorm() < mDomain.searchRadiusSquared)
                            return &entry;
                }
        return nullptr; // TODO17:  This is a classic case for the c++17 feature std::optional<T>. Here we have to
        // introduce pointers to our precious interface.
    }

private:
    int Index(const T& t) const
    {
        int i = (t[0] - mDomain.start[0]) / mDelta[0];
        int j = (t[1] - mDomain.start[1]) / mDelta[1];
        int k = (t[2] - mDomain.start[2]) / mDelta[2];
        return Index(i, j, k);
    }

    int Index(int i, int j, int k) const
    {
        return i + j * mDomain.subs[0] + k * mDomain.subs[1] * mDomain.subs[0];
    }

    Domain mDomain;
    std::vector<std::vector<T>> mData;
    const Eigen::Vector3d mDelta; // mDelta is precalculated since it is used a lot.
};


Eigen::Vector3d To3D(Eigen::VectorXd v)
{
    Eigen::Vector3d v3d = Eigen::Vector3d::Zero();
    v3d.segment(0, v.size()) = v;
    return v3d;
}

struct NodePoint : Eigen::Vector3d // to inherit operator[]
{
    NodePoint(Eigen::VectorXd v, NuTo::NodeSimple& node)
        : Eigen::Vector3d(To3D(v))
        , mNode(node)
    {
    }
    NuTo::NodeSimple& mNode;
};

SubBoxes<NodePoint>::Domain SetupSubBoxDomain(const NuTo::MeshFem& mesh, int numBoxesPerDirection, double eps)
{
    Eigen::VectorXd start = mesh.Elements[0].CoordinateElement().GetNode(0).GetValues();
    Eigen::VectorXd end = start;

    unsigned long numNodesTotal = 0;
    for (auto& elementCollection : mesh.Elements)
    {
        const int numNodes = elementCollection.CoordinateElement().GetNumNodes();
        numNodesTotal += numNodes;

        for (int iNode = 0; iNode < numNodes; ++iNode)
        {
            Eigen::VectorXd coords = elementCollection.CoordinateElement().GetNode(iNode).GetValues();
            start = start.array().min(coords.array()); // array view (.array()) allows coeffwise operations
            end = end.array().max(coords.array());
        }
    }
    const int dimension = start.rows();

    SubBoxes<NodePoint>::Domain d;
    d.start = To3D(start) - Eigen::Vector3d::Constant(eps);
    d.end = To3D(end) + Eigen::Vector3d::Constant(eps);

    d.subs = Eigen::Vector3i::Constant(1);

    for (int i = 0; i < dimension; ++i)
        d.subs[i] = numBoxesPerDirection;

    d.searchRadiusSquared = eps * eps;
    return d;
}


void NuTo::AddDofInterpolation(MeshFem* rMesh, DofType dofType, const InterpolationSimple& interpolation)
{
    // Setup subbox. These argument values are quite arbibrary and should maybe be chosen based on
    // the dimensions of the mesh.
    SubBoxes<NodePoint> subBoxes(SetupSubBoxDomain(*rMesh, /*numBoxesPerDirection=*/300, /*eps=*/1.e-10));

    for (auto& elementCollection : rMesh->Elements)
    {
        std::vector<NodeSimple*> nodesForTheNewlyCreatedElement;

        const auto& coordinateElement = elementCollection.CoordinateElement();
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            Eigen::Vector3d coord = To3D(Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode)));

            NodePoint* nodePoint = subBoxes.FindAt(coord);
            if (nodePoint)
            {
                nodesForTheNewlyCreatedElement.push_back(&nodePoint->mNode);
            }
            else
            {
                auto& node = rMesh->Nodes.Add(Eigen::VectorXd::Zero(dofType.GetNum()));
                subBoxes.Add(NodePoint(coord, node));
                nodesForTheNewlyCreatedElement.push_back(&node);
            }
        }
        elementCollection.AddDofElement(dofType, ElementFem(nodesForTheNewlyCreatedElement, interpolation));
    }
}

void NuTo::ChangeCoordinateInterpolation(MeshFem* rMesh, const InterpolationSimple& interpolation)
{
    // Setup subbox. These argument values are quite arbibrary and should maybe be chosen based on
    // the dimensions of the mesh.
    SubBoxes<NodePoint> subBoxes(SetupSubBoxDomain(*rMesh, /*numBoxesPerDirection=*/300, /*eps=*/1.e-10));

    // Save old coordinate node refs to remove them later
    auto oldCoordinateNodes = rMesh->NodesTotal();

    for (size_t i = 0; i < rMesh->Elements.Size(); i++)
    {
        auto& elementCollection = rMesh->Elements[i];
        std::vector<NodeSimple*> nodesForTheNewlyCreatedElement;

        const auto& coordinateElement = elementCollection.CoordinateElement();
        for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
        {
            Eigen::Vector3d coord = To3D(Interpolate(coordinateElement, interpolation.GetLocalCoords(iNode)));

            NodePoint* nodePoint = subBoxes.FindAt(coord);
            if (nodePoint)
            {
                nodesForTheNewlyCreatedElement.push_back(&nodePoint->mNode);
            }
            else
            {
                auto& node = rMesh->Nodes.Add(coord.head(coordinateElement.GetDofDimension()));
                subBoxes.Add(NodePoint(coord, node));
                nodesForTheNewlyCreatedElement.push_back(&node);
            }
        }
        elementCollection.CoordinateElement() = ElementFem(nodesForTheNewlyCreatedElement, interpolation);
    }

    for (auto& nd : oldCoordinateNodes)
    {
        auto isIdentical = [&](NodeSimple& n) { return ((&n) == (&nd)); };
        auto& nodes = rMesh->Nodes;
        nodes.Erase(std::find_if(nodes.begin(), nodes.end(), isIdentical));
    }
}
