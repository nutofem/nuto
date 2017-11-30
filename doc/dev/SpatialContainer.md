@page SpatialContainer

# About

- convenient wrapper for ANN
  - [*A Library for Approximate Nearest Neighbor Searching*](https://www.cs.umd.edu/~mount/ANN/)
  - `annkSearch()` - find \f$k\f$ nearest neighbors
  - `annKFRSearch()` - find nearest neighbors within a radius

# Example: find duplicate nodes

```{.cpp}
// vector of relevant objects
std::vector<NuTo::NodeBase*> nodes = GetNodes();

// struct to extract a coordinate vector from that object
struct NodeCoordinate
{
    Eigen::VectorXd operator () (const NuTo::NodeBase* rNode)
    {
        return rNode->Get(NuTo::Node::eDof::COORDINATES);
    }
};

NuTo::SpatialContainer<NuTo::NodeBase*, NodeCoordinate> tree(nodes);

double searchRadius = 0.1;
std::vector<std::vector<NuTo::NodeBase*>> duplicates =  tree.GetAllDuplicateValues(searchRadius);

for (std::vector<NuTo::NodeBase*>& nodesAtSameCoordinate : duplicates)
    // do stuff...
```
