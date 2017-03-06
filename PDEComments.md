# Some comments on the PDE implementation

## Intro

- have a look at `/test/mechanics/cell/Cell.cpp` for an unit test that builds the internal gradient and extracts integration point values


## TODO

- renaming 
    - some classes got the suffix `Simple` to avoid clashes with existing classes/namespaces
- fix resize of DofContainer to "magicContainerSize = 10"
- add automatic assignment of a DofType::mId (currently done manually)
    - caution: correctly handle copy/move operations
- mechanics/nodes/DofVector is currently derived from DofContainer. A more elegant approach is to use a Eigen::Vector as the underlying data structure and a clever use of of the Eigen::VectorXd::segment() method to access the dof types. This makes arithmetics trivial.
- Hessian0/1/2: needs another kind of `DofContainer`, namely `DofMatrix`, should be trivial apart from that
- more complex tests
    - dof numbering - IMO not in `NuTo::NodeSimple` but in another class that accesses `NuTo::NodeSimple`, maybe `NuTo::Assembler`?
    - add constraints
    - `NuTo::Mesh` that only holds `NuTo::ElementSimple` for coordinates
- add the new implementation to `application/Benchmark/BuildInternalGradient.cpp`
- think about the general case of $n$ dimensional interpolation in $n+x$ dimensional space

## more or less trivial tasks

- build a *InterpolationType* library (implements `NuTo::Interpolation`)
- build a *ConstitutiveLaw* library (few base classes, maybe `NuTo::Constitutive::MechanicsBase`)
- build a *Integrand* library for common PDEs
- modify interface of `IntegrationTypeBase` to something based on iterators

~~~{.cpp}
for (const auto& ip : mIntegrationType)
{
    ... = ...GetShapeFunctions(ip.coords);
    ... += ... *ip.weight;
}
~~~

## Visualization

- `NuTo::VisualizationCell` - data structure that saves the geometry of a integration point visualization (one `NuTo::Cell` containts `GetNumIPs()` many `NuTo::VisualizationCell`s)
- extract and save `NuTo::VisualizationCell`s for each cell and for each IP
    - store some kind of identifier for the cell and the integration point
    - optional: run ANN to remove duplicate visualization points
- use the `NuTo::IPValue` to visualize the data only
    - two members: `std::string mName` and `Eigen::MatrixXd mValues`
    - visualize everthing with the same `mName` in the same data set?

## Cache

We should think about caching / memoization approaches that do to clutter the interface. Currently, we store the ShapeFunctions via their index in an integration type. This leads to unneccesary methods like `CalculateShapeFunctions(coordinates)` and `GetShapeFunctions(integrationPointIndex)`. On top of that, the InterpolationType needs to now an IntegrationType to precalculate N and dN for each integration point. And on another top of that, this storage is done in the InterpolationType. This means that ShapeFunctions and DerivativeShapeFunctions are stored instead of the N and B(local). The latter ones mean the 'blown' versions. 

My proposal: Write a caching class. Use it on Element level since the element knows how to build N and B*local* based on its DofType/PDEType. API could look like this.

~~~{.cpp}
const Eigen::MatrixXd& Element::N(Eigen::VectorXd rLocalIpCoord) const
{
    auto cacheIterator = mCacheN.Get(rLocalIPCoord);
    if (cacheIterator == mCacheN.end())
    {
        *cacheIterator = Eigen::MatrixXd(...);
        // build N for rLocalIpCoord;
    }
    return *cacheIterator;
}
~~~

- Advantage: Everything is hidden.
- use `mutable` to keep constness 



## Edge cases

### Nonlocal integral model

- How to pass information of neighboring cells to the integrand?


### Change of time integration scheme

In the current idea, an `Integrand` is derived from `IntegrandNewmark` **or** `IntegrandRungeKutta`. Problem: The history data is stored at the `Integrand` and needs to be copied(?) from Newmark to RungeKutta.

### Concerns of Peter

- IGA
- Contact


## Requests

- please provide some kind of documentation. UnitTests are great for that...
- get familiar with [`fakeit::Mock<...>`](https://github.com/eranpeer/FakeIt/wiki/Quickstart) and try to use it to write *smaller* (less includes) unit test
