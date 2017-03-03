# Some comments on the PDE implementation


## TODO

- fix resize of DofContainer to "magicContainerSize = 10"
- add automatic assignment of a DofType::mId (currently done manually
- The approach with BMatrixLocal was tremendously stupid. The transformation to BMatrixGlobal cannot be done afterwards. This has to be moved to `NuTo::CellIPData`. However, it could be done for the `N` matrix...
    - new typedefs: DerivativeShapeFunctionsLocal --> DerivativeShapeFunctionsGlobal --> BMatrixStrain/BMatrixGradient

## Cache

We should think about caching / memoization approaches that do to clutter the interface. Currently, we store the ShapeFunctions via their index in an integration type. This leads to unneccesary methods like `CalculateShapeFunctions(coordinates)` and `GetShapeFunctions(integrationPointIndex)`. On top of that, the InterpolationType needs to now an IntegrationType to precalculate N and dN for each integration point. And on another top of that, this storage is done in the InterpolationType. This means that ShapeFunctions and DerivativeShapeFunctions are stored instead of the N and B(local). The latter ones mean the 'blown' versions. 

My proposal: Write a caching class. Use it on Element level since the element knows how to build N and B*local* based on its DofType/PDEType. API could look like this.

~~~{.cpp}
const Eigen::MatrixXd& Element::N(Eigen::VectorXd rLocalIpCoord)
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
Advantage: Everything is hidden.

Disadvantage: Method not `const`



## Edge cases

### Nonlocal integral model

- How to pass information of neighboring cells to the integrand?


### Change of time integration scheme

In the current idea, an `Integrand` is derived from `IntegrandNewmark` **or** `IntegrandRungeKutta`. Problem: The history data is stored at the `Integrand` and needs to be copied(?) from Newmark to RungeKutta.

### Concerns of Peter

- IGA
- Contact

