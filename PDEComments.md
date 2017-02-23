# Some comments on the PDE implementation


## TODO

- fix resize of DofContainer to "magicContainerSize = 10"
- add automatic assignment of a DofType::mId (currently done manually)


## Edge cases

### Nonlocal integral model

- How to pass information of neighboring cells to the integrand?


### Change of time integration scheme

In the current idea, an `Integrand` is derived from `IntegrandNewmark` **or** `IntegrandRungeKutta`. Problem: The history data is stored at the `Integrand` and needs to be copied(?) from Newmark to RungeKutta.
