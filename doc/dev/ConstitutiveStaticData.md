@page ConstitutiveStaticData ConstitutiveStaticData: Some design choices



### Setup
- each constitutive law (NuTo::ConstitutiveBase) is coupled with its corresponding NuTo::IPConstitutiveLawBase
- there is one instance of NuTo::IPConstitutiveLawBase at each integration point
- it is created via the constitutive law with the method
~~~{.cpp}
virtual std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> 
   NuTo::ConstitutiveBase::CreateIPLaw() = 0;
~~~

- CASE A: constitutive laws __with__ static data:
    - they return a NuTo::Constitutive::IPConstitutiveLaw\<SpecificLaw\> in CreateIPLaw() 
    - they define the public typedef _StaticDataType_
~~~{.cpp}
typedef <specific static data type> StaticDataType
// e.g. NuTo::MoistureTransport:
typedef NuTo::Constitutive::StaticData::DataMoistureTransport StaticDatatType
// e.g. NuTo::GradientDamageEngineeringStress:
typedef double StaticDataType
~~~

- CASE B: constitutive laws __without__:
    - they return a NuTo::Constitutive::IPConstitutiveLawWithoutData\<SpecificLaw\> in CreateIPLaw() 
    - they do not need the typedef

Benefits:
- the IPConstitutiveLawBase automatically allocates, deallocates and provides the correct static data, if the constitutive law needs them.

### Evaluate()
The element calls IPConstitutiveLawBase::Evaluate(inout, output). This calls is forwarded to the constitutive law providing an additional vector of _StaticDataType_, if needed.
 
 
### Access
via 
~~~{.cpp}
template <typename TLaw>
DataContainer<TLaw> NuTo::Contitutive::IPConstitutiveLawBase::GetData();
~~~
The class DataContainer<TLaw> contains accessors for the data. 

E.g:
~~~{.cpp}
ElementBase* e;
unsigned ipIndex = 2;
e->GetIPData().GetIPConstitutiveLaw(ipIndex).GetData<NuTo::PhaseField>.SetData(42.);
~~~

### Serialization

see @ref NuToSerializeStream

Make sure that your _StaticDataType_ is
 - serializable by default, that is
    - Eigen::Matrix/Vector
    - NuTo::FullMatrix/FullVector
    - double, bool, int
 - provides a NuToSerializeSave(...) and NuToSerializeLoad(...)