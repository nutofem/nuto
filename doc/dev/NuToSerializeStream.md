@page NuToSerializeStream NuToSerializeStream: How to use it.

# NuTo::SerializeStream

### Purpose
NuTo::SerializeStream provides a little framework to serialize data apart from boost::serialize. Its main usage will be the serialization of "restart data". This includes

 - nodal values (NuTo::BlockVector or NuTo::StructureOutputBlockVector)
 - static data (history data)
 - information from the time integration scheme

### Functionality

- output in binary / plain text via flag in constructor of NuTo::SerializeStreamIn / NuTo::SerializeStreamOut
- serialization "primitive types"
  - NuTo::SerializeStream::NuToSerializeNumber for double, int, enum, bool
  - NuTo::SerializeStream::NuToSerializeMatrix for Eigen::Matrix and NuTo::FullMatrix/NuTo::FullVector
- serialization of a class "A" (with the operator<< and operator>>) requires the implementation of
  - A::NuToSerializeWrite(NuTo::SerializeStreamOut&)
  - A::NuToSerializeRead(NuTo::SerializeStreamIn&)
  - these methods may call the serialization of other classes/numbers/matrices

### Warning
boost::serialize provides the luxury feature of automatically calling the __DerivedClass::serialize__ when serializing the __BaseClass__. This is __not__ the case for NuTo::SerializeStream. You will have to define the NuToSerializeWrite/Read as a virtual method in the base class and override it properly. When dealing with multiple inheritance, this involves explictly calling the base classes NuToSerializeWrite/Read.

~~~{.cpp}
class Base
{
public:
    virtual NuToSerializeWrite(SerializeStreamOut& rStream) const
    {
        rStream.NuToSerializeNumber(mData);
    }

    virtual NuToSerializeRead(SerializeStreamIn& rStream)
    {
        rStream.NuToSerializeNumber(mData);
    }
private:        
    double mData;
};

class Derived : public Base
{
public:
    virtual NuToSerializeWrite(SerializeStreamOut& rStream) const override
    {
        Base::NuToSerializeWrite(rStream);
        rStream.NuToSerializeMatrix(mData);
    }

    virtual NuToSerializeRead(SerializeStreamIn& rStream) override
    {
        Base::NuToSerializeRead(rStream);
        rStream.NuToSerializeMatrix(mData);
    }
private:
    Eigen::Matrix3d mData;
};
~~~
