@page NuToSerializeStream NuToSerializeStream: How to use it.

# NuTo::SerializeStream

### Purpose

NuTo::SerializeStream provides a little framework to serialize data apart from boost::serialize. Its main usage will be the serialization of "restart data". This includes

- nodal values (NuTo::BlockVector or NuTo::StructureOutputBlockVector)
- static data (history data)
- information from the time integration scheme

### Functionality

- output in binary / plain text via flag in constructor of NuTo::SerializeStreamIn / NuTo::SerializeStreamOut
- serialization with the **Serialize(rData)** method for both saving and loading
  - You should use this method as seen in the example. This ensures saving and loading in the same way.
  - Downside: The rData argument cannot be constant anymore... 
- additional wrapping with **operator &lt;&lt;** for saving and **operator >>** for restoring
- **Serialize()** is implemented for the types **double** and **Eigen::Matrix&lt;...>**
- every other class needs to implement the methods
  - NuToSerializeSave(NuTo::SerializeStreamOut&) const
  - NuToSerializeLoad(NuTo::SerializeStreamIn&)
- (both are implemented for NuTo::FullMatrix&lt;...>)

### Warning

boost::serialize provides the luxury feature of automatically calling the **DerivedClass::serialize** when serializing the **BaseClass**. This is **not** the case for NuTo::SerializeStream. You will have to define the NuToSerializeWrite/Read as a virtual method in the base class and override it properly. When dealing with multiple inheritance, this involves explictly calling the base classes NuToSerializeWrite/Read.

```{.cpp}
class Base
{
public:
    virtual NuToSerializeSave(SerializeStreamOut& rStream)
    {
        Serialize(rStream);
    }

    virtual NuToSerializeLoad(SerializeStreamIn& rStream)
    {
        Serialize(rStream);
    }
private:    
    
    template <typename TStream>
    void SerializeMe(TStream& rStream)
    {
        rStream.Serialize(mData);
    }
    
    double mData;
};

class Derived : public Base
{
public:
    virtual NuToSerializeSave(SerializeStreamOut& rStream) override
    {
        Base::NuToSerializeSave(rStream);
        SerializeMe(rStream);
    }

    virtual NuToSerializeLoad(SerializeStreamIn& rStream) override
    {
        Base::NuToSerializeLoad(rStream);
        SerializeMe(rStream);
    }
private:
    template <typename TStream>
    void SerializeMe(TStream& rStream)
    {
        rStream.Serialize(mData);
    }


    Eigen::Matrix3d mData;
};
```

### Separator

The stream continuously reads and writes numbers. If there are more numbers in the stream than the object needs, this will not be noticed. Therefore, you can write a separator of type char at certain points in the serialization, e.g. after serializing the data of one element. An NuTo::Exception is thrown, if the expected separator is not at the right position in the stream.
