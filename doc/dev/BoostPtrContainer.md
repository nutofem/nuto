@page BoostPtrContainer `boost::ptr_container` vs. `std::container` of `std::unique_ptr`

**Note: If you need shared semantics, neither is the right solution. Use a
`std::container` of `std::shared_ptr`s**

### Advantages of `boost::ptr_container`

#### Value semantics:

```
boost::ptr_vector<MyClass> vec;
vec.push_back(new MyClass);
vec[0].foo();
```

vs. pointer semantics with unique pointers

```
std::vector<std::unique_ptr<MyClass>> vec;
vec.push_back(std::make_unique<MyClass>());
vec[0]->foo();
```

#### Propagates constness

```
void f(const boost::ptr_vector<MyClass> vec)
{
    vec[0].SetMember(42.0); // Compile time error
}
```

vs.

```
void g(const std::vector<std::unique_ptr<MyClass>>& vec)
{
    vec[0]->SetMember(42.0); // will compile and run
}
```

#### Built-in support for deep copy semantics

For a `std::container` of unique pointers, you would have to implement a deep
copy yourself. In Boost, it is built in, using the
[Clonable Concept](http://www.boost.org/doc/libs/1_61_0/libs/ptr_container/doc/reference.html#the-clonable-concept)

#### Lower memory footprint and usually faster

See [Boost Ptr Container Documentation](http://www.boost.org/doc/libs/1_61_0/libs/ptr_container/doc/ptr_container.html#motivation).

### Advantages of `std::container<std::unique_ptr<MyClass>>`

#### No shallow copy possible

You can't make a shallow copy, because the `unique_ptr` is not copy
constructible. This is good. With boost pointer containers, you could have copy
of the container with dangling pointers after the original container went out
of scope.

**Don't make shallow copies of owning containers, regardless of how they are
implemented.**
