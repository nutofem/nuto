@page ForwardDeclarations When to use forward declarations

What is a forward declaration?
------------------------------

Have a look at the following piece of Code

    class MyMatrix;
    class MatrixGenerator
    {
    public:
        MyMatrix GetAMatrix();
    };

The first line is a forward declaration. It tells the compiler, that there is
some class named `MyMatrix` somewhere in our program files, nothing more. We
can now use `MyMatrix` in some contexts without including its header file. This
is beneficial because every time we include something we increase the amount of
code in the current file, which needs to be compiled too. In addition, the
included header-file might also have some includes themselves. In the worst
case we end up including all the header files of our program just to declare
`MyMatrix` as the return type of one function. The result is a dramatically
increased compile time of our program.

Therefore you should always avoid includes in header-files and move them to the
implementation files (.cpp) . But there are some cases, when you can't use
forward declarations (see next sections).

**IMPORTANT: If your class is encapsulated in a special namespace, then your
forward declaration needs to be in this namespace too!**


When do I need includes?
------------------------

You need an `#include` if the compiler needs to know the size or functionality
of the class you want to use. This is usually the case if you want to use it as
member (NOT ALWAYS. Have a look at section: Pointer and References) of another
class, if you want to use it as base class of another class or if you want to
use functions or members of the class. The next piece of code shows you those
cases where you definitely need an include:

~~~{.cpp}
#include <iostream>
#include "B.h" // <--- parent
#include "C.h" // <--- Member
#include "D.h" // <--- used in function SomeFunction which is defined in the header

class A : public B
{
    A(); // ctor <--- defined in .cpp
    A(const A& rOther); // copy ctor <--- defined in .cpp
    A(A&& rOther); // move ctor <--- defined in .cpp
    ~A(); // destructor <--- defined in .cpp

     A& operator =(const A& rOther);     // copy assignment <--- defined in .cpp
     A& operator =(A&& rOther);          // move assignment <--- defined in .cpp

    D SomeFunction(D rParameter)
    {
        std::cout << rParameter.ToString() << std::endl;
        return rParameter;  
    }

    C mMember;
}
~~~

Here `B` is the parent of the class. So we need the include.  `mMember` is an
instance of `C `and not a pointer or reference.  Therefore we also need to
include C.h. `SomeFunction` is defined in the header and uses a member function
of `D`. The compiler needs to know if `D` really has a function called
`ToString` and therefore needs is definition. So `D.h` must be included.

These are basically the only 3 cases, where you need an include in your header
file. In every other case you can use forward declarations. **One exception
here are the standard containers.** The standard forbids to define stuff in the
`std` namespace. For this reason you can't use forward declarations for `std`
containers.

The next lines of code is the previous example with some minor changes,
to show you, when you don't need an include:

~~~{.cpp}
#include <iostream>
#include "B.h" // <--- parent

class C;
class D;

class A : public B
{
    A(); // ctor <--- defined in .cpp
    A(const A& rOther); // copy ctor <--- defined in .cpp
    A(A&& rOther); // move ctor <--- defined in .cpp
    ~A(); // destructor <--- defined in .cpp

     A& operator =(const A& rOther); // copy assignment <--- defined in .cpp
     A& operator =(A&& rOther); // move assignment <--- defined in .cpp

    D SomeFunction(D rParameter); // defined in .cpp

    C* mMember;
}
~~~

The definition of `B` still needs to be included, because it is the parent of
`A`. But `mMember` is now a pointer to a `C` object. So `C` can be forward
declared. The reason for this is, that a pointer has always the same size. It
doe not matter which object type it is pointing to, it just stores a memory
address. We also moved the implementation of `SomeFuntion` into the `A.cpp`. So
the only thing the compiler needs to know about `D` is that it exists, even
though it is passed and returned by value. The reason for this is that the
function interface does not need size or functionality informations from our
class `D`.

Standard container
------------------

As stated before, you can not use forward declarations with members of
the `std` namespace. So every container that appears anywhere in your
header needs to be included. But what you don't need to include are the
classes that the containers store. For example:

~~~{.cpp}
#include <vector>

class B;

class A
{
    A();                                
    A(const A& rOther);                  
    A(A&& rOther);                       
    ~A();

     A& operator =(const A& rOther);     
     A& operator =(A&& rOther);          

    std::vector<B> myBVec;
}
~~~

Here you need to include the `vector` (even if it wasn't in the
namespace `std` you to include the header because it is a member!). But
what you don't need is the definition of `B`, even though it is not
marked as a pointer or reference! This is true because the `std::vector`
internally uses pointers and therefore we can get away with a forward
declaration. But there are some pitfalls (See next section!).  While
this is the case for most standard containers, there might be some
exceptions. I have not tested it, but I am pretty sure that it does not
work for the `std::array` because its size is a compile time constant!

Pitfalls
--------

Read this section again, and again and again until you understand the
problem, otherwise you wont get it, when it occurs!

~~~{.cpp}
#include <vector>

class B;

class A
{
    A();                                
    A(const A& rOther);                  
    A(A&& rOther);                       
    ~A(){}

     A& operator=(const A& rOther);     
     A& operator=(A&& rOther) = default;          

    std::vector<B> myBVec;
}
~~~

The code above has two minor changes when compared to the example of the
standard container section. But those two changes will give you some very
cryptic linker errors. You will see a lot of template function calls and the
compile output will give you a lot of references to the standard headers.  So
what happens?

I implemented the destructor of our class `A` directly in the header (replaced
`;` by `{}` ) and I told the compiler to generate a default move assignment
operator. Both changes generate the same problem: Now we need to now something
about our class `B` in our header file! 

For example our destructor `~A()` automatically calls the destructor of
`myBVec`, because it is a class member of `A`. This is the destructor of the
`std::vector` which wants to delete all the current members of specified type
`B` dynamically!  But for dynamic memory allocations and deallocations you need
the full class definition to determine the memory footprint of the class.
Therefore you need the include.

If you use the default keyword for constructors, destructors or assignment
operators, it is the same thing. It works like writing the default
implementation directly into the header. The only way to avoid this linker
error is to implement the definitions in the .cpp file of your class.

**IMPORTANT: The compiler automatically generates copy/move constructors and
assignment operators as well as the destructor even though you don't ask him to
do so by using the default keyword.**

Some of them aren't always generated automatically (move/copy/ctor - ask google
when it is the case and when not) but if they are generated, they trigger the
error described above. So you have to implement their definitions in the .cpp
to solve this problem!  This such a hard to track error, because it is not
about what you have written. It is about what you have not written into your
code and the compiler messages aren't very helpful. 

NOTE: If you don't need the automatically generated functionality, the `delete`
keyword should also do the trick, because it suppresses the automatic code
generation.

Enum classes
------------

Enum classes can be forward declared like normal classes. Just write:

    enum class MyEnumClassName;

That's it. As long as you don't use any actual values of this class
(`MyEnumClassName::VALUE1`) you can do nearly everything with it without
including it, because the compiler knows, that is basically just an integer
data type. Only if the compiler needs to know which keyword is related with
which number, you need to include the definition.

If you changed the underlying data type of your enum class, you also have to do
that in the forward declaration:

    enum class
    MyEnumClassName : unsigned char;

NOTE: If you want to set a default enum class value in a function you should
probably implement the function twice (one function calls the other with
default value) instead of including the enum class definition in the header.

Templates
---------

Forward declarations for templates are also possible. You just need to
add the template stuff to the forward declaration:

    template <typename T, int rows, int cols> class MyMatrix;

Functions
---------

It's rare in NuTo, but there are sometimes functions, that do not belong to a
class. You can also use forward declarations on them:

    double MyFunction(int rParam);

In this case it shouldn't even be necessary to include the header in the .cpp
file, as long as the implementation is linked in later.

NuTo Specific Stuff
-------------------

Some template classes in NuTo (I think only the math classes) have an
additional header with `_Def.h` suffix. This is because of the problem that
very general templates usually can't be implemented in `.cpp` files (ask
google). The `_Def.h` file now contains the stuff of a normal header while the
`.h` file contains the definitions (usually the job of a .cpp). So if you can't
get away with a forward declaration and need a definition of such a class,
always use the `_def.h` file in your headers and not the `.h`. Otherwise you
also compile the definitions of the classes member functions, which is not
necessary in most cases and cost you a lot of time.

**But here is another pitfall:** Lets say you have included the
`FullMatrix_Def.h` in your class header. Now you write some implementations in
your classes `.cpp` that involve some operations of the FullMatrix. Because you
have the full definition included in your classes header (which should be
included in your `.cpp`), the compiler won't complain about anything, because
it now knows the FullMatrix-Class and it's interface. But until you are very,
very lucky, the linker will crash with some undefined reference
`FullMatrix.someFunction(...)` error. The reason is, that the definitions of
the FullMatrix class are not compiled. They are not defined in a `.cpp` file
and there is no object to link in. They are defined in  a `.h` file. That's why
the linker crashes! In order to get those definitions, you need to include the
`FullMatrix.h` in your `.cpp` file! 
