@page AvoidUnsigned Using `unsigned int` is almost never a good idea

To quote Stroustrup (C++ Programming Language, 4.ed, p 73):

> Using an `unsigned` instead of ant `int` [...] is almost never a good idea.
> Attemps to ensure that some values are positive by declaring variables
> `unsigned` will typically be defeated by the implicit conversion rules.

Scott Meyers agrees (obviously), and has written a short paper that better
explains the pitfalls: "Signed and unsigned types in interfaces," C++ Report,
September 1995.
