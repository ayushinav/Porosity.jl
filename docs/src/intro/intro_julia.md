# Julia 101

If this is the first time you've come across [Julia](https://julialang.org/), this page will provide a basic tutorial into Julia. Much of the syntax is similar to MATLAB and Python.

## Basics

The basics of operation are very similar. 

```@repl
x = 10

y = x + 3
y = x * 3
y = x / 3
y = x ÷ 3
y = x % 3

for i in 1:5
    println("2 x i = ", 2i)
end
```

## Arrays

Like MATLAB, array indices start with 1

```@repl
x = [10,20,30,40,50]

a = x[1]
b = x[2]
```

Like MATLAB, to perform element-wise operations on arrays, use dot in front of the operator. This performs broadcasting and will be useful for arrays of different sizes

```@repl
x = [10,20,30,40,50]

y = x .+ 5
y = x./10
y = x .+ x'
```

When functions and complex expressions are needed to be performed, `@.` can be used to get a cleaner code.

```@repl
x = [10,20,30,40,50]
y = exp.(sin.(x./10))
y = @. exp(sin(x/10))
```

## Functions

Usage of function is similar to, in the above example, we already call functions. Julia also has a bunch of mutating functions. The name of such functions mostly end with a `!` and mutate the value stored in the first (or more) variables
```@repl

x = [30, 40, 10, 50, 20]
sort!(x)
x
```

Such functions are helpful when you want to reduce allocations to speed up the code.

## Named Tuple

Named Tuples are data structures that hold values assigned to them by a name. These are the most important data structures in the context of the package and are used a lot.

```@repl

x = (; a = 30, b = 2., c= "Hello")

x.a
x.b
x.c
```
Named-tuples are immutable, that is, once created, their values do not change. However, we can create another named-tuple of the same name that can feel like mutating.

```@repl

x = (; a = 30, b = 2., c= "Hello")
x.a
x = (; x..., a = 50)
x.a
```

While Julia has lot more features that we can cover here, this should help you get started with this package.
