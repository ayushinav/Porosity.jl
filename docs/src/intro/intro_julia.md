# Julia 101

```@setup fn_demo
```

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

Usage of function is similar to MATLAB .In the above example, we already call functions. Julia also has a bunch of mutating functions. The name of such functions mostly end with a `!` and mutate the value stored in (usually) the first variables
```@repl
x = [30, 40, 10, 50, 20]
sort!(x)
x
```

Such functions are helpful when you want to reduce allocations to speed up the code.

### Creating new functions
Creating new functions is also similar to MATLAB, and is done via the following syntax :
```@repl
function f(x,y)
    a = x/(1 + x^2)
    b = y/(1 + y^2)
    z = a + b
    return z
end

f(2., 3.)
```

You can also do a single line definition : 
```@repl
f(x,y) = 2x + y

f(2., 3.)
```

### Optional and keyword arguments
A lot of times, we want default arguments or named arguments. Consider the following function definition: 
```@example fn_demo
function f(x,y,a=2; operation_type = 2)
    if operation_type == 1
        return x+y/a
    elseif operation_type == 2
        return x/(x^2 + a^2)
    elseif operation_type == 3
        return x/(x^2 + a^2)
    else
        return y/a
    end
end
```

In the above definition, `a` is called an *optional* argument. It assumes a value if it is not explicitly passed

!!! note

    Note that optional arguments are positional. If a function is defined as 
    ```julia
    function fn(x,y,a=1.,b=2.,c=3.)
        # some operation
    end
    ```
    Then the function call `fn(2., 3., 2., 4., 5.)` assumes `a = 2.`, `b = 4.` and `c = 5.`. Similarly, a function call `fn(2., 3., 4.)` assumes `a = 4.` and default values for `b` and `c`.

```@example fn_demo
f(2., 3.)
```
A different value for `a` can be passed without referring to it explicitly:
```@example fn_demo
f(2., 3., 1.)
```
`operation_type` is a *keyword* argument. In the above case, it also assumed a value but this is not necessary. If you want to pass a value to it, you have to explicitly refer to it in the function call:
```@example fn_demo
f(2., 3., operation_type = 2)
```
and with a different optional argument as :
```@example fn_demo
f(2., 3., 1., operation_type = 2)
```

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

While Julia has lot more features than we can cover here, this should help you get started with this package.
