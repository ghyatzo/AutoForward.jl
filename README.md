# Generalised API

## Basic Usage
the wrapper will transparently behave as the specified type (if a single symbol)
or as a tuple of the the specified types that automatically splats:
```julia
@forward <T> struct W <: MaybeSubtyped
	...
	fields
	...
end (<M>)
```
Here `<T>` can be different things:

1. A single type:
	1. `P => :p ::Pair{DataType, Symbol}`
	2. `P::DataType`: In this case it is required that the struct `W` has a **single field of type** `P`, then `P` can be univocally made equivalent to `P => :p`
2. Multiple types surrounded by `{}`:
	1. a sequence of `Pair{DataType, Symbol}` `{T => :t, P => :p, Q => :q}`
	2. a sequence of type `{T, P, Q}`: In this case, like above it is required that the struct `W`, for every type in the tuple, has a single field of that type the it is possible to have an unique mapping.
	3. multiples of the same type in the sequence: `{T, T, T}`, so long as there are exactly that many fields in the struct with that type, otherwise the pair syntax has to be used.

Here `<M>` is an (optional) tuple with two possible different elements

1. method names: `sin`, `cos`, `exp`, ... to extend a specific set of methods
2. module names: `Base`, `LinearAlgebra`, ... to extend all matching methods withing the specified modules.
3. if not specified all matching methods in the current module are extended

### Single Type Forwarding Examples:
```julia
@forward P struct W <: MaybeSubtyped
	p::P 	#= p is the only field of type P =#
	x 	#= untyped fields are simply ignored =#
  	...
end (method1, method2, method3)
```
We can also have multiple `P` in `W` but then we need the extended syntax:
then, it is clear over which element we want to expand
```julia
@forward P => :p1 struct W <: MaybeSubtyped
	p1::P
	p2::P
	...
end (Module)
```

or a more concrete example:
```jl
abstract type AbstractPolygon end

mutable struct Polygon <: AbstractPolygon
    x::Vector{Float64}
    y::Vector{Float64}
end

# Retrieve the number of vertices, and their X and Y coordinates
vertices(p::Polygon) = length(p.x)
coords_x(p::Polygon) = p.x
coords_y(p::Polygon) = p.y

# Move the polygon
function move!(p::Polygon, dx::Real, dy::Real)
    p.x .+= dx
    p.y .+= dy
end

@forward Polygon,
mutable struct RegularPolygon <: AbstractPolygon
    p::Polygon
    radius::Float64

	function RegularPolygon(n::Integer, radius::Real)
		@assert n >= 3
		θ = range(0, stop=2pi - (2pi / n), length=n)
		c = radius .* exp.(im .* θ)
		return RegularPolygon(Polygon(real(c), imag(c)), radius)
	end
end


square = RegularPolygon(4, 5)

# we have all the methods for Polygon automatically available also for RegularPolygon
@assert vertices(square) == 4
```

### Splat Forwarding:
same thing with `<T>` being a sequence of types in braces
```julia
@forward {T, P, Q} struct W <: MaybeSubtyped
	t::T
	...
	p::P
	...
	q::Q
	...
end
```

Then the macro will define the struct and then generate forwarding methods
for the methods that have the specified argument patterns,
consuming the patterns from left to right
(in each pattern the ".." represent arguments WITHOUT the pattern `<T>`)
```
#= method signature =#                 #= generated methods =#
m(.., <T>, ..)          --> m(.., W, ..) --> <T> = P:         m(.., W.p, ..)
                                         |-> <T> = P => p2:   m(.., W.p2, ..)
                                         |-> <T> = {T, P, Q}: m(.., W.t, W.p, W.q, ..)

m(.., <T>, .., <T>, ..) --> m(.., <T>, .., W, ..)
                        |-> m(.., W, .., <T>, ..)
                        |-> m(.., W, .., W, ..)
```
Assume the pattern `<T> = {T, T}` and that there is a method `m` with a signature with `m(.., T, T, T, ..)`.
Then we will only generate the methods:

* `m(.., W, T, ..)`
* `m(.., T, W, ..)`.

while methods with a `m(.., T, .., T, ..)` will be ignored. This is because if we decouple the ordering
we will end up with an explosion of methods and avoiding conflicting methods is much harder.

#### Use case:
```julia
@forward {Int, Int} struct Point2
	x::Int
	y::Int
end

@forward {Int, Int, Int} struct Point3
	x::Int
	y::Int
	x::Int
end

method1(a::Int, b::Int)             --> method1(Point2) = method1(Point2.x, Point2.y)
method2(a::Int, b::Int, depth::Int) --> method2(Point2, depth) = method2(Point2.x, Point2.y, depth)
                                    |-> method2(a, Point2) = method2(a, Point2.x, Point2.y)
                                    |-> method2(Point3) = method2(Point3.x, Point3.y, Point3.z)
```
More complex and arbitrary calls such as `method2(Point2.x, b, Point2.y)` can be defined manually on a case by case basis.
Defining a rule for such arbitrary cases would require making assumptions on the signature of any method, which in principle can have an arbitrary signature, and it is not really possible to enforce it.

## Multiple Forwarding
I've breifly explored the idea of having simultaneously forward a type as two different types. But even
simple cases like:
```julia
#= we want W to forward both as an Int and as a Float64 =#
@forward [Int, Float64] struct W
	p::Int
	q::Float64
end
```
But:
```
m1(Int, Int) -> m1(W, W) -----|
                              |--> m1(W.p, W.p) or m1(W.p, W.q)? Undecidable.
m1(Int, Float64) -> m1(W, W) -|
```

<!-- # ==[Multiple Forwarding]== (FOR THE FUTURE.)
# The wrapper can be made to cover two different objects at the same time.
# Opposed to the single forwarding, when defining multiple forwarding needs to be
# much more strict on the assumptions to avoid conflicting methods.
#
# 		 @derive (<T1>, <T2>, ...) struct W <: MaybeSubtyped
# 		 	...
# 		 	fields
# 		 	...
# 		 end (<M>)
#
# Here <T1> and <T2> rappresents two argument patterns, with the same rules explained
# for the single forwarding. But we will need to introduce some extra constraints on
# all patterns <T*> to avoid method collision as well as on the methods generated.
#
# Pattern constraint: No pattern can strictly contain any other pattern, we still allow intersections.
#   that means that we can allow two patterns such as {T,P} and {P,Q} but we can't allow {T,P} and P at the same time.
#   That is ok, since if we have disjoint patterns in a method arguments such as
#  		m(.., T, P, .., P, Q, ..) we can univocally map them one by one.
#   or if we have joint patterns
#  		m(.., T, P, Q, ..) it can be mapped to m(.., W, Q, ..) and m(.., T, W, ..) while having clear
#   what both these expressions means.
#   In the latter case the rule would disallow {T,P,Q} to coexsists with either {T,P} or {P,Q}
#   due to it raising uncertainty in how to interpret the above pattern.

# Method constraint:
#	We can only forward one argument at a time for each matching pattern.
#	For example:
#
# 		 @derive (P, Q) struct W
# 		   p::P
# 		   q::Q
# 		 end

# m(.., P, ..) -> m(.., W, ..) = m(.., W.p, ..)
# m(.., Q, ..) -> m(.., W, ..) = m(.., W.q, ..)
#
# m(.., P, .., P, ..) -> m(.., P, .., W, ..) = m(.., P, .., W.p, ..)
# 					  -> m(.., W, .., P, ..) = m(.., W.p, .., P, ..)
#					    ** Notice how we are not defining all values: why?
#					    ** because now W has to act both as P and Q
# m(.., P, .., Q, ..) -> m(.., P, .., W, ..) = m(.., P, .., W.q, ..)
# 				 	  -> m(.., W, .., Q, ..) = m(.., W.p, .., Q, ..)
#
# m(.., P, .., Q, .., P, ..) -> m(.., W, .., Q, .., P) = m(.., W.p, .., Q, .., P)
# 				 	  		 -> m(.., P, .., W, .., P) = m(.., P, .., W.q, .., P)
# 				 	  		 -> m(.., P, .., Q, .., W) = m(.., P, .., Q, .., W.p)

# The wrapper type can't behave like P and Q at the same time. But it can behave like P
# or like Q one at a time.
# @derive (Int, Float64) struct W
#	p::Int
#	q::Float64
# end
#
# m1(Int, Int) -> m1(W, W) -----|
#								|--> m1(W.p, W.p) or m1(W.p, W.q)? Undecidable.
# m1(Int, Float64) -> m1(W, W) -| -->