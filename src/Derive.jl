module Derive
using InteractiveUtils
using MacroTools
using Combinatorics

export AbstractPolygon, Polygon, RegularPolygon, vertices, coords_x, coords_y,
    side, area, @derive

abstract type AbstractPolygon end

mutable struct Polygon <: AbstractPolygon
    x::Vector{Float64}
    y::Vector{Float64}
end

# Retrieve the number of vertices, and their X and Y coordinates
vertices(p::Polygon) = length(p.x)
coords_x(p::Polygon) = p.x
coords_y(p::Polygon) = p.y

mutable struct RegularPolygon <: AbstractPolygon
    p::Polygon
    radius::Float64
end

function RegularPolygon(n::Integer, radius::Real)
    @assert n >= 3
    θ = range(0, stop=2pi - (2pi / n), length=n)
    c = radius .* exp.(im .* θ)
    return RegularPolygon(Polygon(real(c), imag(c)), radius)
end

# Extended methods only applicable to a Regular Polygon
# Compute length of a side and the polygon area
side(p::RegularPolygon) = 2 * p.radius * sin(pi / vertices(p))
area(p::RegularPolygon) = side(p)^2 * vertices(p) / 4 / tan(pi / vertices(p))

### EXTENSIONS TO AUTOMATE
vertices(p::RegularPolygon) = vertices(p.p)
coords_x(p::RegularPolygon) = coords_x(p.p)
coords_y(p::RegularPolygon) = coords_y(p.p)
###

# ---- ADDING ANOTHER LAYER!!

struct StrangeName
    name::String
end
name(s::StrangeName) = "Spooooky: $s"

mutable struct StrangePolygon <: AbstractPolygon
    P::RegularPolygon
    name::StrangeName
end

### METHODS TO EXTEND!
vertices(p::StrangePolygon) = vertices(p.p)
coords_x(p::StrangePolygon) = coords_x(p.p)
coords_y(p::StrangePolygon) = coords_y(p.p)

side(p::StrangePolygon) = side(p.p)
area(p::StrangePolygon) = area(p.p)

name(s::StrangePolygon) = name(s.name)

## Useful stuff
# methodswith(type, [Module], supertype=false) -> returns a list of methods with a specified type
# names()
# parentmodule()
# fieldtypes(T::Type)
# fieldnames(x::DataType)
# nameof(t::DataType) -> Symbol

## Generalised API

# ==[Single Forwarding]==
# the wrapper will transparently behave as the specified type (if a single symbol)
# or as a tuple of the the specified types that automatically splats:
#
# 		 @derive <T> struct W <: MaybeSubtyped
# 		 	...
# 		 	fields
# 		 	...
# 		 end (<M>)
#
# Here <T> can be different things:
#	1. a Pair{DataType, Symbol}:  P => :p
#	1b. a Symbol: P
# 		** In this case it is required that the struct W has a single P type **
#		** then P can be univocally made equivalent to P => :p **
#	2. a sequence of Pair{DataType, Symbol} surrounded by braces {}: {T => :t, P => :p, Q => :q}
#	2b. a sequence of types surrounded by braces {}: {T, P, Q}
#		** In this case, like above it is required that the struct W, **
#		** for every type in the tuple, has a single field of that type **
#		** the it is possible to have an unique mapping **
#	2c. you can have multiple of the same type in the sequence: {T, T, T}, so long as
#	 	there are exactly that many fields in the struct with that type, in any other case the
#		pair syntax has to be used.
#
# Here <M> is an (optional) tuple with two possible different elements
#	1. method names: sin, cos, exp, ... to extend a specific set of methods
#	2. module names: Base, LinearAlgebra, ... to extend all matching methods withing the specified modules.
#	3. if not specified all matching methods are extended

# Examples:
#
#        @derive P struct W <: MaybeSubtyped
#           p::P 	** p is the only field of type P **
#           x 		** untyped fields are simply ignored **
#           ...
#        end (method1, method2, method3)

# We can also have multiple P in W but then we need the extended syntax:
# then, it is clear over which element we want to expand

#        @derive P => :p1 struct W <: MaybeSubtyped
#        	p1::P
#        	p2::P
#        	...
#        end (method1, method2, method3)

# same thing with <T> being of type (2)

#        @derive {T, P, Q} struct W <: MaybeSubtyped
#        	t::T
#       	..
#       	p::P
#       	..
#       	q::Q
#       	..
#        end (method1, method2, method3)


# Then the macro will define the struct and then generate forwarding methods
# for the methods that have the specified argument patterns,
# consuming the patterns from left to right
# (in each pattern the ".." represent arguments WITHOUT the pattern <T>)
# m(.., <T>, ..) 		  -> m(.., W, ..) = --> <T> = P: 	  	   m(.., W.p, ..)
#				  						    |-> <T> = P => p2:   m(.., W.p2, ..)
#										    |-> <T> = (T, P, Q): m(.., W.t, W.p, W.q, ..)

# m(.., <T>, .., <T>, ..) -> m(.., <T>, .., W, ..)
# 					  	  -> m(.., W, .., <T>, ..)
#					  	  -> m(.., W, .., W, ..)

# Assume the pattern <T> = {T, T}
# and there is a method with a m(.., T, T, T, ..) in it's signature:
# then we will only generate m(.., W, T, ..) and m(.., T, W, ..) methods.
# while methods with a m(.., T, .., T, ..) will be ignored. This is because if we decouple the ordering
# we will end up with an explosion of methods and avoiding conflicting methods is much harder.
#
# Use case:
# @derive {Int, Int} struct Point2 		** notice the double parenthesis
# 	x::Int
# 	y::Int
# end
# @derive {Int, Int, Int} struct Point3	** notice the double parenthesis
# 	x::Int
# 	y::Int
# 	x::Int
# end

# method1(a::Int, b::Int)			  -> method1(Point2)
# method2(a::Int, b::Int, depth::Int) -> method2(Point2, depth)
# 									  -> method2(a, Point2)
#									  -> method2(Point3)

# More complex and arbitrary calls such as method2(Point2.x, b, Point2.y)
# can be defined manually on a case by case basis.
# Defining a rule for such arbitrary cases would require making assumptions on the
# signature of any method, which in principle can have arbitrary signature, and
# it is not really possible to enforce it.

# ==[Multiple Forwarding]== (FOR THE FUTURE.)
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
# m1(Int, Float64) -> m1(W, W) -|


# ASSUMPTION (FOR NOW)
# All types to be derived must be present in the struct.
# All types to be derived must be unique within the struct. No two fields with the same derived type.

# possible inputs:
# 	:Symbol ( must be automatically turned into a Pair)
#	:Sym => :sym
#	{:Sym, ...}
# 	{:Sym => :sym, ...}
#	tuple of all of the above...
# output:
#	Tuple of :Sym => :sym
function isvalid_pair(e::Expr)
    Base.isexpr(e, :call) &&
        e.args[1] == :(=>) &&
        e.args[2] isa Symbol &&
        e.args[3] isa QuoteNode
end

panic(msg) = throw(ArgumentError(msg))
type_not_in_struct(type) = panic("The type $type is not present in the struct.")
field_not_in_struct(field) = panic("The field $field is not present in the struct.")

"""
    parse_braces(T)

This helper function is to allow the `{...}` syntax for the macro. It translates any brace
into a tuple, recursively. Additionally it checks that all inputs are of valid type.
You should not call this method manually.

T can be a `:Symbol`, a `{ ... }` or a `Pair{Type, Symbol}` (`MyType => :field`)

"""
function parse_braces(T)
    if T isa Symbol || isvalid_pair(T)
        return (T,)
    end

    if Base.isexpr(T, :braces)
        return ntuple(length(T.args)) do i
            T.args[i]
        end
    end

    throw(ArgumentError("Invalid pattern."))
end

"""
    expand_types(T, fieldnames, fieldtypes)

Helper function to perform the automatic expansion of types where there is no ambiguity.
For example, a `struct` with a single field `p::P` and a deriving type defined as `@derive P`
then we can expand `P` automatically into `P => :p`.
Checks validity of the inputs.

`T` is a tuple obtained from the `parse_braces` helper function.
`fieldnames` is a vector of the fieldnames of the struct.
`fieldtypes` is a vector of Symbols for the MATCHING field types of the struct.

"""
function expand_to_pairs(T, fieldnames, fieldtypes)

    # check for ambiguity
    type_count = Dict{Symbol,Int}()
    for e in T
        e isa Symbol || continue
        e in fieldtypes || type_not_in_struct(e)
        type_count[e] = get(type_count, e, 0) + 1
    end

    for k in keys(type_count)
        if type_count[k] != count(isequal(k), fieldtypes)
            throw(ArgumentError("""
            	Mismatch between number of implicit `$k` to be derived compared to explicit fields of type `$k` in the struct.
            	Specify the derive types by using the explicit notation: $k => <Symbol of the field name>
            	or match the number of fields in the struct.
            """))
        end
    end

    # match fields
    type_indexes = Dict(k => 1 for k in keys(type_count))
    expanded_pairs = ntuple(length(T)) do i
        e = T[i]

        if e isa Symbol
            e in fieldtypes || type_not_in_struct(e)
            # match fields in order of appearance in the struct
            type_indexes[e] = findnext(isequal(e), fieldtypes, type_indexes[e]) + 1
            matching_field = fieldnames[type_indexes[e]-1]
            return Expr(:call, :(=>), e, QuoteNode(matching_field))
        elseif isvalid_pair(e)
            (e.args[2] in fieldtypes) || type_not_in_struct(e.args[2])
            (e.args[3].value in fieldnames) || field_not_in_struct(e.args[3])
            return e
        end

    end

    expanded_pairs
end

macro derive(args...)
    derive(__module__, args...)
end

function derive(_module_, T, S)

    @capture(S, struct Sdef_
        Sfields__
    end | mutable struct Sdef_
        Sfields__
    end) || error("Needs to be a struct.")

    if @capture(Sdef, S_t_ <: S_sup_)
        Stype = S_t
    else
        Stype = Sdef
    end

    @show Sfields
    fieldtypes = []
    fieldnames = []
    for expr in Sfields
        @capture(expr, fn_::ft_ | fn_) || error("unsupported field.")

        ft = isnothing(ft) ? :Any : ft
        push!(fieldtypes, ft)
        push!(fieldnames, fn)
    end
    @show fieldtypes
    @show fieldnames
    derive_pairs = expand_to_pairs(parse_braces(T), fieldnames, fieldtypes)
    @show derive_pairs
    evaldpairs = [Core.eval(_module_, d) for d in derive_pairs]
    @info evaldpairs
    sig = first.(evaldpairs)
    @show sig

    # builds a set of methods that contain our signature:
    # 1. first, get a set of methods that contain at least all our types singularly
    meths_1 = intersect(Set.(methodswith.(sig))...)

    # 2. filter all the methods that have less than the number of types in our signature
    meths_2 = filter(m -> m.nargs > length(sig), meths_1)

    # 3. Scan the signature and discard all methods that do not contain the correct
    # order. At the same time, match the matching positions with the method.
    meths_3 = Dict()
    for m in meths_2
        msig = fieldtype.(m.sig, collect(2:m.nargs))
        sigwidth = length(sig) - 1

        positions = []
        for i in 1:length(msig)-sigwidth
            if all(typeintersect.(msig[i:i+sigwidth], sig) .!= Union{})
                push!(positions, i:i+sigwidth)
            end
        end

        meths_3[m] = positions
    end

    # Use: Base.arg_decl_parts(method) to deconstruct the method signature
    # such that to keep the typevar information: x<:T y<:T where T<:...

    # generate the new signatures
    genmethods = []
    for (m, swap_positions) in meths_3
        tv, decl, _... = Base.arg_decl_parts(m)
        argnames = [Symbol(d[1]) for d in decl[2:end]]
        argtypes = [Symbol(d[2]) for d in decl[2:end]]

        for positions in combinations(swap_positions)
            # newsiglen = (length(decl) - 1) - sum(length.(positions) .- 1)

            # TODO: deal with typevars
            argnameswaps = [gensym(Stype) for _ in 1:length(positions)]
            argtypesswaps = fill(Stype, length(positions))

            sendargnames = swapat(argnames, positions, argnameswaps)
            sendargtypes = swapat(argtypes, positions, argtypesswaps)

            senddecl = collect(zip(sendargnames, sendargtypes))

            # For each element of our new type, we expand into the pattern, following the derive_pairs
            # @info gennewsig(m, senddecl)
            # @info gensplatsig(m, Stype, senddecl, derive_pairs)
            ex = :($(gennewsig(m, senddecl)) = $(gensplatsig(m, Stype, senddecl, evaldpairs)))
            @info ex
            push!(genmethods, ex)
        end
    end

    retblk = Expr(:block)
    push!(retblk.args, S)
    for gm in genmethods
        push!(retblk.args, gm)
    end

    @show esc(retblk)
    # return esc(retblk)
end

function swapat(base, positions, swaps)
    @assert length(positions) == length(swaps)
    # positions must be non overllapping. TODO
    swapped = []
    for ip in eachindex(positions)
        if ip == 1
            startpos = firstindex(base)
            endpos = positions[ip][begin] - 1
        else
            startpos = positions[ip-1][end] + 1
            endpos = positions[ip][begin] - 1
        end
        push!(swapped, base[startpos:endpos]...)
        push!(swapped, swaps[ip])
    end

    if last(positions)[end] < length(base)
        startpos = last(positions)[end] + 1
        endpos = length(base)
        push!(swapped, base[startpos:endpos]...)
    end

    return swapped
end

function gennewsig(method, decl)
    mexpr = Expr(:call, method.name)
    for (argn, argt) in decl
        push!(mexpr.args, Expr(:(::), argn, argt))
    end
    return mexpr
end

function gensplatsig(method, deriveT, decl, derivepairs)
    mname = :($(method.module).$(method.name))
    mexpr = Expr(:call, mname)
    fields = last.(derivepairs)
    for (argn, argt) in decl
        if argt == deriveT
            for f in fields
                getfieldex = Expr(:call, :getfield, argn, QuoteNode(f))
                push!(mexpr.args, getfieldex)
            end
        else
            push!(mexpr.args, argn)
        end
    end
    return mexpr
end


# macro forward(ex, fs)
#     @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
#     T = esc(T)
#     fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
#     :($([:($f(x::$T, args...; kwargs...) = (Base.@_inline_meta; $f(x.$field, args...; kwargs...)))
#          for f in fs]...);
#     nothing)
# end

function forward(sender::Tuple{Type,Symbol}, receiver::Type, method::Method;
    withtypes=true, allargs=true)
    function newmethod(sender_type, sender_symb, argid, method, withtypes, allargs)
        s = "p" .* string.(1:method.nargs-1)
        (withtypes) && (s .*= "::" .* string.(fieldtype.(Ref(method.sig), 2:method.nargs)))
        s[argid] .= "p" .* string.(argid) .* "::$sender_type"
        if !allargs
            s = s[1:argid[end]]
            push!(s, "args...")
        end

        # Module where the method is defined
        ff = fieldtype(method.sig, 1)
        if isabstracttype(ff)
            # costructors
            m = string(method.module)
            # m = string(method.module.eval(:(parentmodule($(method.name)))))  # Constructor
        else
            # all methods except constructors
            m = string(parentmodule(ff))
        end
        m *= "."
        l = "$m:(" * string(method.name) * ")(" * join(s, ", ") * "; kw..."
        m = string(method.module) * "."
        l *= ") = $m:(" * string(method.name) * ")("
        s = "p" .* string.(1:method.nargs-1)
        if !allargs
            s = s[1:argid[end]]
            push!(s, "args...")
        end
        s[argid] .= "getfield(" .* s[argid] .* ", :$sender_symb)"
        l *= join(s, ", ") * "; kw...)"
        l = join(split(l, "#"))
        return l
    end

    @assert isstructtype(sender[1])
    @assert sender[2] in fieldnames(sender[1])
    sender_type = string(parentmodule(sender[1])) * "." * string(nameof(sender[1]))
    sender_symb = string(sender[2])
    code = Vector{String}()

    # Search for receiver type in method arguments
    foundat = Vector{Int}()
    for i in 2:method.nargs
        argtype = fieldtype(method.sig, i)
        (sender[1] == argtype) && (return code)
        if argtype != Any
            (typeintersect(receiver, argtype) != Union{}) && (push!(foundat, i - 1))
        end
    end
    (length(foundat) == 0) && (return code)
    if string(method.name)[1] == '@'
        @warn "Forwarding macros is not yet supported."  # TODO
        display(method)
        println()
        return code
    end

    for ii in combinations(foundat)
        push!(code, newmethod(sender_type, sender_symb, ii, method, withtypes, allargs))
    end

    tmp = split(string(method.module), ".")[1]
    code = "@eval " .* tmp .* " " .* code .*
           " # " .* string(method.file) .* ":" .* string(method.line)
    if (tmp != "Base") &&
       (tmp != "Main")
        pushfirst!(code, "using $tmp")
    end
    code = unique(code)
    return code
end


end # module Derive
