module Derive
using InteractiveUtils
using MacroTools
using Combinatorics

export @derive

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

@noinline panic(msg) = throw(ArgumentError(msg))
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
        if e isa Symbol
            e in fieldtypes || type_not_in_struct(e)
            type_count[e] = get(type_count, e, 0) + 1
        end
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
        ex = T[i]

        if ex isa Symbol
            S = ex
            type_indexes[S] = findnext(isequal(S), fieldtypes, type_indexes[S]) + 1
            f = fieldnames[type_indexes[S]-1]
        elseif isvalid_pair(ex)
            S = ex.args[2]
            f = ex.args[3].value
        end
        S in fieldtypes || type_not_in_struct(S)
        f in fieldnames || field_not_in_struct(f)
        return Expr(:call, :(=>), S, QuoteNode(f))

    end

    return expanded_pairs
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

    Stype = @capture(Sdef, S_t_ <: S_sup_) ? S_t : Sdef

    fieldtypes = []
    fieldnames = []
    for expr in Sfields
        @capture(expr, fn_::ft_ | fn_) || error("unsupported field.")

        ft = isnothing(ft) ? :Any : ft
        push!(fieldtypes, ft)
        push!(fieldnames, fn)
    end

    derive_pairs = expand_to_pairs(parse_braces(T), fieldnames, fieldtypes)
    evaldpairs = [Core.eval(_module_, d) for d in derive_pairs]
    sig = first.(evaldpairs)

    # builds a set of methods that contain our signature:
    # 1. first, get a set of methods that contain at least all our types singularly
    # !! We limit the search to the current module for now
    candidate_methods = intersect(Set.(methodswith.(sig, (_module_,)))...)

    # 2. filter all the methods that have less than the number of types in our signature
    filter!(m -> m.nargs > length(sig), candidate_methods)

    # 3. remove constructors from the methods to derive
    filter!(m -> getfield(_module_, m.name) isa Function, candidate_methods)

    # 3. Scan the signature and discard all methods that do not contain the correct
    # order. At the same time, match the matching positions with the method.
    allmethods = Dict()
    for m in candidate_methods
        msig = fieldtype.(m.sig, collect(2:m.nargs))
        sigwidth = length(sig) - 1

        positions = []
        for i in 1:length(msig)-sigwidth
            if all(typeintersect.(msig[i:i+sigwidth], sig) .!= Union{})
                push!(positions, i:i+sigwidth)
            end
        end
        if !isempty(positions)
            allmethods[m] = positions
        end
    end

    # Use: Base.arg_decl_parts(method) to deconstruct the method signature
    # such that to keep the typevar information: x<:T y<:T where T<:...

    # generate the new signatures
    methods_to_generate = []
    for (m, swap_positions) in allmethods
        tv, decl, _... = Base.arg_decl_parts(m)
        argnames = [Symbol(d[1]) for d in decl[2:end]]
        argtypes = [Symbol(d[2]) for d in decl[2:end]]

        for positions in combinations(swap_positions)
            # newsiglen = (length(decl) - 1) - sum(length.(positions) .- 1)

            argnameswaps = [gensym(Stype) for _ in 1:length(positions)]
            argtypesswaps = fill(Stype, length(positions))

            newargnames = swapat(argnames, positions, argnameswaps)
            newargtypes = swapat(argtypes, positions, argtypesswaps)

            newdecl = collect(zip(newargnames, newargtypes))

            newsignature = generate_signature(m, newdecl)
            if !isempty(tv)
                # filter the tv to remove the typevars we substituted
                newtv = filter(tvar -> tvar.name in newargtypes, tv)
                newsignature = Expr(:where, newsignature, newtv...)
            end
            methodforwardcall = generate_forward_call(m, Stype, newdecl, evaldpairs)

            ex = Expr(:(=), newsignature, methodforwardcall)
            # @info ex
            push!(methods_to_generate, ex)
        end
    end

    retblk = Expr(:block)
    push!(retblk.args, S)
    for gm in methods_to_generate
        push!(retblk.args, gm)
    end

    @show esc(retblk)
    return esc(retblk)
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

function generate_signature(method, decl)
    mexpr = Expr(:call, method.name)
    for (argn, argt) in decl
        push!(mexpr.args, Expr(:(::), argn, argt))
    end
    return mexpr
end

function generate_forward_call(method, deriveT, decl, derivepairs)
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

end # module Derive
