module AutoForward
using InteractiveUtils
using MacroTools
using Combinatorics

export @forward

# possible inputs:
# 	:Symbol ( must be automatically turned into a Pair)
#	:Sym => :sym
#	{:Sym, ...}
# 	{:Sym => :sym, ...}
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
        S == :Any && panic("Can't forward fields of type `Any`.")
        f in fieldnames || field_not_in_struct(f)
        return Expr(:call, :(=>), S, QuoteNode(f))

    end

    return expanded_pairs
end

macro forward(exargs...)
    if length(exargs) == 1 && isexpr(exargs[1], :tuple)
        TSM = exargs[1].args
        if length(TSM) == 3
            T, S, M = TSM
        elseif length(TSM) == 2
            T, S = TSM
            M = Expr(:tuple, Symbol(__module__))
        end
    elseif length(exargs) == 2
        TS, SM = exargs
        if isexpr(TS, :tuple)
            T, S = TS.args
            M = SM
        elseif isexpr(SM, :tuple)
            T = TS
            S, M = SM.args
        else
            T, S = TS, SM
            M = Expr(:tuple, Symbol(__module__))
        end
    else
        T, S, M = exargs
    end
    if M isa Symbol
        M = Expr(:tuple, M)
    end

    # @info "Z" T S M
    forward(__module__, T, S, M)
end

function forward(_module_, @nospecialize(T), @nospecialize(S), @nospecialize(M))

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

    @capture(M, (filters__,)) || panic("Specific functions or modules must be in a tuple.")
    materialized_filters = []
    for filter in filters
        if isexpr(filter)
            @capture(filter, modsym_.func_) || panic("unsupported filter: $filter")

            mod = isdefined(_module_, modsym) ? getfield(_module_, modsym) :
                  throw(UndefVarError(modsym, _module_))

            f = isdefined(mod, func) ? getfield(mod, func) :
                throw(UndefVarError(func, mod))

            f isa Function ? push!(materialized_filters, f) : panic("$f must be a Function")

        elseif filter isa Symbol
            mf = isdefined(_module_, filter) ? getfield(_module_, filter) :
                 throw(UndefVarError(filter, _module_))

            mf isa Union{Module,Function} ?
            push!(materialized_filters, mf) :
            panic("$mf must be a Module or a Function")
        else
            panic("unsupported filter: $filter")
        end
    end

    derive_pairs = expand_to_pairs(parse_braces(T), fieldnames, fieldtypes)
    evaldpairs = [Core.eval(_module_, d) for d in derive_pairs]
    sig = first.(evaldpairs)

    # builds a set of methods that contain our signature:
    # 1. first, get a set of methods that contain at least all our types singularly
    candidate_methods = Set()
    for mod_or_func in materialized_filters
        union!(candidate_methods, intersect(Set.(methodswith.(sig, (mod_or_func,)))...))
    end

    # 2. filter all the methods that have less than the number of types in our signature
    filter!(m -> m.nargs > length(sig), candidate_methods)

    # 3. remove constructors from the methods to derive
    filter!(m -> getfield(m.module, m.name) isa Function, candidate_methods)

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
            allmethods[m] = positions # we rely on the sortedness later so better ensure it
        end
    end

    # Use: Base.arg_decl_parts(method) to deconstruct the method signature
    # such that to keep the typevar information: x<:T y<:T where T<:...

    # generate the new signatures
    methods_to_generate = []
    for (m, swap_positions) in allmethods
        tv, decl, _... = Base.arg_decl_parts(m) # internal
        argnames = [Symbol(d[1]) for d in decl[2:end]]
        argtypes = [Symbol(d[2]) for d in decl[2:end]]

        for positions in combinations(swap_positions)
            ranges_overlap_pairwise(sort!(positions)) && continue

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
            push!(methods_to_generate, ex)
        end
    end

    retblk = Expr(:block)
    push!(retblk.args, S)
    for gm in methods_to_generate
        push!(retblk.args, gm)
    end

    # @show esc(retblk)
    return esc(retblk)
end

function swapat(base, positions, swaps)
    @assert length(positions) == length(swaps)

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

function generate_forward_call(method, forward_t, decl, derivepairs)
    mname = :($(method.module).$(method.name))
    mexpr = Expr(:call, mname)
    fields = last.(derivepairs)
    for (argn, argt) in decl
        if argt == forward_t
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

function ranges_overlap_pairwise(positions)
    any(do_overlap(positions[i], positions[i+1]) for i = 1:length(positions)-1)
end
do_overlap(range1, range2) = max(range1[begin], range2[begin]) <= min(range1[end], range2[end])

end # module AutoForward
