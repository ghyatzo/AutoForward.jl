module MethodForwarding
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
function isvalid_type(e::Expr)
    isexpr(e, :curly) || isexpr(e, :where)
end
isvalid_type(s::Symbol) = true

function isvalid_pair(e::Expr)
    Base.isexpr(e, :call) &&
        e.args[1] == :(=>) &&
        isvalid_type(e.args[2]) &&
        e.args[3] isa QuoteNode
end


@noinline panic(msg) = throw(ArgumentError(msg))
type_not_in_struct(type) = panic("The type $type is not present in the struct.")
field_not_in_struct(field) = panic("The field $field is not present in the struct.")

function parse_braces(T)
    if isvalid_type(T) || isvalid_pair(T)
        return (T,)
    end

    if Base.isexpr(T, :braces)
        return ntuple(length(T.args)) do i
            T.args[i]
        end
    end
    display(T)
    throw(ArgumentError("Invalid pattern."))
end

function expand_to_pairs(T, fieldnames, fieldtypes)
    @info T fieldnames fieldtypes
    # check for ambiguity
    type_count = Dict{Any,Int}()
    for e in T
        if isvalid_type(e)
            key = isexpr(e, :where) ? e.args[1] : e
            key in fieldtypes || type_not_in_struct(key)
            type_count[key] = get(type_count, key, 0) + 1
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
        @show ex
        if isvalid_type(ex)
            S = ex
            key = isexpr(ex, :where) ? ex.args[1] : ex
            type_indexes[key] = findnext(isequal(key), fieldtypes, type_indexes[key]) + 1
            f = fieldnames[type_indexes[key]-1]
        elseif isvalid_pair(ex)
            S = ex.args[2]
            notunionall_S = isexpr(S, :where) ? S.args[1] : S
            f = ex.args[3].value
            notunionall_S in fieldtypes || type_not_in_struct(S)
        end
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
    @info derive_pairs
    evaldpairs = [Core.eval(_module_, d) for d in derive_pairs]
    @info evaldpairs
    sig = first.(evaldpairs)
    if any(s == Any for s in sig)
        panic("Can't forward over Any.")
    end

    # builds a set of methods that contain our signature:
    # 1. first, get a set of methods that contain at least all our types singularly
    candidate_methods = Set()
    for mod_or_func in materialized_filters
        union!(candidate_methods, intersect(Set.(methodswith.(sig, (mod_or_func,)))...))
    end

    filter!(m -> begin
            m.nargs > length(sig) &&                # has enough arguments
                !startswith(string(m.name), '@') &&  # is not a macro
                fieldtype(m.sig, 1) <: Function     # is not a constructor
        end, candidate_methods)
    # Scan the signature and discard all methods that do not contain the correct
    # order. At the same time, match the matching positions with the method.
    allmethods = Dict()
    for m in candidate_methods
        msig = fieldtype.(m.sig, collect(2:m.nargs))
        sigwidth = length(sig) - 1

        positions = []
        for i in 1:length(msig)-sigwidth
            msig_window = msig[i:i+sigwidth]
            if all(sig .<: msig_window) && !any(msig_window .== Any)
                push!(positions, i:i+sigwidth)
            end
        end

        if !isempty(positions)
            allmethods[m] = positions # we rely on the sortedness later so better ensure it
        end
    end

    ## Useful decomposition:
    # sig is a unionall T{<:P, S} where {S}
    # unwrap the unionall which is nested:
    # ua = unwrap_unionall(sig)
    # T = nameof(ua)
    # typevars: sig.var
    # parameters = ua.parameters (which are typevars type with lowerbound, name and upperbound)

    # generate the new signatures
    argnametag = isexpr(Stype, :curly) ? Stype.args[1] : Stype
    methods_to_generate = []
    for (m, swap_positions) in allmethods
        msig = m.sig
        tv = []
        while msig isa UnionAll
            push!(tv, msig.var)
            msig = msig.body
        end
        # tv, decl, _... = Base.arg_decl_parts(m)

        argnames = Base.method_argnames(m)[2:end]
        argtypes = fieldtype.(m.sig, collect(2:m.nargs))
        for positions in combinations(swap_positions)
            ranges_overlap_pairwise(sort!(positions)) && continue

            argnameswaps = [gensym(argnametag) for _ in 1:length(positions)]
            argtypesswaps = fill(Stype, length(positions))

            newargnames = swapat(argnames, positions, argnameswaps)
            newargtypes = swapat(argtypes, positions, argtypesswaps)

            newdecl = collect(zip(newargnames, newargtypes))

            if !isempty(tv)
                # filter the tv to remove the typevars we substituted
                newtv = filter(tvar -> tvar.name in newargtypes, tv)
            else
                newtv = []
            end

            # we need to add the tv variants for the
            # remaining parameters of the struct
            if isexpr(Stype, :curly)
                params = Stype.args[2:end]
                #TODO: multiple signatures
                if sig[1] isa UnionAll
                    sigtvs = [p.name for p in Base.unwrap_unionall(sig[1]).parameters]
                else
                    sigtvs = []
                end
                paramtv = [TypeVar(p) for p in params if p ∉ sigtvs]
                push!(newtv, paramtv...)
            end
            # for each unionall in the signature we desire
            # extract and save the typevariables
            for typ in sig
                typ isa UnionAll || continue
                typ = Base.unwrap_unionall(typ)
                tvs = typ.parameters
                for tv in tvs
                    tv isa TypeVar && push!(newtv, tv)
                end
            end
            newsignature = generate_signature(m, newdecl, newtv)
            @info m newsignature
            # Meta.dump(newsignature)
            methodforwardcall = generate_forward_call(m, Stype, newdecl, evaldpairs)

            methodforwardcall == Symbol("#skip#") && continue

            ex = Expr(:(=), newsignature, methodforwardcall)
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
    mname = :($(method.module).$(method.name))
    mexpr = Expr(:call, mname)
    for (argn, argt) in decl
        ex = argn == Symbol("#unused#") ?
             Expr(:(::), argt) : Expr(:(::), argn, argt)

        push!(mexpr.args, ex)
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
            if argn == Symbol("#unused#")
                defaultconstructor = methods(argt, Tuple{})
                isempty(defaultconstructor) && return Symbol("#skip#")
                push!(mexpr.args, Expr(:call, argt))
            else
                push!(mexpr.args, argn)
            end
        end
    end
    return mexpr
end

function ranges_overlap_pairwise(positions)
    any(do_overlap(positions[i], positions[i+1]) for i = 1:length(positions)-1)
end
do_overlap(range1, range2) = max(range1[begin], range2[begin]) <= min(range1[end], range2[end])

end # module MethodForwarding
