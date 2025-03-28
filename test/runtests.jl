using Derive
using InteractiveUtils
using Test

# =-=-= ~~~~~~~~~~~ Setup ~~~~~~~~~~~~

abstract type AbstractPolygon end

mutable struct Polygon <: AbstractPolygon
    x::Vector{Float64}
    y::Vector{Float64}
end

# Retrieve the number of vertices, and their X and Y coordinates
vertices(p::Polygon) = length(p.x)
coords_x(p::Polygon) = p.x
coords_y(p::Polygon) = p.y

# Move, scale and rotate a polygon
function move!(p::Polygon, dx::Real, dy::Real)
    p.x .+= dx
    p.y .+= dy
end

function scale!(p::Polygon, scale::Real)
    m = mean(p.x)
    p.x = (p.x .- m) .* scale .+ m
    m = mean(p.y)
    p.y = (p.y .- m) .* scale .+ m
end

function rotate!(p::Polygon, angle_deg::Real)
    θ = float(angle_deg) * pi / 180
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    x = p.x .- mean(p.x)
    y = p.y .- mean(p.y)
    (x, y) = R * [x, y]
    p.x = x .+ mean(p.x)
    p.y = y .+ mean(p.y)
end

# =-=-= ~~~~~~~~~~ Extension

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


# Forward methods from `RegularPolygon` to `Polygon`
vertices(p1::RegularPolygon) = vertices(getfield(p1, :polygon))
coords_x(p1::RegularPolygon) = coords_x(getfield(p1, :polygon))
coords_y(p1::RegularPolygon) = coords_y(getfield(p1, :polygon))
move!(p1::RegularPolygon, p2::Real, p3::Real) = move!(getfield(p1, :polygon), p2, p3)
rotate!(p1::RegularPolygon, p2::Real) = rotate!(getfield(p1, :polygon), p2)
function scale!(p::RegularPolygon, scale::Real)
    scale!(p.polygon, scale) # call "super" method
    p.radius *= scale        # update internal state
end


# =-=-= Testing =-=-=
# create dummy types
struct A end
struct B end
struct T end
struct S end
struct P end
struct Q end

# test parsing
@testset "Brace Parsing" begin
    patterns = [
        :Q,
        :(P => :p),
        :({A, B}),
        :({T, T}),
        :({P => :p, A => :a, B}),
        :({P => :p, T, T})
    ]
    parse_result = [
        (:Q,),
        (:(P => :p),),
        (:A, :B),
        (:T, :T),
        (:(P => :p), :(A => :a), :B),
        (:(P => :p), :T, :T)
    ]
    bad_patterns = [
        :([Q, A, B]), # vec instead of bracers
        :((Q, A, B) => :a), # bad pair
        :(P => A) # bad pair
    ]
    @assert length(patterns) == length(parse_result)
    @testset for (i, pattern) in enumerate(patterns)
        @test Derive.parse_braces(pattern) == parse_result[i]
    end

    for i = 1:length(bad_patterns)
        @test_throws ArgumentError Derive.parse_braces(bad_patterns[i])
    end
end
# test expanding
@testset "Expand Types" begin
    parse_result = [
        (:Q,),
        (:(P => :p),),
        (:A, :B),
        (:T, :T),
        (:(P => :p), :(A => :a), :B),
        (:(P => :p), :T, :T)
    ]
    struct Wrap
        a::A
        b::B
        t1::T
        t2::T
        q::Q
        p::P
    end
    fnames = fieldnames(Wrap)
    ftypes = Symbol.(fieldtypes(Wrap))

    expand_results = [
        (:(Q => :q),),
        (:(P => :p),),
        (:(A => :a), :(B => :b)),
        (:(T => :t1), :(T => :t2)),
        (:(P => :p), :(A => :a), :(B => :b)),
        (:(P => :p), :(T => :t1), :(T => :t2))
    ]

    bad_parses = [
        (:T, :T, :T), # mismatching number of implicit fields
        (:(Noexist => :a),),# non existant field type in pair
        (:(P => :noexists),),# non existant field name in pair
        (:Noexists,) # non existant field type as symbol
    ]

    @assert length(expand_results) == length(parse_result)
    @testset for (i, pattern) = enumerate(parse_result)
        @test Derive.expand_to_pairs(pattern, fnames, ftypes) == expand_results[i]
    end

    @testset for badparse in bad_parses
        @test_throws ArgumentError Derive.expand_to_pairs(badparse, fnames, ftypes)
    end
end


