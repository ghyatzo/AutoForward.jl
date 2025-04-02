using Test
using MethodForwarding
using Statistics

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
        @test MethodForwarding.parse_braces(pattern) == parse_result[i]
    end

    for i in eachindex(bad_patterns)
        @test_throws ArgumentError MethodForwarding.parse_braces(bad_patterns[i])
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
        @test MethodForwarding.expand_to_pairs(pattern, fnames, ftypes) == expand_results[i]
    end

    @testset for badparse in bad_parses
        @test_throws ArgumentError MethodForwarding.expand_to_pairs(badparse, fnames, ftypes)
    end
end


# Automatic Derivations
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
mutable struct TestRegularPolygon <: AbstractPolygon
    p::Polygon
    radius::Float64
end

function TestRegularPolygon(n::Integer, radius::Real)
    @assert n >= 3
    θ = range(0, stop=2pi - (2pi / n), length=n)
    c = radius .* exp.(im .* θ)
    return TestRegularPolygon(Polygon(real(c), imag(c)), radius)
end

# Extended methods only applicable to a Regular Polygon
# Compute length of a side and the polygon area
side(p::TestRegularPolygon) = 2 * p.radius * sin(pi / vertices(p))
area(p::TestRegularPolygon) = side(p)^2 * vertices(p) / 4 / tan(pi / vertices(p))

# Forward methods from `RegularPolygon` to `Polygon`
# Manually implemented to check
vertices(p::TestRegularPolygon) = vertices(getfield(p, :p))
coords_x(p::TestRegularPolygon) = coords_x(getfield(p, :p))
coords_y(p::TestRegularPolygon) = coords_y(getfield(p, :p))
move!(p::TestRegularPolygon, p2::Real, p3::Real) = move!(getfield(p, :p), p2, p3)
rotate!(p::TestRegularPolygon, p2::Real) = rotate!(getfield(p, :p), p2)
function scale!(p::TestRegularPolygon, scale::Real)
    scale!(p.p, scale) # call "super" method
    p.radius *= scale        # update internal state
end

# Automatic derivation
@forward Polygon,
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

testregpoly = TestRegularPolygon(4, 5.0)
derivedregpoly = RegularPolygon(4, 5.0)

@test vertices(testregpoly) == vertices(derivedregpoly)
@test coords_x(testregpoly) == coords_x(derivedregpoly)
@test coords_y(testregpoly) == coords_y(derivedregpoly)
move!(testregpoly, 1, 1)
move!(derivedregpoly, 1, 1)
@test coords_x(testregpoly) == coords_x(derivedregpoly)
@test coords_y(testregpoly) == coords_y(derivedregpoly)
rotate!(testregpoly, 1)
rotate!(derivedregpoly, 1)
@test coords_x(testregpoly) == coords_x(derivedregpoly)
@test coords_y(testregpoly) == coords_y(derivedregpoly)
scale!(testregpoly, 1)
scale!(derivedregpoly, 1)
@test coords_x(testregpoly) == coords_x(derivedregpoly)
@test coords_y(testregpoly) == coords_y(derivedregpoly)

## Multitype Forward
method1(a::Int, b::Int) = a + b
method2(a::Int, b::Int) = a - b
method3(a::Int, b::Int, c::Int) = a + b + c

@forward {Int, Int} struct Point
    x::Int
    y::Int
end

p = Point(1, 1)

@test method1(p) == 2
@test method3(1, p) == 3
@test method3(p, 1) == 3

module AnotherModule
export testmethod #also work with public
testmethod(a::Int, b::Int) = a + b + 100
privatemethod(a::Int, b::Int) = a + b + 200
end

@forward {Int, Int}, struct Point2
    x::Int
    y::Int
end, (AnotherModule,)

p2 = Point2(1, 1)
@test testmethod(p2) == 102
@test_throws MethodError method1(p2)

@forward {Int, Int} struct Point3
    x::Int
    y::Int
end (
    method1,
    method3,
    AnotherModule,
    AnotherModule.privatemethod
)
p3 = Point3(1, 1)

@test method1(p3) == 2
@test_throws MethodError method2(p3)
@test method3(1, p3) == 3
@test method3(p3, 1) == 3
@test testmethod(p3) == 102
@test privatemethod(p3) == 202


@forward T,
struct StackedCall
    t::T
end (AnotherModule)

@forward T struct CallingSequence
    t::T
end, (AnotherModule,)


