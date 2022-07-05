module FitEllipse

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end FitEllipse

export fit_ellipse, ellipse_from_parametric

"""
    fit_ellipse(x, y) → a, b, θ, x₀, y₀, p
Fit an ellipse into the 2D data defined by `(x, y)` using ordinary least squares.

Return: semi-major axis length, semi-minor axis length, ellipse rotation,
center coordinates and a parameter container for quadratic form of ellipse.
The first five parameters represent the parametric form of an ellipse,
which may have an arbitrary rotation and translation w.r.t. the origin:
```julia
x(t) = cos(θ)*a*cos(t) - sin(θ)*b*sin(t) + x₀
y(t) = sin(θ)*a*cos(t) + cos(θ)*b*sin(t) + y₀
t ∈ [0, 2π)
```
Use `ellipse_from_parametric(a, b, θ, x₀, y₀) → x, y` to construct points
along the ellipse. Notice that always two possible ellipses can fit the data,
one with `a, b, θ` and one with `b, a, θ-π/2`.

The quadratic form of the ellipse is created using `p`,
and is given by the form `(x^2, x*y, y^2, x, y) ⋅ p = 1`.

Code modified from:
https://www.matecdev.com/posts/julia-least-squares-qr.html
using a lot of stuff from:
https://en.wikipedia.org/wiki/Ellipse#General_ellipse
"""
function fit_ellipse(x, y)
    @assert length(x) == length(y)
    # design matrix with columns the quadratic form of ellipse
    M = hcat(x.^2, x.*y, y.^2, x, y)
    p = M\ones(length(x)) # best fit parameters for the ellipse
    A, B, C, D, E = p
    F = -1.0
    # Calculate the "useful" parametric-form parameters from quadratic (Wikipedia):
    Δ = B^2 - 4*A*C
    Λ = (A-C)^2 + B^2
    # Notice that we clamp the square root because I've noticed extreme cases
    # where we get something like `sqrt(-7.559810648707566e-26)` (its obviously zero)
    a, b = [-sqrt(
            clamp(
            2 * (A*E^2 + C*D^2 - B*D*E + (Δ)*F) *
            ( (A+C) + op(sqrt(Λ)) ), 0, Inf)) / (B^2 - 4A*C)   for op in (+, -)
        ]
    θ = atan((C - A - sqrt(Λ))/B)
    x₀ = (2*C*D - B*E)/Δ
    y₀ = (2*A*E - B*D)/Δ
    return a, b, θ, x₀, y₀, p
end

function ellipse_from_parametric(a, b, θ, x₀, y₀, N = 1000)
    fx = (t) -> cos(θ)*a*cos(t) - sin(θ)*b*sin(t) + x₀
    fy = (t) -> sin(θ)*a*cos(t) + cos(θ)*b*sin(t) + y₀
    ts = range(0, prevfloat(2π); length = N)
    return fx.(ts), fy.(ts)
end


end # module
