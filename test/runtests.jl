using FitEllipse
using Test

@testset "FitEllipse.jl" begin

function test_fit_ellipse(;
        θ = π/3, a = 3, b = 1.5, x_0 = 3, y_0 = -1,
        N = 10000, ε = 1e-1, plot = false, ξ = 0.001,
    )
    @assert a > b
    x, y = ellipse_from_parametric(a, b, θ, x_0, y_0, N)
    # add noise
    xξ = x .+ randn(N)*ξ;
    yξ = y .+ randn(N)*ξ;

    af, bf, θf, x_0f, y_0f, pf = fit_ellipse(xξ, yξ)
    xo, yo = ellipse_from_parametric(af, bf, θf, x_0f, y_0f)

    @show a, b, θ, x_0, y_0
    @show af, bf, θf, x_0f, y_0f

    if plot # requires loaded GLMakie
        fig, ax = lines(x, y; label = "generator")
        scatter!(ax, xξ, yξ; label = "input")
        lines!(ax, xo, yo; label = "output")
        axislegend(ax)
        display(fig)
    end

    if bf > af
        af, bf, θf = bf, af
    end

    @test isapprox(af, a; atol = ε, rtol = ε)
    @test isapprox(bf, b; atol = ε, rtol = ε)
    @test isapprox(x_0f, x_0; atol = ε, rtol = ε)
    @test isapprox(y_0f, y_0; atol = ε, rtol = ε)
end

test_fit_ellipse()
test_fit_ellipse(; a = 2.0)
test_fit_ellipse(; θ = 0)
test_fit_ellipse(; θ = π/2, ξ = 0.1, ε = 0.1)

end
