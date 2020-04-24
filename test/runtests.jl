using BakSneppenSimulation
using BakSneppenSimulation: n, width, diameter
using Test
using Random

@testset "BakSneppenSimulation.jl" begin
    p0 = 0.6

    Random.seed!(1234)
    xlist = rand(0:1, 10)
    sim = Simulation(xlist, p0, 3, 10^6)

    @testset "3 mutations" begin
        drawball(sim)

        @test sim.xlist == [0, 1, 1, 1, 0, 0, 0, 0, 1, 1]
        @test sim.x0list == findall(x->x==0, xlist)

        @test diameter(sim.x0list, n(sim)) == 7
        @test width(sim.xlist) == 8

        drawball(sim)

        @test sim.xlist == [0, 1, 1, 1, 0, 0, 1, 1, 1, 1]
        @test sim.x0list == findall(x->x==0, xlist)

        @test diameter(sim.x0list, n(sim)) == 6
        @test width(sim.xlist) == 6

        sim = Simulation(xlist, p0, 2, 10^6)
    end

    @testset "2 mutations" begin
        xlist = rand(0:1, 10)
        sim = Simulation(xlist, p0, 2, 10^6)

        drawball(sim)
        @test sim.xlist == [0, 1, 1, 1, 1, 1, 0, 1, 0, 1]

        drawball(sim)
        @test sim.xlist == [0, 1, 1, 1, 1, 1, 0, 1, 0, 1]

        drawball(sim)
        @test sim.xlist == [1, 0, 1, 1, 1, 1, 0, 1, 0, 1]

        drawball(sim)
        @test sim.xlist == [1, 0, 1, 1, 1, 1, 0, 0, 0, 1]
    end
end
