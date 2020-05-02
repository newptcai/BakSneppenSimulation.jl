using BakSneppenSimulation
using BakSneppenSimulation: n, diameter
using Test
using Random

@testset "BakSneppenSimulation.jl" begin

    @testset "3 mutations" begin
        Random.seed!(1234)

        config = SimConfig(0.6, 3)
        state = SimState(rand(0:1, 10))
        sim = SimCircular(config, state)

        @test state.xlist == [0, 1, 1, 1, 0, 0, 0, 1, 0, 1]
        @test state.x0list == [1, 5, 6, 7, 9]

        mutate(sim)

        @test state.xlist == [0, 1, 1, 1, 1, 0, 0, 1, 0, 1]
        @test state.x0list == findall(x->x==0, state.xlist)
        @test diameter(sim) == 6

        mutate(sim)

        @test state.xlist == [0, 1, 1, 1, 1, 1, 0, 1, 0, 1]
        @test state.x0list == findall(x->x==0, state.xlist)
        @test diameter(sim) == 5

        mutate(sim)

        @test state.xlist == [0, 1, 1, 1, 1, 1, 0, 0, 0, 0]
        @test sort!(state.x0list) == findall(x->x==0, state.xlist)
        @test diameter(sim) == 5
    end

    @testset "2 mutations" begin
        Random.seed!(1234)

        config = SimConfig(0.6, 2)
        state = SimState(rand(0:1, 10))
        sim = SimLinear(config, state, true)

        @test state.xlist == [0, 1, 1, 1, 0, 0, 0, 1, 0, 1]
        @test state.x0list == [1, 5, 6, 7, 9]

        mutate(sim)

        @test state.xlist == [0, 1, 1, 1, 1, 1, 0, 1, 0, 1]
        @test state.x0list == findall(x->x==0, state.xlist)
        @test diameter(sim) == 9

        mutate(sim)

        @test state.xlist == [0, 1, 1, 1, 1, 1, 0, 1, 0, 1]
        @test state.x0list == findall(x->x==0, state.xlist)
        @test diameter(sim) == 9

        mutate(sim)

        @test state.xlist == [1, 0, 1, 1, 1, 1, 0, 1, 0, 1]
        @test sort!(state.x0list) == findall(x->x==0, state.xlist)
        @test diameter(sim) == 8
    end
end
