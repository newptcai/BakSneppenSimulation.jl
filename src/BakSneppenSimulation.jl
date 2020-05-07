module BakSneppenSimulation

export applymutations!,
    mutate,
    test,
    runsim,
    plotsim,
    plotzero,
    simlinear,
    Simulation,
    SimCircular,
    SimLinear,
    SimConfig,
    SimState

import Base: show

ENV["GKS_ENCODING"] = "utf8"
using Random, Plots

struct DataPoint
    time::Int
    zeros::Int
    diameter::Int
    death::Int
end

abstract type Simulation end

show(io::IO, sim::Simulation) = print(
    io,
    "Simulation for p=$(sim.config.p) and $(sim.config.mutnum) mutations",
)

struct SimConfig
    n::Int
    p::Float64
    mutnum::Int
    stop::Int
    pointnum::Int
    step::Int
    shift::Int
end

function SimConfig(
    n::Int,
    p::Float64,
    mutnum::Int,
    stop::Int = 10^6,
    pointnum::Int = 10^3,
)
    step = max(1, floor(stop / pointnum))
    shift = convert(Int, ceil(mutnum / 2))
    SimConfig(n, p, mutnum, stop, pointnum, step, shift)
end

struct SimState
    xlist::Array{Int,1}
    x0list::Array{Int,1}
end

function SimState(xlist::Array{Int,1})
    x0list = findall(x -> x == 0, xlist)
    SimState(xlist, x0list)
end

struct SimCircular <: Simulation
    config::SimConfig
    state::SimState
end

function SimCircular(config::SimConfig)
    xlist = zeros(Int, config.n)
    state = SimState(xlist)
    SimCircular(config, state)
end

struct SimLinear <: Simulation
    config::SimConfig
    state::SimState
    bounded::Bool
end

function SimLinear(config::SimConfig, bounded = false)
    xlist = ones(Int, config.n)
    xlist[1] = 0
    state = SimState(xlist)
    SimLinear(config, state, bounded)
end

n(sim::Simulation) = length(sim.state.xlist)
zeros(sim::Simulation) = length(sim.state.x0list)

function circularind(ind, len)
    if ind < 1
        len + ind
    elseif len < ind
        ind - len
    else
        ind
    end
end

function circulardist(dist, len)
    if dist < 0
        dist + len
    else
        dist
    end
end

function diameter(sim::SimCircular)
    num0 = zeros(sim)
    if num0 == 0
        return 0
    end

    if num0 == 1
        return 1
    end

    x0list = sim.state.x0list
    sort!(x0list)
    wlist = similar(x0list)
    for i = 1:length(wlist)
        prev = circularind(i - 1, num0)
        next = circularind(i + 1, num0)

        w1 = circulardist(x0list[next] - x0list[i] - 1, n(sim))
        w2 = circulardist(x0list[i] - x0list[prev] - 1, n(sim))
        wlist[i] = n(sim) - max(w1, w2)
    end
    findmin(wlist)[1]
end

function diameter(sim::SimLinear)
    l = findfirst(x -> x == 0, sim.state.xlist)
    r = findlast(x -> x == 0, sim.state.xlist)
    r - l + 1
end

function runsim(sim::Simulation)
    data = Array{DataPoint,1}()
    death = 0
    for t = 1:sim.config.stop
        mutate(sim)
        if zeros(sim) == 0
            death += 1
        end
        if (t - 1) % sim.config.step == 0
            z = zeros(sim)
            d = 0
            if z != 0
                d = diameter(sim)
            end
            push!(data, DataPoint(t, z, d, death))
        end
    end
    data
end

function simlinear(config::SimConfig, bounded = false)
    sim = SimLinear(config, bounded)
    data = runsim(sim)
    (data, sim)
end

function mutate(sim::Simulation)
    x0len = zeros(sim)
    x0ind = 0
    ind = 0
    if x0len > 0
        x0ind = rand(1:x0len)
        ind = sim.state.x0list[x0ind]
    else
        ind = 1
    end

    mutations = map(x -> (x < sim.config.p ? 1 : 0), rand(sim.config.mutnum))

    applymut!(sim, mutations, ind)
end


function applymut!(sim::SimLinear, mutations::Vector{Int}, ind::Int)
    @debug "applymut!" mutations, ind
    for i = 1:sim.config.mutnum
        mutationind = ind + i - sim.config.shift
        if mutationind < 1
            continue
        elseif mutationind == n(sim) + 1
            # extend the array
            if !sim.bounded
                push!(sim.state.xlist, 1)
            else
                continue
            end
        end

        update!(sim.state, mutations[i], mutationind)
    end
end

function applymut!(sim::SimCircular, mutations::Array{Int}, ind::Int)
    @debug "applymut!" mutations, ind
    for i = 1:sim.config.mutnum
        mutationind = circularind(ind + i - sim.config.shift, n(sim))

        update!(sim.state, mutations[i], mutationind)
    end
end

function update!(state::SimState, new, mutationind)
    old = state.xlist[mutationind]

    if new == old
        return
    end

    state.xlist[mutationind] = new
    if new == 1 && old == 0
        filter!(a -> a != mutationind, state.x0list)
    else
        push!(state.x0list, mutationind)
    end
end

function plotsim(data::Array{DataPoint,1}, sim)
    gr()
    tl = map(p -> p.time, data)
    zl = map(p -> p.zeros, data)
    dl = map(p -> p.diameter, data)
    ylim = max(findmax(dl)[1], sim.config.n)
    p = plot(tl, zl, ylims = (-1, ylim), label = "zeros")
    plot!(tl, dl, label = "diameter", title = string(sim))
end
plotsim((data, sim)) = plotsim(data, sim)

function plotzero(data, sim)
    gr()
    tl = map(p -> p.time, data)
    pl = map(p -> p.zeros / p.diameter, data)
    p = plot(
        tl,
        pl,
        ylims = (0, 1),
        label = "percentage of zeros",
        title = string(sim),
    )
end
plotzero((data, sim)) = plotzero(data, sim)

end # module
