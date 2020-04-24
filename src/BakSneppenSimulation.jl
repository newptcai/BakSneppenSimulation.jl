module BakSneppenSimulation

export drawball, swtichball, test, runsim, plotdata, Simulation, simdirected, simcircular, simdirectedbounded

ENV["GKS_ENCODING"]="utf8"
using Random, Plots

struct DataPoint
    time::Int
    zeros::Int
    width::Int
    diameter::Int
    death::Int
end

struct Simulation
    xlist::Array{Int, 1}
    x0list::Array{Int, 1}
    p::Float64
    mutnum::Int
    stop::BigInt
    pointnum::Int
    circular::Bool
    n::Int
    step::Int
    bounded::Bool
end

n(sim::Simulation) = length(sim.xlist)

function Simulation(xlist::Array{Int,1}, p::Float64, mutnum::Int, stop::Int, pointnum::Int=10^3,
                    bounded::Bool=false, circular::Bool=false)
    x0list = findall(x->x==0, xlist)
    step = max(1, floor(stop/pointnum))
    Simulation(xlist, x0list, p, mutnum, stop, 10^3, circular, length(xlist), step, bounded)
end

function circularind(ind, len)
    if ind < 1
        len + ind
    elseif len < ind
        ind - len
    else
        ind
    end
end

function diameter(x0list, n)
    zeros = length(x0list)
    if zeros == 0
        return 0
    end

    if zeros == 1
        return 1
    end

    sort!(x0list)
    wlist = similar(x0list)
    for i in 1:length(wlist)
        prev = circularind(i-1, zeros)
        next = circularind(i+1, zeros)

        w1 = mod(x0list[next] - x0list[i] - 1, n)
        w2 = mod(x0list[i] - x0list[prev] - 1, n)
        wlist[i] = n - max(w1, w2)
    end
    findmin(wlist)[1]
end

function simdirected(n::Int, p::Float64, mutnum::Int, stop::Int, returnsim::Bool=false)
    xlist = ones(Int, n)
    xlist[1] = 0
    sim = Simulation(xlist, p, mutnum, stop, 10^3, false, false)
    data = runsim(sim)
    if returnsim
        return (sim, data)
    else
        return data
    end
end

function simdirectedbounded(n::Int, p::Float64, mutnum::Int, stop::Int, returnsim::Bool=false)
    xlist = ones(Int, n)
    xlist[1] = 0
    sim = Simulation(xlist, p, mutnum, stop, 10^3, true, false)
    data = runsim(sim)
    if returnsim
        return (sim, data)
    else
        return data
    end
end

function simcircular(n::Int, p::Float64, mutnum::Int, stop::Int, returnsim::Bool=false)
    xlist = zeros(Int, n)
    sim = Simulation(xlist, p, mutnum, stop, 10^3, false, true)
    data = runsim(sim)
    if returnsim
        return (sim, data)
    else
        return data
    end
end

function width(xlist)
    l = findfirst(x->x==0, xlist)
    r = findlast(x->x==0, xlist)
    r-l+1
end

function runsim(sim::Simulation)
    data = Array{DataPoint, 1}()
    death = 0
    for t in 1:sim.stop
        drawball(sim)
        if length(sim.x0list) == 0
            death+=1
        end
        if (t-1) % sim.step == 0
            z = length(sim.x0list)
            w = 0
            d = 0
            if z != 0
                if !sim.circular
                    w = width(sim.xlist)
                    d = w
                else
                    d = diameter(sim.x0list, n(sim))
                end
            end
            push!(data, DataPoint(t, z, w, d, death))
        end
    end
    data
end

#runsim(n, p) = runsim(n, p, 2*n)

function drawball(sim::Simulation)
    x0len = length(sim.x0list)
    x0ind = 0
    ind = 0
    if x0len> 0
        x0ind = rand(1:x0len)
        ind = sim.x0list[x0ind]
    else
        ind = 1
    end
    newball = map(x-> (x < sim.p ? 1 : 0), rand(sim.mutnum))

    swtichball(sim, newball, ind)
end


function swtichball(sim::Simulation, newball::Array{Int}, ind::Int)
    shift = convert(Int, ceil(sim.mutnum/2))
    for i in 1:sim.mutnum
        ballind = ind + i - shift
        if !sim.circular
            if ballind < 1
                continue
            elseif  ballind == n(sim)+1
                # extend the array
                if !sim.bounded
                push!(sim.xlist, 1)
                else
                    continue
                end
            end
        else
            if ballind < 1 
                ballind = n(sim) + ballind
            elseif ballind >n(sim)
                ballind = ballind -n(sim)
            end
        end

        new = newball[i]
        old = sim.xlist[ballind]

        if new == old
            continue
        end

        sim.xlist[ballind] = new
        if new == 1 && old == 0
            filter!(a -> a != ballind, sim.x0list)
        else
            push!(sim.x0list, ballind)
        end
    end
end

function plotdata(data::Array{DataPoint, 1}, n::Number=Inf, showwidth::Bool=false)
    gr()
    tl = map(p->p.time, data)
    zl = map(p->p.zeros, data)
    wl = map(p->p.width, data)
    dl = map(p->p.diameter, data)
    p = plot(tl, zl, ylims=(-1,n), label="zeros")
    if showwidth
        plot!(tl, wl, label="width")
    end
    plot!(tl, dl, label="diameter")
end

end # module
