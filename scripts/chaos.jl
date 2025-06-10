using DrWatson
@quickactivate

using ACME, JLD2, Random, Statistics, PyPlot

include(srcdir("circs.jl"))
include("../src/circs.jl")
include(srcdir("cluster.jl"))
include("../src/cluster.jl")

params = Dict(
    "circ" => "memriki",
    # "circ" => "nomemriki",
    "nr" => 1e+1,
    "c1" => 1e-4,
    "r1" => 1e-5,
    "xtau" => 10,
    "tau" => τ,
    "co" => 1e-3,
    "lo" => 1e-3,
    "ro" => 1e-3,
    "f" => 700.0,
    "dt" => 1e-4,
    "Nloops" => 200.0,
    "amp" => 0.25,
    "noise" => 0.0)

function simulate(params, save=nothing)
    f = params["f"]
    dt = params["dt"] / f
    ts = 0:dt:params["Nloops"]/f

    input = [(1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2f * t) for _ in 1:1, t in ts]

    if params["circ"] == "memriki"
        base_circ = shinriki_element
    elseif params["circ"] == "nomemriki"
        base_circ = nomemriki_element
    end

    circ = @circuit begin
        memriki = base_circ(params)
        Vin = voltagesource()
        Vout = voltageprobe()
        memriki["in+"] == Vin[+]
        memriki["out+"] == Vout[+]
        memriki["-"] == Vin[-] == Vout[-]
    end
    model = DiscreteModel(circ, dt)
    y = run!(model, input, showprogress=false)
    if save == true || (save === nothing && params["Nloops"] > 999)
        @tagsave(
            datadir("shinriki", savename(params, "jld2")),
            @strdict y input ts params circ
        )
    end
    y, input, ts
end

function phase_plot(fs, rs, recs; vmax=7, ax=nothing, kwargs...)
    if ax === nothing
        fig, ax = subplots()
    end
    map = ax.imshow(length.(recs), origin="lower", extent=[log10.(extrema(rs))..., log10.(extrema(fs))...]; rasterized=true, vmin=0.5, vmax=vmax + 0.5, kwargs...)
    statepoints = [Dict("f" => 120, "r1" => 0.425, "sym"=>"<", "color"=>"deeppink"), Dict("f" => 400, "r1" => 0.425, "sym"=>"^", "color"=>"tab:green"), Dict("f" => 400, "r1" => 1.0, "sym"=>">", "color"=>"tab:orange")]
    for st in statepoints
        ax.plot(log10(st["r1"]), log10(st["f"]), marker=st["sym"], color=st["color"])
    end
    ax,map
end


fs = 10 .^ range(1, 3, length=100)
rs = 10 .^ range(-1, 2, length=100)

# simulate

for c in ["memriki", "nomemriki"]
    params["circ"] = c
    recs = fr2recs(params, fs, rs)
    @tagsave(
        datadir("chaos-sweeps", savename(params, "jld2")),
        @strdict fs rs params recs
    )
end

# plot

skip = 30
vmax = 20

cmap = get_cmap("viridis",vmax)
fig, axs = subplots(2, 1, figsize=(5, 7), layout="compressed", sharex=true, sharey=true)
gmap = nothing
for (ax, circ) in zip(axs, ["memriki", "nomemriki"])
    params["circ"] = circ
    file = jldopen(datadir("chaos-sweeps", savename(params, "jld2")))
    _, map = phase_plot(fs, rs, file["recs"]; vmax, ax, cmap)
    global  gmap = map
end
norm = plt.matplotlib.colors.BoundaryNorm(1:vmax+1,cmap.N)
sm = plt.matplotlib.cm.ScalarMappable(norm=norm)
for ax in axs
    ax.set_ylabel(L"\mathrm{log}_{10}(f /\mathrm{Hz})")
end
axs[end].set_xlabel(L"\mathrm{log}_{10}(R_1 g_0)")
fig.colorbar(gmap; ticks=1:vmax+1, ax=axs, orientation="horizontal", location="top", aspect=40, extend="max", fraction=0.05,pad=0.005)
fig.savefig(plotsdir("phase.pdf"), dpi=600)
