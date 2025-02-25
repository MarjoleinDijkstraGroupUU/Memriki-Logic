using DrWatson
@quickactivate

using ACME, PyPlot

include(srcdir("circs.jl"))
include(srcdir("combined_gates.jl"))

params = Dict(
    "c1" => 1e-4,
    "r1" => 0.4,
    "tau" => τ,
    "co" => 1e-3,
    "lo" => 1e-3,
    "ro" => 1e-3,
    "f" => 0.3,
    "dt" => 1e-5,
    "Nloops" => 1,
    "amp" => 5.0,
    "vh" => 1.0,
    "vl" => -1.0,
    "noise" => 0.003)


vertical = false

if vertical
    fig, ax = subplots(2, layout="compressed", sharex=true, figsize=(5, 4))
else
    fig, ax = subplots(1, 2, layout="compressed", sharex=true, sharey=true, figsize=(5, 2))
end

ys, ins, ts = and_gate(params, datadir("combined_gates", "and.jld2"))
plot_sequential(last(ys), ins[1][1:2, :], ts; ax=ax[1], skip=100)

ys, ins, ts = or_gate(params, datadir("combined_gates", "or.jld2"))
plot_sequential(last(ys), ins[1][1:2:3, :], ts; ax=ax[2], skip=100)

for a in ax
    a.set_xlim(0, params["Nloops"] / params["f"])
    if vertical
        a.set_ylabel("voltage [V]")
    else
        a.set_xlabel("time [s]", labelpad=-5)
    end
end

if vertical
    ax[2].set_xlabel("time [s]", labelpad=-5)
    else
    ax[1].set_ylabel("voltage [V]")
end

fig.savefig(plotsdir("new-andor-$(vertical ? "v" : "h").pdf"), bbox_inches="tight")
