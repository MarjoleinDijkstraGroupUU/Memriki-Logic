using DrWatson
@quickactivate
using ACME, PyPlot,JLD2

include(srcdir("circs.jl"))

function simulate(params, input, dt)
    circ = @circuit begin
        memriki = shinriki_element(params)
        Vin = voltagesource()
        Vout = voltageprobe()
        memriki["in+"] == Vin[+]
        memriki["out+"] == Vout[+]
        memriki["-"] == Vin[-] == Vout[-]
    end
    model = DiscreteModel(circ, dt)
    run!(model, input, showprogress=true)
end

ogparams = Dict(
    "circ" => "qr",
    "nr" => 1e+1,
    "c1" => 1e-4,
    "tau" => τ,
    "co" => 1e-3,
    "lo" => 1e-3,
    "ro" => 1e-3,
    "dt" => 1e-5,
    "Nloops" => 1100.0,
    "amp" => 0.25,
    "noise" => 0.0)

overwrides = [Dict("f" => 120, "r1" => 0.425), Dict("f" => 400, "r1" => 0.425), Dict("f" => 400, "r1" => 1.0)]


fig, axs = subplots(2, 3, figsize=(10,3))
ogpreskip = 100
skip = 1000
for (i, ow) in enumerate(overwrides)
    params = merge(ogparams, ow)
    f = params["f"]
    dt = params["dt"] / f
    ts = 0:dt:params["Nloops"]/f

    input = [(1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2f * t) for _ in 1:1, t in ts]
    name = datadir("dynamics", savename(ow, "jld2"))
    if isfile(name)
        file = jldopen(name)
        ys = file["ys"]
        input = file["input"]
        ts = file["ts"]
    else
        ys = simulate(params, input, dt)
        mkpath(dirname(name))
        jldsave(name; ys, params, input, ts)
    end

    preskip = round(Int, ogpreskip / params["dt"]) + 1
    axs[1, i].set_title(latexstring("R_1=", params["r1"], "\\,\\mathrm{g_0}^{-1}, f=", params["f"], "\\,\\mathrm{Hz}"))
    tss = ts[preskip:skip:end]
    ins = input[1, preskip:skip:end]
    y = ys[preskip:skip:end]
    axs[1, i].plot(tss .* params["f"], y, label=L"U_\mathrm{out}", zorder=-2)
    axs[1, i].plot(tss .* params["f"], ins, label=L"U_\mathrm{in}", zorder=-3)
    axs[1, i].set_xlim(100, 125)
    axs[2, i].plot(ins, y, linewidth=20 / params["Nloops"], zorder=-2)
    axs[2, i].set_xlabel(L"U_\mathrm{in}~[\mathrm{V}]")
    axs[1, i].set_xlabel(L"tf", labelpad=-5.5)
end
axs[1, 1].set_ylabel(L"U~[\mathrm{V}]")
axs[1, 1].legend(loc="upper right", frameon=true, framealpha=0.8)
axs[2, 1].set_ylabel(L"U_\mathrm{out}~[\mathrm{V}]")
axs[2, 2].set_rasterization_zorder(0)

for (ax, c) in zip(axs, "adbecf")
    ax.text(-0.1, 0.95, string("\\textbf ", c), transform=ax.transAxes)
end

savefig(plotsdir("dynamics.pdf"), dpi=1200)
