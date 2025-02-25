using DrWatson
@quickactivate

using ACME, PyPlot, JLD2

include(srcdir("circs.jl"))


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

f = params["f"]
dt = params["dt"] / f
ts = -1/f:dt:params["Nloops"]/f # add one additional loop for equilibration

input = [params["noise"] * randn() + (1.0 + params["noise"] * randn()) * clamp((1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2params["f"] * t + p / 2), 0.0, params["vh"]) for p in [1, 2], t in ts]

function sim_gate(gate_circ)
    circ = @circuit begin
        gate = gate_circ
        Vin1 = voltagesource()
        Vin2 = voltagesource()
        Vout = voltageprobe()
        gate["in1"] == Vin1[+]
        gate["in2"] == Vin2[+]
        gate["out+"] == Vout[+]
        gate["out-"] == Vout[-]
        gate["in-"] == Vin1[-] == Vin2[-]
    end
    model = DiscreteModel(circ, dt)
    y = run!(model, input)
end

vertical = true

if vertical
    fig, axs = subplots(2, layout="compressed", figsize=(5, 4), sharex=true)
else
    fig, axs = subplots(1, 2, layout="compressed", figsize=(5, 2), sharex=true, sharey=true)
end
for (ax, gate, name) in zip(axs, [xor_element(params), nand_element(params)], ["xor", "nand"])
    fullname = datadir("single-gates", name * ".jld2")
    if isfile(fullname)
        y = jldopen(fullname)["y"]
    else
        y = sim_gate(gate)
        mkpath(dirname(fullname))
        jldsave(fullname; y, input, params)
    end
    for (j, ini) in pairs(eachrow(input))
        ax.plot(ts, ini, label="input $j")
    end
    ax.plot(ts, y', label="output")
    ax.legend()
    ax.set_xlabel("time [s]")
    ax.set_xlim(0, params["Nloops"] / params["f"])
    if vertical
        ax.set_ylabel("voltage [V]")
    end

end
if !vertical
    axs[1].set_ylabel("voltage [V]")
    axs[2].yaxis.set_label_position("right")
    axs[2].yaxis.tick_right()
end

fig.savefig(plotsdir("nand-xor-$(vertical ? "v" : "h").pdf"))
