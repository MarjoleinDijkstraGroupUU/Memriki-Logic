using DrWatson
@quickactivate

using ACME, PyPlot, JLD2

include(srcdir("circs.jl"))

function nand1(params)
    @circuit begin
        nand = logic_nand(params)
        in1 = voltagesource()
        in2 = voltagesource()
        out = voltageprobe()
        in1[+] == nand["in1"]
        in2[+] == nand["in2"]
        out[+] == nand["out+"]
        in1[-] == in2[-] == out[-] == nand["gnd"]
    end
end

function sequential_nands(params, input, name=nothing)
    if isfile(name)
        file = jldopen(name)
        ys = file["ys"]
        ins = file["ins"]
        ts = file["ts"]
        ys, ins, ts
    else
        circ = nand1(params)

        f = params["f"]
        dt = params["dt"] / f
        ts = -1/f:dt:params["Nloops"]/f # add one additional loop for equilibration

        model = DiscreteModel(circ, dt)
        y1 = run!(model, input[1:2, :])
        y2 = run!(model, input[3:4, :])

        in3 = hcat(y1[end, :] + params["noise"] .* randn(size(y1[end, :])), y2[end, :] + params["noise"] .* randn(size(y2[end, :])))'
        y = run!(model, in3)
        ys = (y1, y2, y)
        ins = (input, input, in3)
        mkpath(dirname(name))
        jldsave(name; ys, ins, ts, params)
        ys, ins, ts
    end
end

function plot_sequential(y, input, ts; ax=nothing, fig=nothing, skip=1)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    ax.plot(ts[1:skip:end], input[1, 1:skip:end], label="input 1")
    ax.plot(ts[1:skip:end], input[2, 1:skip:end], label="input 2")
    ax.plot(ts[1:skip:end], y[end, 1:skip:end], label="output")
    ax.legend()
    fig, ax
end

function or_gate(params, name)
    ts = (-1:params["dt"]:params["Nloops"]) ./ params["f"]
    input = [params["noise"] * randn() + (1.0 + params["noise"] * randn()) * clamp((1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2params["f"] * t + p / 2), 0.0, params["vh"]) for p in [1, 1, 2, 2], t in ts]
    sequential_nands(params, input, name)
end

function and_gate(params, name)
    ts = (-1:params["dt"]:params["Nloops"]) ./ params["f"]
    input = [params["noise"] * randn() + (1.0 + params["noise"] * randn()) * clamp((1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2params["f"] * t + p / 2), 0.0, params["vh"]) for p in [1, 2, 1, 2], t in ts]
    sequential_nands(params, input, name)
end
