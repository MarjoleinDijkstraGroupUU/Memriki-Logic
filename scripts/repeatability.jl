using DrWatson
@quickactivate

using ACME, Statistics, JLD2

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
    "Nloops" => 100.0,
    "amp" => 5.0,
    "vh" => 1.0,
    "vl" => -1.0,
    "noise" => 0.003)

Vhigh = params["vh"]

is_high(v) = v > 2Vhigh/3
is_low(v) = v < Vhigh/3
is_invalid(v) = Vhigh/3 < v < 2Vhigh/3

function test_gate(gate_circ, params, expected, name)
    circ = @circuit begin
        gate = gate_circ
        Vin1 = voltagesource()
        Vin2 = voltagesource()
        Vout = voltageprobe()
        gate["in1"] == Vin1[+]
        gate["in2"] == Vin2[+]
        gate["out+"] == Vout[+]
        gate["gnd"] == Vin1[-] == Vin2[-] == Vout[-]
    end

    model = DiscreteModel(circ, dt)
    fullname = datadir("repeatability", name*".jld2")
    if isfile(fullname)
        y = jldopen(fullname)["y"]
    else
        y = run!(model, input)'
        mkpath(dirname(fullname))
        jldsave(fullname; y, input, ts, params)
    end

    compare_results(y[preskip:end], expected, params)
end

function compare_results(y, expected, params)
    Ndl =  1 / (4params["f"] * dt)
    vals = [mean(y[round(Int,1+(i-1)*Ndl):round(Int,i*Ndl)]) for i in 1:Int(4params["Nloops"])]
    expected_all = repeat(expected, Int(params["Nloops"]))
    report(vals,expected_all)
end

function report(vals, expected)
    println("high correct in ", sum(is_high.(vals[expected])), " of ", sum(expected))
    println("low correct in ", sum(is_low.(vals[.!expected])), " of ", sum(.!expected))
    println("invalid in ", sum(is_invalid.(vals)), " of ", size(vals,1))
end


f = params["f"]
dt = params["dt"] / f
ts = -1/f:dt:params["Nloops"]/f
preskip = round(Int,1/params["dt"])

input = [params["noise"] * randn() + (1.0 + params["noise"] * randn()) * clamp((1.0 + params["noise"] * randn()) * params["amp"] * sinpi(2params["f"] * t + p / 2), 0.0, params["vh"]) for p in [1, 2], t in ts]


function test_nand()
    test_gate(logic_nand(params), params, [true,true,true,false], "nand")
end

function test_xor()
    test_gate(logic_xor(params), params, [true,false,true,false], "xor")
end

function test_and()
    y, _, _ = and_gate(params, datadir("repeatability", "and.jld2"))
    compare_results(last(y)'[preskip:end], [false, false, false, true], params)
end

function test_or()
    y, _, _ = or_gate(params, datadir("repeatability", "or.jld2"))
    compare_results(last(y)'[preskip:end], [true, false, true, true], params)
end

test_nand()
test_xor()
test_and()
test_or()

