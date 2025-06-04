using DrWatson
@quickactivate

using ACME, StaticArrays

include(srcdir("memristor.jl"))

memristor(; τ=τ, g0=1.0) = ACME.Element(
    mv=[-1; 0; 0; 0], mi=[0; -1; 0; 0], mq=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], mx=[0; 0; -1; 0], mxd=[0; 0; 0; -1];
    nonlinear_eq=function (q)
        v, i, g, gd = q
        g∞, g∞′ = g0 .* gs(v)
        scale = 1 / g0
        res = @SVector [scale * (v * g - i), scale * (g∞ - g - τ * gd)]
        J = @SMatrix [scale*g -scale scale*v 0; scale*g∞′ 0 -scale -scale*τ]
        return (res, J)
    end
)
memristor(τ=τ, g0=1.0) = memristor(; τ, g0)

qubic_conductance(g1, g3) = ACME.Element(
    mv=[1; 0], mi=[0; 1], mq=[-1 0; 0 -1],
    nonlinear_eq = function (q)
        v,i = q
        res = @SVector [-g1*v-g3*v^3+i]
        J = @SMatrix [-g1-3g3*v^2 1]
        return res, J
    end
)

impedance_changer() = ACME.Element(
    mv=[1 -1; 0 0], mi=[0 0; 1 0],
    ports=["in+" => "in-", "out+" => "out-"]
)

function get_c(c=1e-3)
    circ = @circuit begin
        C1 = capacitor(c)
        C2 = capacitor(c)
    end
    composite_element(circ, pinmap=Dict("in+" => (:C1, 1), "in-" => (:C2, 1), "out+" => (:C1, 2), "out-" => (:C2, 2)))
end

function get_rect(is=1e-3)
    rectifyer_circ = @circuit begin
	      D1 = diode(; is)
        D2 = diode(; is)
        D3 = diode(; is)
        D4 = diode(; is)

        D1[+] == D2[+]
        D3[-] == D4[-]
        D1[-] == D3[+]
        D2[-] == D4[+]
    end
    composite_element(rectifyer_circ, pinmap=Dict("in+" => (:D1, :-), "in-" => (:D4, :+), "out+" => (:D4, :-), "out-" => (:D1, :+)))
end

function shinriki_circ(params)
    @circuit begin
        nr = qubic_conductance(-3, 1)
        c1 = capacitor(params["c1"])
        r1 = resistor(params["r1"])

        m = memristor(τ=params["tau"])
        mr = memristor(τ=params["tau"])

        co = capacitor(params["co"])
        lo = inductor(params["lo"])
        ro = resistor(params["ro"])

        nr[2] == c1[1] == m[1] == mr[2] == r1[2]
        m[2] == mr[1] == co[1] == lo[1]
        lo[2] == ro[1]
        co[2] == ro[2] == c1[2] == r1[1]
    end
end
function shinriki_element(params)
    composite_element(shinriki_circ(params), pinmap=Dict("in+" => (:nr, 1), "-" => (:co, 2), "out+" => (:nr, 2)))
end

function nomemriki_circ(params)
    @circuit begin
        nr = qubic_conductance(-3, 1)
        c1 = capacitor(params["c1"])
        r1 = resistor(params["r1"])

        m = resistor(0.25)
        mr = resistor(0.25)

        co = capacitor(params["co"])
        lo = inductor(params["lo"])
        ro = resistor(params["ro"])

        nr[2] == c1[1] == m[1] == mr[2] == r1[2]
        m[2] == mr[1] == co[1] == lo[1]
        lo[2] == ro[1]
        co[2] == ro[2] == c1[2] == r1[1]
    end
end
function nomemriki_element(params)
    composite_element(nomemriki_circ(params), pinmap=Dict("in+" => (:nr, 1), "-" => (:co, 2), "out+" => (:nr, 2)))
end

function xor_circ(params)
    @circuit begin
        sh1 = shinriki_element(params)
        sh2 = shinriki_element(params)
        sh3 = shinriki_element(params)

        first = get_c()
        sh1["out+"] == first["in+"]
        sh2["out+"] == first["in-"]
        first["out+"] == sh3["in+"]
        first["out-"] == sh3[-]

        sh1[-] == sh2[-]
    end
end

function xor_element(params)
    composite_element(xor_circ(params), pinmap=Dict("in1" => (:sh1, "in+"), "in2" => (:sh2, "in+"), "out+" => (:sh3, "out+"), "out-" => (:sh3, "-"), "in-" => (:sh1, "-")))
end

function nand_circ(params)
    @circuit begin
        sh1 = shinriki_element(params)
        sh2 = shinriki_element(params)
        sh3 = shinriki_element(params)

        first = impedance_changer()
        second = get_c()
        sh1["out+"] == first["in+"]
        sh2["out+"] == first["in-"]
        first["out+"] == second["in+"]
        first["out-"] == second["in-"]
        second["out+"] == sh3["in+"]
        second["out-"] == sh3[-]

        sh1[-] == sh2[-]
    end
end

function nand_element(params)
    composite_element(nand_circ(params), pinmap=Dict("in1" => (:sh1, "in+"), "in2" => (:sh2, "in+"), "out+" => (:sh3, "out+"), "in-" => (:sh1, "-"), "out-"=>(:sh3,"-")))
end

function coupler(params)
    vh = params["vh"] + params["noise"] * params["amp"] / 2
    vl = params["vl"] - params["noise"] * params["amp"] / 2
    circ = @circuit begin
        imp = opamp(Val{:macak}, 15, vl, vh)
        imp2 = opamp(Val{:macak}, 35, vl, vh)
        rect = get_rect(3e-1)
        c = capacitor(5.0e-2)
        r = resistor(1e-2)
        imp["out+"] == rect["in+"]
        imp["out-"] == rect["in-"]
        rect["out+"] == c[1] == r[1]
        rect["out-"] == imp2["in-"] == c[2]
        imp2["in+"] == r[2]
    end
    composite_element(circ, pinmap=Dict("in+"=>(:imp,"in+"), "in-"=>(:imp,"in-"), "out+"=>(:imp2,"out+"), "out-"=>(:imp2,"out-")))
end

function logic_nand(params)
    vh = params["vh"]
    vl = params["vl"]
    circ = @circuit begin
        nand = nand_element(params)
        coup = coupler(params)
        nand["out+"] == coup["in+"]
        nand["out-"] == coup["in-"]
        nand["in-"] == coup["out-"]
    end
    composite_element(circ, pinmap=Dict("in1" => (:nand, "in1"), "in2" => (:nand, "in2"), "out+" => (:coup, "out+"), "gnd" => (:coup, "out-")))
end

function logic_xor(params)
    vh = params["vh"]
    vl = params["vl"]
    circ = @circuit begin
        nand = xor_element(params)
        coup = coupler(params)
        nand["out+"] == coup["in+"]
        nand["out-"] == coup["in-"]
        nand["in-"] == coup["out-"]
    end
    composite_element(circ, pinmap=Dict("in1" => (:nand, "in1"), "in2" => (:nand, "in2"), "out+" => (:coup, "out+"), "gnd" => (:coup, "out-")))
end
