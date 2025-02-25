using DrWatson
@quickactivate

using ACME, PyPlot, JLD2

include(srcdir("circs.jl"))

circ = @circuit begin
    m = memristor()
    V = voltagesource()
    I = currentprobe()
    V[+] == m[1]
    m[2] == I[+]
    I[-] == V[-]
end

const use_euler = false

fig, ax = subplots(1, 2, layout="compressed", sharex=true, figsize=(5,2))
for f in [0.4, 40, 4000]
    dt = 1e-5 / f
    ts = 0.0:dt:200/f
    skip = round(Int, 199 / (f * dt))

    V(t) = sin(2π * f * t)
    input = V.(ts)'

    if use_euler
        gss = zeros(size(ts))
        gss[1] = g∞(V(ts[1]))
        for i in 2:length(gss)
            gss[i] = gss[i-1] + dt / τ * (g∞(V(ts[i])) - gss[i-1])
        end

        ax[1].plot(V.(ts[skip:end]), gss[skip:end], label="f=$f\\,Hz")
        ax[2].plot(V.(ts[skip:end]), gss[skip:end] .* V.(ts[skip:end]), label="f=$f\\,Hz")
    else
        name = datadir("memristor", "data$f.jld2")
        if isfile(name)
            y = jldopen(name)["y"]
        else
            model = DiscreteModel(circ, dt)
            y = run!(model, input)
            mkpath(dirname(name))
            jldsave(name; y)
        end
        ax[1].plot(input'[skip:end], y'[skip:end] ./ input'[skip:end], label="f=$f\\,Hz")
        ax[2].plot(input'[skip:end], y'[skip:end], label="f=$f\\,Hz")
    end
end

ax[1].set_xlabel("voltage [V]")
ax[1].set_ylabel("conductance [\$\\mathrm{g_0}\$]")
ax[2].set_xlabel("voltage [V]")
ax[2].set_ylabel("current [\$\\mathrm{g_0}\$V]")
ax[2].yaxis.set_label_position("right")
ax[2].yaxis.tick_right()
ax[2].legend()

fig.savefig(plotsdir("fig1sub.pdf"))
