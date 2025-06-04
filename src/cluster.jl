using DrWatson
@quickactivate

using Folds

function naive_cluster(input, tol=1e-2)
    clusters = Set()
    dist(x,y) = abs(x-y)
    for i in input
        broke = false
        for c in clusters
            if dist(i,first(c)) < tol
                push!(c,i)
                broke = true
                break
            end
        end
        if !broke
            push!(clusters,[i])
        end
    end
    clusters
end

function recurrences(params, n_last=100)
    y, _, _ = simulate(params)
    skip = round(Int, 1 / params["dt"])
    recs = y[:, (end-n_last*skip):skip:end]
    naive_cluster(recs)
end


function fr2recs(params, fs, rs)
    fr = [(f,r) for f in fs, r in rs]
    Folds.map(fr) do locfr
        myparams = deepcopy(params)
        myparams["f"] = locfr[1]
        myparams["r1"] = locfr[2]
        recurrences(myparams)
    end
end
