using HDF5
using Plots
using Schro

function sho(Li=-5e0, Lf=5e0, xpts=101)
    xpts= 201
    Nstates = 3
    x = linspace(Li, Lf, xpts)
    potential = 5e-1 * x.^2
    dx = x[2] - x[1]
    evals, evecs = schro1d(dx, Nstates, potential)
    println("evals=", evals)
    plot(x, evecs)
    # Create an HDF5 file, a group called N3_X201
    h5open("sho.h5", "w") do file
        g = g_create(file, "N$(Nstates)_X$(xpts)")
        # Save the data from this run
        g["x"] = Array(x)
        g["evals"] = evals
        g["evecs"] = evecs
        attrs(g)["Description"] = "$(Nstates) states, $(xpts) points"
    end
end

sho()
