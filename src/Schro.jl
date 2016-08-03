module Schro
    import Base.eigs, Base.linspace, Base.spdiagm
    export schro1d

    function schro1d(dx, Nstates, potential)
        xpts = length(potential)
        offdiag = -1 / (2 * dx^2) + zeros(xpts - 1)
        diag = 1 / dx^2 + potential
        H = spdiagm((offdiag, diag, offdiag), (-1, 0, 1))
        evals, evecs = eigs(H, nev=Nstates, which=:SM)
        norms = 1 ./ sqrt(sum(dx * evecs.^2, (1)))
        println("evecs shape=", size(evecs))
        evecs .*= norms
        #print("norms=", sum(dx * evecs.^2, (1)))
        return evals, evecs
    end
end
