"""
function YoungsMod!(P, Info)
Obtain an estimate of youngs modulus using linear fitting for experiments conducted in the Murrell gas apparatus, UCL, UK.
Inputs:
    * P: a dictionary containing processed mechanical data from the Murrell
    * Info: a dictionary containing information about experimental conditions, including four indices P[:IE] = [1, 2, 3, 4], where 1 is a manual pick start of elastic loading and 2 the end of elastic loading, 3 the start of elastic unloading and 4 the end of unloading.
Stores output of fitting in input dictionary P as P[:E] = [E1,E2], where E1 = loading modulus and E2 unloading modulus in GPa
"""

function YoungsMod!(P, Info)
    x = P[:ε][Info[:IE][1]:Info[:IE][2]]
    y = P[:σ_MPa_j][Info[:IE][1]:Info[:IE][2]]
    (E1, ~) = linfit(x, y./1e3)
    x = P[:ε][Info[:IE][3]:Info[:IE][4]]
    y = P[:σ_MPa_j][Info[:IE][3]:Info[:IE][4]]
    (E2, ~) = linfit(x, y./1e3)
    P[:E] = [E1,E2]
end

"""
function Hmod!(P, Info, εrange, Δε)
Obtain an estimate of tangent moduli using linear fitting for experiments conducted in the Murrell gas apparatus, UCL, UK.
Inputs:
    * P: a dictionary containing processed mechanical data from the Murrell
    * Info: a dictionary containing information about experimental conditions, including two indices P[:I] = [1, 2], where 1 is the hit point and 2 the start of unloading.
Stores output of fitting in input dictionary P as P[:H] the tangent modulus in GPa, and corresponding strain interval as P[:ε_H].
"""

function Hmod!(P, Info, εrange, Δε)
    x = P[:ε][Info[:I][1]:Info[:I][2]]
    y = P[:σ_MPa_j][Info[:I][1]:Info[:I][2]]./1e3
    ε =  εrange[1]:Δε:εrange[2]
    I = zeros(Int64,length(ε))
    for i = 1:length(ε)
        I[i] = findlast(x .< ε[i])
    end
    H = zeros(length(ε)-1)
    for i = 1:length(ε)-1
        (out, ~) = linfit(x[I[i]:I[i+1]],y[I[i]:I[i+1]])
        H[i] = out
    end
    # (H[end], ~) = linfit(x[I[1]:I[end]],y[I[1]:I[end]])
    P[:H]   = H
    P[:ε_H] = εrange[1]+Δε/2:Δε:εrange[2]-Δε/2
end

"""
function YieldStress!(P, Info, Δε)
Obtain an estimate of tangent moduli using linear fitting for experiments conducted in the Murrell gas apparatus, UCL, UK.
Inputs:
    * P: a dictionary containing processed mechanical data from the Murrell
    * Info: a dictionary containing information about experimental conditions, including two indices P[:I] = [1, 2], where 1 is the hit point and 2 the start of unloading.
    * Δε: the deviation to use to define yieldstress
"""

function YieldStress!(P, Info, Δε)
    εpl = P[:ε][Info[:I][1]:Info[:I][2]] .-P[:σ_MPa_j]./(P[:E][1]*1e3)
    σ = P[:σ_MPa_j][Info[:I][1]:Info[:I][2]]
    εpl .-= εpl[Info[:IE][1]]
    P[:σy] = σ[findfirst(εpl .> Δε)]
end
