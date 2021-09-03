function YoungsMod!(P, Info)
    x = P[:ε][Info[:IE][1]:Info[:IE][2]]
    y = P[:σ_MPa_j][Info[:IE][1]:Info[:IE][2]]
    (E1, ~) = linfit(x, y./1e3)
    x = P[:ε][Info[:IE][3]:Info[:IE][4]]
    y = P[:σ_MPa_j][Info[:IE][3]:Info[:IE][4]]
    (E2, ~) = linfit(x, y./1e3)
    P[:E] = [E1,E2]
end

function Hmod!(P, Info, εrange, Δε)
    x = P[:ε][Info[:I][1]:Info[:I][2]]
    y = P[:σ_MPa_j][Info[:I][1]:Info[:I][2]]./1e3
    ε =  εrange[1]:Δε:εrange[2]
    I = zeros(length(ε))
    for i = 1:length(ε)
        I[n] = findlast(x .< ε[i])
    end
    H = zeros(length(ε))
    for i = 1:length(ε)-1
        (out, ~) = linfit(x[I[i]:I[i+1]],y[I[i]:I[i+1]])
        H[i] = out
    end
    (H[end], ~) = linfit(x[I[1]:I[end]],y[I[1]:I[end]])
    P[:H]   = H
    P[:ε_H] = εrange[1]+Δε/2:Δε:εrange[2]-Δε/2
end
