using FdeSolver
using Plots

# ---------------------------- inputs ----------------------------
tSpan = [0, 130]          # run long enough to see the collapse
h      = 0.01
μ      = [1, 1]
X0     = [0.2, 0.0]       # [N0, N1] – start with no resistance

# ------------------------- parameters ---------------------------
par = [ 0.2,   # r0  : growth rate of N0 (sensitive)
        0.2,   # r1  : baseline growth rate of N1 (resistant)
        0.3,   # k   : kill rate on N0 during pulse
        5.0,   # pulse length
        10.0,  # period between pulses
        0.03 ] # c   : cost added to N1 growth per pulse     ← NEW

# memory‑dependent conversion rate
alpha_memory(T) = 0.01*T           # α(T)

function F(t, x, par)
    N0, N1              = x
    r0, r1, k, pulse,
    period,  c          = par      # c extracted from the vector
    T       = Int(floor(t/period)) # memory = number of past pulses
    α       = alpha_memory(T)

    # --------- N0 dynamics (unchanged) ---------
    if mod(t, period) < pulse
        dN0 = r0*N0*(1-(N0+N1)) - k*N0 - α*N0
    else
        dN0 = r0*N0*(1-(N0+N1))        - α*N0
    end

    # --------- N1 dynamics with metabolic cost ---------
    r1_eff = r1 - c*T                  # cost line
    dN1    = r1_eff*N1*(1-(N0+N1)) + α*N0

    return [dN0, dN1]
end

# -------------------------- solve -------------------------------
t, x = FDEsolver(F, tSpan, X0, μ, par, h = h)

# ------------------------- plotting ----------------------------
p=plot(    xlim=(0, tSpan[2]) )

gr()  # set the GR backend
# 1) Shade intervals [0,5], [10,15], [20,25], ... [90,95] in one color (e.g., lightblue)
for i in 0:22  # i=0 => [0,5], i=1 => [10,15], ... i=9 => [90,95]
    vspan!(p, [i*10, i*10+5], fillalpha=0.1, color=:darkorange2, label="")
end
plot!(t, x[:,1] .+ x[:,2],  lw=4,  color=:gray25,
     label="Total",   legend=:topleft,
     xlabel="Time", ylabel="Abundance")
plot!(t, x[:,1], lw=2, color=:royalblue1,  label="S (sensitive)")
plot!(t, x[:,2], lw=2, color=:darkorange2, label="R (resistant)")

savefig("Memory_rise_and_collapse.png")
