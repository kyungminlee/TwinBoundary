using Base.Test

include("Tightbinding.jl")
using Tightbinding

Base.ishermitian{T <: Number}(m ::Matrix{T}; atol=0.0) = isapprox(m, m'; atol=atol)

t = 1.0# + 0.5im
tb = TightbindingModel{2, Complex128}([("s", [0.0, 0.0])])

add_hopping!(tb, ( 1, 0), "s", "s", t)
#add_hopping!(tb, ( 0, 1), "s", "s", t)

@test !ishermitian(tb)

#add_hopping!(tb, (-1, 0), "s", "s", conj(t))
#add_hopping!(tb, ( 0,-1), "s", "s", conj(t))

#@test ishermitian(tb)

tb2 = double_unitcell(tb)
@show tb
#dump(tb)
@show tb2
#dump(tb2)