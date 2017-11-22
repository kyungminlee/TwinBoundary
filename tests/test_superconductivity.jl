using Base.Test

include("Tightbinding.jl")
using Tightbinding

Base.ishermitian{T <: Number}(m ::Matrix{T}; atol=0.0) = isapprox(m, m'; atol=atol)

tb = TightbindingModel{2, Complex128}(
  [("orb1", [0.0, 0.0]),
   ("orb2", [0.0, 0.0])])

t = 1.0
#tb = TightbindingModel{2, Complex128}([("s", [0.0, 0.0])])
const n_nambu = 2
const n_basis = 2

#add_hopping!(tb, ( 1, 0), "s", "s", t)
#add_hopping!(tb, ( 0, 1), "s", "s", t)

#add_hopping!(tb, (-1, 0), "s", "s", conj(t))
#add_hopping!(tb, ( 0,-1), "s", "s", conj(t))

tb2 = double_unitcell(tb)

@show tb
@show tb2
#=

nx, ny = 4, 4
kxs = linspace(0, 2*pi, nx+1)[1:end-1]
kys = linspace(0, 2*pi, ny+1)[1:end-1]

let
  hamiltonian = zeros(Complex128, (n_nambu, n_orbital, n_basis, nx, ny, n_nambu, n_orbital, n_basis, nx, ny))
  t = double_daghofer_realspace_dense
  hamiltonian[1, :, :, :, 1, :, :, :] = reshape( t, (n_basis, nx, ny, n_basis, nx, ny))
  hamiltonian[2, :, :, :, 2, :, :, :] = reshape(-t, (n_basis, nx, ny, n_basis, nx, ny))
end

=#