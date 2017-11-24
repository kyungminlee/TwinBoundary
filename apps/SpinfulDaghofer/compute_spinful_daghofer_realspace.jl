#
# Compute Real Space Spectrum With an Edge
#
#

using ProgressMeter
using JLD
using PyCall
using PyPlot

@pyimport numpy.linalg as npl

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("spinful_daghofer.jl")

if true
  # no param confirmed
  const λ = 0.0 # confirmed
  const Δd = 0.05
  const Δs = 0.025
  const Δp = 1.0
  const ξ₀ = 1E-8
  const ξ₁ = 1E-8

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const nx = 32
  const ny = 32

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  result_realspace = compute_realspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁, periodic=false)
  eigenvalues_realspace = result_realspace["eigenvalues"]
  sort!(eigenvalues_realspace)

  plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
  save("spinful_daghofer_realspace.jld", result_realspace)
  savefig("foo.png", dpi=300)
end
