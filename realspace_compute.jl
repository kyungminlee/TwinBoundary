#
# Compute Real Space Spectrum With an Edge
#
#

using ProgressMeter
using PyCall
using PyPlot
using ArgParse
using JLD

@pyimport numpy.linalg as npl

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("spinful_daghofer.jl")


function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "OutputFile"
      arg_type = String
      required = true
    "--IsingSpinOrbitCoupling"
      arg_type = Float64
      required = true
    "--PairingDwave"
      arg_type = Float64
      required = true
    "--PairingSwave"
      arg_type = Float64
      required = true
    "--PairingPwave"
      arg_type = Float64
      required = true
    "--CorrelationLengthDwave"
      arg_type = Float64
      required = true
    "--CorrelationLengthSwave"
      arg_type = Float64
      required = true
    "--Nx"
      arg_type = Int
      required = true
    "--Ny"
      arg_type = Int
      required = true
  end
  return parse_args(s)

end


function main()

  args = parse_commandline()
  @show args
  const output_filename = args["OutputFile"]
  const λ = args["IsingSpinOrbitCoupling"]
  const Δd = args["PairingDwave"]
  const Δs = args["PairingSwave"]
  const Δp = args["PairingPwave"]
  const ξ₀ = args["CorrelationLengthDwave"]
  const ξ₁ = args["CorrelationLengthSwave"]
  const nx = args["Nx"]
  const ny = args["Ny"]

  # no param confirmed
  #const λ = 0.1
  #const Δd = 0.3
  #const Δs = 0.1
  #const Δp = 1.0
  #const ξ₀ = 1E-8
  #const ξ₁ = 1E-8

  #const λ = 0.1
  #const Δd = 0.3
  #const Δs = 0.1
  #const ξ₀ = 1E-8
  #const ξ₁ = 1E-8


  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  const UP = 1
  const DN = 2

  const kxs = linspace(0, 2π, nx+1)[1:end-1]
  const kys = linspace(0, 2π, ny+1)[1:end-1]

  result_realspace = compute_realspace_sparse(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁; periodic=false)
  eigenvalues = result_realspace["eigenvalues"]
  eigenvectors = result_realspace["eigenvectors"]

  n_normal = n_basis * n_orbital * n_spin * n_spin
  n_reduced = length(eigenvalues)
  n_full, n_reduced2 = size(eigenvectors)

  @assert n_reduced == n_reduced2

  eigenvectors_view = reshape(eigenvectors, (n_nambu, n_orbital, n_spin, n_basis, nx, ny, n_reduced))
  ρ = abs2.(eigenvectors_view[1, :, :, :, :, :, :])

  save(output_filename,
    "IsingSpinOrbitCoupling", λ,
    "PairingDwave", Δd,
    "PairingSwave", Δs,
    "PairingPwave", Δp,
    "CorrelationLength0", ξ₀,
    "CorrelationLength1", ξ₁,
    "SystemSize", (nx, ny),
    "Eigenvalues", eigenvalues,
    "Densities", ρ)
  #plot(eigenvalues_realspace .+ 0, ".-", alpha=0.5, label="realspace", linewidth=2)
end

main()
