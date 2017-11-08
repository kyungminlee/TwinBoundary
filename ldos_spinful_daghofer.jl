using ArgParse
using JLD
using HDF5

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("spinful_daghofer.jl")

function parse_commandline()
  argsetting = ArgParseSettings("Spinful Daghofer One Dimension",
                                version="Version 0.1",
                                add_version=true)

  @add_arg_table argsetting begin
    "--nx"
      arg_type = Int
      required = true
    "--ny"
      arg_type = Int
      required = true
    "--soc"
      arg_type = Float64
      required = true
    "--xi0"
      arg_type = Float64
      required = true
    "--xi1"
      arg_type = Float64
      required = true
    "--pairing_d"
      arg_type = Float64
      required = true
    "--pairing_s"
      arg_type = Float64
      required = true
    "--pairing_p"
      arg_type = Float64
      required = true
    "--outfile", "-o"
      arg_type = String
      required = true
      help = "output filename"
  end
  parsed_args = parse_args(ARGS, argsetting)
  return parsed_args
end


function main()
  args = parse_commandline()

  const n_nambu = 2
  const n_basis = 2
  const n_orbital = 3
  const n_spin = 2

  nx = args["nx"]
  ny = args["ny"]
  λ = args["soc"]
  Δd = args["pairing_d"]
  Δs = args["pairing_s"]
  Δp = args["pairing_p"]
  ξ₀ = args["xi0"]
  ξ₁ = args["xi1"]

  r = compute_mixedspace(nx, ny, λ, Δd, Δs, Δp, ξ₀, ξ₁)
  n = length(r["eigenvalues"])
  ψ = reshape(r["eigenvectors"], (n_nambu, n_orbital, n_spin, n_basis, nx, n))
  ρ = reshape( sum( (abs.(ψ[1,:,:,:,:,:]).^2), (1, 2, 3)), (nx, n))


  # JLD
  #=
  save(args["outfile"],
       "VERSION", 2.0,
       "TYPE", "LDOS",
       "HILBERT/DIMENSION ORDERING", ["NAMBU", "ORBITAL", "SPIN", "BASIS", "NX"],
       "HILBERT/NUMBER OF NAMBU", 2,
       "HILBERT/NUMBER OF ORBITALS", 3,
       "HILBERT/NUMBER OF SPINS", 2,
       "HILBERT/NUMBER OF BASIS SITES", 2,
       "LATTICE/NX", nx,
       "LATTICE/NY", ny,
       "HAMILTONIAN/CORRELATION LENGTH OF S-WAVE", args["xi0"],
       "HAMILTONIAN/CORRELATION LENGTH OF P-WAVE", args["xi1"],
       "HAMILTONIAN/PAIRING GAP D-WAVE", args["pairing_d"],
       "HAMILTONIAN/PAIRING GAP S-WAVE", args["pairing_s"],
       "HAMILTONIAN/PAIRING GAP P-WAVE", args["pairing_p"],
       "SPECTRUM/MOMENTUMS", r["momentums"],
       "SPECTRUM/EIGENVALUES", r["eigenvalues"],
       "SPECTRUM/EIGENVECTORS", r["eigenvectors"],
       "SPECTRUM/SPECTRAL WEIGHTS", ρ,
       )
  =#
  h5open(args["outfile"], "w") do file
    file["VERSION"] = "2.0"
    file["TYPE"] = "LDOS"
    file["HILBERT/DIMENSION ORDERING"] = ["NAMBU", "ORBITAL", "SPIN", "BASIS", "NX"]
    file["HILBERT/NUMBER OF NAMBU"] = n_nambu
    file["HILBERT/NUMBER OF ORBITALS"] = n_orbital
    file["HILBERT/NUMBER OF SPINS"] = n_spin
    file["HILBERT/NUMBER OF BASIS SITES"] = n_basis
    file["LATTICE/NX"] =  nx
    file["LATTICE/NY"] =  ny
    file["HAMILTONIAN/CORRELATION LENGTH OF S-WAVE"] = args["xi0"]
    file["HAMILTONIAN/CORRELATION LENGTH OF P-WAVE"] = args["xi1"]
    file["HAMILTONIAN/PAIRING GAP D-WAVE"] = args["pairing_d"]
    file["HAMILTONIAN/PAIRING GAP S-WAVE"] = args["pairing_s"]
    file["HAMILTONIAN/PAIRING GAP P-WAVE"] = args["pairing_p"]
    file["HAMILTONIAN/SPIN ORBIT COUPLING"] = args["soc"]
    file["SPECTRUM/MOMENTUMS"] = r["momentums"]
    file["SPECTRUM/EIGENVALUES"] = r["eigenvalues"]
    file["SPECTRUM/EIGENVECTORS REAL"] = real( r["eigenvectors"] )
    file["SPECTRUM/EIGENVECTORS IMAG"] = imag( r["eigenvectors"] )
    file["SPECTRUM/SPECTRAL WEIGHTS"] = ρ
  end
end

main()
