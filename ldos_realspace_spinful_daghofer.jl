using ArgParse
using JLD
using HDF5

version="0.1"

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("spinful_daghofer.jl")

function parse_commandline()
  argsetting = ArgParseSettings("Spinful Daghofer One Dimension",
                                version=version,
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
  ρ = reshape( sum( (abs.(ψ[1,:,:,:,:,:]).^2), (1, 2)), (n_basis, nx, n))


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
  n_eigen = n_nambu * n_orbital * n_spin * n_basis * nx
  n_total = ny * n_eigen
  h5open(args["outfile"], "w") do file
    file["VERSION"] = version
    file["TYPE"] = "LDOS"

    g = g_create(file, "HILBERT")
    g["DIMENSION ORDERING"] = ["NAMBU", "ORBITAL", "SPIN", "BASIS", "NX"]
    g["NUMBER OF NAMBU"] = n_nambu
    g["NUMBER OF ORBITALS"] = n_orbital
    g["NUMBER OF SPINS"] = n_spin
    g["NUMBER OF BASIS SITES"] = n_basis

    g = g_create(file, "LATTICE")
    g["NX"] =  nx
    g["NY"] =  ny

    g = g_create(file, "HAMILTONIAN")
    g["CORRELATION LENGTH OF S-WAVE"] = args["xi0"]
    g["CORRELATION LENGTH OF P-WAVE"] = args["xi1"]
    g["PAIRING GAP D-WAVE"] = args["pairing_d"]
    g["PAIRING GAP S-WAVE"] = args["pairing_s"]
    g["PAIRING GAP P-WAVE"] = args["pairing_p"]
    g["SPIN ORBIT COUPLING"] = args["soc"]

    g = g_create(file, "SPECTRUM")
    g["MOMENTUMS",         "chunk", (n_eigen,), "shuffle", (), "compress", 9] = r["momentums"]
    g["EIGENVALUES",       "chunk", (n_eigen,), "shuffle", (), "compress", 9] = r["eigenvalues"]
    #g["EIGENVECTORS REAL", "chunk", (n_eigen, n_eigen), "shuffle", (), "compress", 9] = real( r["eigenvectors"] )
    #g["EIGENVECTORS IMAG", "chunk", (n_eigen, n_eigen), "shuffle", (), "compress", 9] = imag( r["eigenvectors"] )
    g["SPECTRAL WEIGHTS",  "chunk", (n_basis, nx, n_eigen), "shuffle", (), "compress", 9] = ρ
  end
end

main()
