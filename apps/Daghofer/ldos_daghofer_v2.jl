using ArgParse
using JLD

include("src/Tightbinding.jl")
include("src/DaghoferModel.jl")

include("src/PairingModel.jl")
include("src/NaivePairingModel.jl")

include("daghofer.jl")

function parse_commandline()
  argsetting = ArgParseSettings("Daghofer One Dimension",
                                version="Version 0.1",
                                add_version=true)

  @add_arg_table argsetting begin
    "--nx"
      arg_type = Int
      required = true
    "--ny"
      arg_type = Int
      required = true
    "--xi0"
      arg_type = Float64
      nargs = '+'
    "--xi1"
      arg_type = Float64
      nargs = '+'
    "--pairing_d"
      arg_type = Float64
      nargs = '+'
    "--pairing_s"
      arg_type = Float64
      nargs = '+'
    "--pairing_p"
      arg_type = Float64
      nargs = '+'
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

  nx = args["nx"]
  ny = args["ny"]

  Δd = args["pairing_d"]
  Δs = args["pairing_s"]
  Δp = args["pairing_p"]
  ξ₀ = args["xi0"]
  ξ₁ = args["xi1"]
    
  r = compute_mixedspace(nx, ny, Δd, Δs, Δp, ξ₀, ξ₁)
  n = length(r["eigenvalues"])
  ψ = reshape(r["eigenvectors"], (n_nambu, n_orbital, n_basis, nx, n))
  ρ = reshape( sum( (abs(ψ[1,:,:,:, :]).^2), (1, 2)), (nx, n))

  save(args["outfile"],
       "VERSION", 2.0,
       "TYPE": "LDOS",
       "HILBERT/DIMENSION ORDERING" => ["NAMBU", "ORBITAL", "BASIS", "NX"],
       "HILBERT/NUMBER OF NAMBU" => 2,
       "HILBERT/NUMBER OF ORBITALS" => 3,
       "HILBERT/NUMBER OF BASIS SITES" => 2,
       "LATTICE/NX" => nx,
       "LATTICE/NY" => ny,
       "HAMILTONIAN/CORRELATION LENGTH OF S-WAVE" => args["xi0"],
       "HAMILTONIAN/CORRELATION LENGTH OF P-WAVE" => args["xi1"],
       "HAMILTONIAN/PAIRING GAP D-WAVE" => args["pairing_d"],
       "HAMILTONIAN/PAIRING GAP S-WAVE" => args["pairing_s"],
       "HAMILTONIAN/PAIRING GAP P-WAVE" => args["pairing_p"],
       "SPECTRUM/MOMENTUMS", r["momentums"],
       "SPECTRUM/EIGENVALUES", r["eigenvalues"],
       "SPECTRUM/EIGENVECTORS", r["eigenvectors"],
       "SPECTRUM/SPECTRAL WEIGHTS", ρ,
       )
       
end

main()

